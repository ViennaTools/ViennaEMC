#ifndef EMC_NUMERIC_BAND_STRUCTURE_HPP
#define EMC_NUMERIC_BAND_STRUCTURE_HPP

#include <array>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <limits>
#include <algorithm>
#include <stdexcept>
#include <type_traits>
#include <string>
#include <vector>

#include <hdf5.h>

/*! \brief Full-band numeric band structure loaded from a .matpkg material
 * package (see FullBandMC spec/matpkg-spec.md, v0.1).
 *
 * Represents E_n(k) on a tetrahedral mesh of the full Brillouin zone in
 * crystal (fractional) coordinates. Provides energy and group velocity by
 * barycentric interpolation inside tetrahedra, located via an adjacency-walk
 * (the package is REQUIRED to ship precomputed adjacency; this class never
 * derives it - spec rule R2). The package is opened read-only and never
 * written to (spec rule R1).
 *
 * Wave vectors at the public interface are Cartesian [1/m]; energies are [eV];
 * velocities are [m/s], matching ViennaEMC conventions.
 */
template <class T> class emcNumericBandStructure {
  static_assert(std::is_floating_point<T>::value, "T must be floating point");

public:
  using Vec3 = std::array<T, 3>;
  using SizeTypeOrInt64 = std::size_t;

  /// opens and loads the package; throws std::runtime_error on any
  /// missing dataset (hard error by design - no silent fallbacks)
  explicit emcNumericBandStructure(const std::string &fileName) {
    hid_t file = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0)
      throw std::runtime_error("matpkg: cannot open " + fileName);

    readMatrix3(file, "/lattice/reciprocal_basis", recBasis);
    invert3x3(recBasis, recBasisInv);

    read2D(file, "/kmesh/points", points, 3);
    read2D(file, "/kmesh/tetrahedra", tets, 4);
    read2D(file, "/kmesh/adjacency", adjacency, 4);
    if (tets.size() != adjacency.size())
      throw std::runtime_error("matpkg: tetrahedra/adjacency size mismatch");

    readBands(file);
    readDOS(file);
    readIbzMap(file);
    H5Fclose(file);

    precomputeTetGeometry();
    buildBuckets();
  }

  emcNumericBandStructure(const emcNumericBandStructure &) = delete;
  emcNumericBandStructure &operator=(const emcNumericBandStructure &) = delete;

  SizeTypeOrInt64 getNrPoints() const { return points.size(); }
  SizeTypeOrInt64 getNrTetrahedra() const { return tets.size(); }
  SizeTypeOrInt64 getNrBands() const { return energies.size(); }
  T getBandMinimum(std::size_t band = 0) const { return bandMin[band]; }

  /// returns energy [eV] at Cartesian wave-vector k [1/m] for given band;
  /// tetHint caches the last containing tetrahedron of the caller (particle)
  /// and is updated - pass the same variable across a particle's lifetime.
  T getEnergy(const Vec3 &kCart, std::size_t band, std::int64_t &tetHint) const {
    Vec3 lam;
    const std::int64_t t = locate(toFrac(kCart), tetHint, lam);
    return tetEnergy(t, band, lam);
  }

  /// energy inside tet t at barycentric lam (lam0 = 1 - sum).
  /// With per-vertex velocities available, adds the quadratic Hermite term
  ///     -1/2 sum_{i<j} lam_i lam_j (grad_j - grad_i).(x_j - x_i)
  /// which is EXACT for a quadratic band (verified in 1D: reproduces x^2)
  /// and raises the interpolation order from O(h^2) to O(h^3). It matters
  /// because linear interpolation of a convex band overestimates E inside
  /// every tetrahedron: the equilibrium <E> came out 13.4 meV too hot at
  /// 32^3 and 2.3 meV at 64^3, which with steeply energy-dependent rates
  /// drives both the residual mesh dependence and the low-field cubic
  /// anisotropy. See validation/m1-si-production/engine_transport_audit.md
  T tetEnergy(std::int64_t t, std::size_t band, const Vec3 &lam) const {
    const auto &tet = tets[t];
    const auto &E = energies[band];
    const T e0 = E[tet[0]];
    T e = e0 + lam[0] * (E[tet[1]] - e0) + lam[1] * (E[tet[2]] - e0) +
          lam[2] * (E[tet[3]] - e0);
    if (!quadEnergy || band >= vertexVel.size() || vertexVel[band].empty())
      return e;
    const T l[4] = {1 - lam[0] - lam[1] - lam[2], lam[0], lam[1], lam[2]};
    const auto &V = vertexVel[band];
    for (int i = 0; i < 4; i++)
      for (int j = i + 1; j < 4; j++) {
        const T w = l[i] * l[j];
        if (w == 0)
          continue;
        // grad_cart E = hbar * v / qe  [eV m];  x = B^T f  [1/m]
        Vec3 df{points[tet[j]][0] - points[tet[i]][0],
                points[tet[j]][1] - points[tet[i]][1],
                points[tet[j]][2] - points[tet[i]][2]};
        for (int c = 0; c < 3; c++)   // shortest image (periodic mesh)
          df[c] -= std::nearbyint(df[c]);
        const Vec3 dx = fracToCart(df);
        T dot = 0;
        for (int c = 0; c < 3; c++)
          dot += (V[tet[j]][c] - V[tet[i]][c]) * dx[c];
        e -= static_cast<T>(0.5) * w * dot * hbar / qe;
      }
    return e;
  }

  /// Parameter t in [0,1] where the band reaches `energy` along the edge
  /// va->vb of tet `t`, CONSISTENT with tetEnergy(). With the quadratic term
  /// the edge profile is
  ///   E(t) = Ea + t (Eb - Ea) - (c/2) t (1-t),  c = (grad_b - grad_a).(x_b - x_a)
  /// i.e. a quadratic; solving it (instead of the linear crossing) keeps the
  /// final-state sampler on the SAME energy surface the engine reads. Without
  /// that, every scattering event drops the carrier a little (linear over-
  /// estimates a convex band) - a fake cooling channel that collapsed the
  /// elastic test package to the band bottom.
  /// Returns false if the edge does not cross `energy`.
  bool edgeCrossing(std::int64_t tIdx, std::size_t band, int ia, int ib,
                    T energy, T &tOut) const {
    const auto &tet = tets[tIdx];
    const T ea = energies[band][tet[ia]], eb = energies[band][tet[ib]];
    T c = 0;
    if (quadEnergy && band < vertexVel.size() && !vertexVel[band].empty()) {
      const auto &V = vertexVel[band];
      Vec3 df{points[tet[ib]][0] - points[tet[ia]][0],
              points[tet[ib]][1] - points[tet[ia]][1],
              points[tet[ib]][2] - points[tet[ia]][2]};
      for (int d = 0; d < 3; d++)
        df[d] -= std::nearbyint(df[d]);
      const Vec3 dx = fracToCart(df);
      for (int d = 0; d < 3; d++)
        c += (V[tet[ib]][d] - V[tet[ia]][d]) * dx[d];
      c *= hbar / qe;                       // eV
    }
    if (c == 0) {                           // linear edge
      if ((ea - energy) * (eb - energy) >= 0)
        return false;
      tOut = (energy - ea) / (eb - ea);
      return true;
    }
    // (c/2) t^2 + (eb - ea - c/2) t + (ea - energy) = 0
    const T A = c / 2, B = eb - ea - c / 2, C = ea - energy;
    const T disc = B * B - 4 * A * C;
    if (disc < 0)
      return false;
    const T sq = std::sqrt(disc);
    const T r1 = (-B + sq) / (2 * A), r2 = (-B - sq) / (2 * A);
    // pick the root inside the edge; if both, the one nearer the linear guess
    const bool ok1 = r1 >= 0 && r1 <= 1, ok2 = r2 >= 0 && r2 <= 1;
    if (!ok1 && !ok2)
      return false;
    if (ok1 && ok2) {
      const T lin = (eb != ea) ? (energy - ea) / (eb - ea) : T(0.5);
      tOut = (std::abs(r1 - lin) <= std::abs(r2 - lin)) ? r1 : r2;
    } else {
      tOut = ok1 ? r1 : r2;
    }
    return true;
  }

  /// quadratic energy interpolation on/off (default on when velocities exist)
  void setQuadraticEnergy(bool on) { quadEnergy = on; }
  bool hasQuadraticEnergy() const {
    return quadEnergy && !vertexVel.empty();
  }

  /// returns group velocity [m/s] = grad_k E / hbar (constant per tetrahedron)
  Vec3 getVelocity(const Vec3 &kCart, std::size_t band,
                   std::int64_t &tetHint) const {
    Vec3 lam;
    const std::int64_t t = locate(toFrac(kCart), tetHint, lam);
    if (band < vertexVel.size() && !vertexVel[band].empty()) {
      // continuous velocity field: barycentric blend of vertex velocities
      const auto &tet = tets[t];
      const auto &V = vertexVel[band];
      const T l0 = 1 - lam[0] - lam[1] - lam[2];
      Vec3 v{};
      for (int c = 0; c < 3; c++)
        v[c] = l0 * V[tet[0]][c] + lam[0] * V[tet[1]][c] +
               lam[1] * V[tet[2]][c] + lam[2] * V[tet[3]][c];
      return v;
    }
    return tetVelocity(t, band);
  }

  bool hasVertexVelocities() const { return !vertexVel.empty(); }

  /// interpolated total density of states [1/(eV m^3)] at energy [eV]
  T getDOS(T energy) const {
    if (energy <= dosGrid.front() || energy >= dosGrid.back())
      return 0;
    const auto it =
        std::upper_bound(dosGrid.begin(), dosGrid.end(), energy) - 1;
    const std::size_t i = it - dosGrid.begin();
    const T f = (energy - dosGrid[i]) / (dosGrid[i + 1] - dosGrid[i]);
    return dosTotal[i] * (1 - f) + dosTotal[i + 1] * f;
  }

  // ---- mesh access for full-band scattering (final-state sampling) ----
  const std::vector<std::array<T, 3>> &getPoints() const { return points; }
  const std::vector<std::array<std::int64_t, 4>> &getTetrahedra() const {
    return tets;
  }
  const std::vector<T> &getBandEnergies(std::size_t band) const {
    return energies[band];
  }
  const std::array<std::array<T, 3>, 3> &getReciprocalBasis() const {
    return recBasis;
  }
  /// fractional -> Cartesian [1/m]  (k = B^T f)
  Vec3 fracToCart(const Vec3 &f) const {
    Vec3 k{};
    for (int i = 0; i < 3; i++)
      k[i] = recBasis[0][i] * f[0] + recBasis[1][i] * f[1] +
             recBasis[2][i] * f[2];
    return k;
  }

  /// public point location: containing tetrahedron + barycentric coords
  /// (lam1..3 of vertices 1..3; vertex-0 weight = 1 - sum). For consumers
  /// that interpolate their own per-vertex data (e.g. k-resolved rates).
  std::int64_t locateTet(const Vec3 &kCart, std::int64_t &tetHint,
                         Vec3 &lam) const {
    return locate(toFrac(kCart), tetHint, lam);
  }

  /// IBZ map (spec /kmesh/ibz/map_full_to_ibz); empty if package lacks it
  const std::vector<std::int64_t> &getIbzMap() const { return ibzMap; }
  /// per-point symmetry-op index (R . k_ibz = k_full) and the integer k-ops
  const std::vector<std::int64_t> &getIbzOpMap() const { return ibzOpMap; }
  const std::vector<std::array<std::array<T, 3>, 3>> &getSymOps() const {
    return symOps;
  }

  /// nr of walk-fallback full scans (diagnostic; large values indicate
  /// particles teleporting, e.g. after BZ wrap - expected occasionally)
  std::size_t getNrFallbackScans() const { return fallbackScans.load(); }
  /// most negative barycentric coordinate ever accepted via shell-snapping
  T getWorstShellDeficit() const { return worstShellDeficit.load(); }

private:
  std::array<std::array<T, 3>, 3> recBasis;    // rows b1,b2,b3 [1/m]
  std::array<std::array<T, 3>, 3> recBasisInv; // columns for cart->frac
  std::vector<std::array<T, 3>> points;        // fractional [0,1)
  std::vector<std::array<std::int64_t, 4>> tets;
  std::vector<std::array<std::int64_t, 4>> adjacency;
  std::vector<std::vector<T>> energies; // [band][point], eV
  // OPTIONAL per-vertex band velocities [band][point][xyz], m/s. The
  // per-tetrahedron gradient of the linear interpolant is piecewise CONSTANT,
  // so the field changes a carrier's velocity only when it crosses a tet
  // face. Between collisions the field moves k by a small fraction of a cell,
  // and that staircase inflates the drift response: measured +160% at 1.6
  // cells per thermal momentum, +27% at 12.8, converging only as O(h).
  // Interpolating analytic vertex velocities barycentrically restores a
  // continuous v(k). See validation/m1-si-production/engine_transport_audit.md
  std::vector<std::vector<std::array<T, 3>>> vertexVel;
  // EXPERIMENTAL, OFF by default (enable with setQuadraticEnergy(true)).
  // The quadratic Hermite energy term is correct in itself and fixes the
  // equilibrium <E> (13.4 -> ~0 meV excess at 32^3), but the final-state
  // machinery around it is not yet consistent: the energy BIN WEIGHTS use
  // the Blochl linear-tetrahedron DOS, which is the wrong measure once a
  // tet can hold an interior extremum. Half-consistent, it moves mu the
  // wrong way on the analytic package (676.5 exact: 678 with linear bins,
  // 737 with quadratic-aware ranges). Turning it on for production needs a
  // proper quadratic-tetrahedron DOS + sampling scheme.
  bool quadEnergy = false;
  std::vector<T> bandMin;
  std::vector<T> dosGrid, dosTotal;
  std::vector<std::array<T, 9>> tetInv; // inverse edge matrix per tet
  // bucket grid over the mesh's frac bounding box: O(1) walk re-seeding
  static constexpr int NB = 32;
  std::array<T, 3> bboxLo{}, bboxSpan{};
  std::vector<std::int64_t> buckets; // NB^3, tet seed per bucket (-1 empty)
  // atomic: the const query path may run from several ensemble threads
  mutable std::atomic<std::size_t> fallbackScans{0};
  mutable std::atomic<T> worstShellDeficit{0};

  static constexpr T eps = static_cast<T>(1e-7); // legacy meshes carry float32-era coords: inter-tet cracks ~1e-7
  static constexpr T hbar = static_cast<T>(1.054571817e-34); // J s
  static constexpr T qe = static_cast<T>(1.602176634e-19);   // J / eV

  /// Cartesian [1/m] -> fractional, wrapped to [0,1)
  Vec3 toFrac(const Vec3 &kCart) const {
    Vec3 f{};
    for (int i = 0; i < 3; i++) {
      // k = B^T f  =>  f = (B^T)^-1 k = (B^-1)^T k
      f[i] = recBasisInv[0][i] * kCart[0] + recBasisInv[1][i] * kCart[1] +
             recBasisInv[2][i] * kCart[2];
      f[i] -= std::nearbyint(f[i]); // reduce to [-0.5, 0.5)
    }
    return f;
  }

  /// barycentric coords of p in tet t (lam1..3; lam0 = 1-sum)
  void barycentric(std::int64_t t, const Vec3 &p, Vec3 &lam) const {
    const auto &v0 = points[tets[t][0]];
    const auto &M = tetInv[t];
    const T d0 = p[0] - v0[0], d1 = p[1] - v0[1], d2 = p[2] - v0[2];
    lam[0] = M[0] * d0 + M[1] * d1 + M[2] * d2;
    lam[1] = M[3] * d0 + M[4] * d1 + M[5] * d2;
    lam[2] = M[6] * d0 + M[7] * d1 + M[8] * d2;
  }

  bool inside(const Vec3 &lam) const {
    return lam[0] >= -eps && lam[1] >= -eps && lam[2] >= -eps &&
           (lam[0] + lam[1] + lam[2]) <= 1 + eps;
  }

  std::int64_t bucketOf(const Vec3 &p) const {
    int ix[3];
    for (int i = 0; i < 3; i++) {
      const T r = (p[i] - bboxLo[i]) / bboxSpan[i];
      if (r < 0 || r >= 1)
        return -1;
      ix[i] = static_cast<int>(r * NB);
    }
    return (ix[0] * NB + ix[1]) * NB + ix[2];
  }

  void buildBuckets() {
    for (int i = 0; i < 3; i++) {
      T lo = points[0][i], hi = points[0][i];
      for (const auto &p : points) {
        lo = std::min(lo, p[i]);
        hi = std::max(hi, p[i]);
      }
      bboxLo[i] = lo;
      bboxSpan[i] = (hi - lo) * (1 + static_cast<T>(1e-9));
    }
    buckets.assign(NB * NB * NB, -1);
    for (std::int64_t t = 0; t < (std::int64_t)tets.size(); t++) {
      Vec3 c{};
      for (int v = 0; v < 4; v++)
        for (int i = 0; i < 3; i++)
          c[i] += points[tets[t][v]][i] / 4;
      const std::int64_t bkt = bucketOf(c);
      if (bkt >= 0 && buckets[bkt] < 0)
        buckets[bkt] = t;
    }
    // fill empty buckets from nearest filled neighbor (sweep)
    for (int pass = 0; pass < NB; pass++) {
      bool changed = false;
      for (int x = 0; x < NB; x++)
        for (int y = 0; y < NB; y++)
          for (int z = 0; z < NB; z++) {
            const int i = (x * NB + y) * NB + z;
            if (buckets[i] >= 0)
              continue;
            for (const int d : {i - 1, i + 1, i - NB, i + NB, i - NB * NB,
                                i + NB * NB}) {
              if (d >= 0 && d < NB * NB * NB && buckets[d] >= 0) {
                buckets[i] = buckets[d];
                changed = true;
                break;
              }
            }
          }
      if (!changed)
        break;
    }
  }

  /// bounded adjacency walk from a seed; returns tet or -1
  std::int64_t walkFrom(std::int64_t seed, const Vec3 &q, Vec3 &lam,
                        int maxSteps) const {
    std::int64_t t = seed;
    for (int step = 0; step < maxSteps; step++) {
      barycentric(t, q, lam);
      if (inside(lam))
        return t;
      const T l0 = 1 - lam[0] - lam[1] - lam[2];
      int worst = 0;
      T wv = l0;
      if (lam[0] < wv) { wv = lam[0]; worst = 1; }
      if (lam[1] < wv) { wv = lam[1]; worst = 2; }
      if (lam[2] < wv) { wv = lam[2]; worst = 3; }
      const std::int64_t nx = adjacency[t][worst];
      if (nx < 0)
        return -1;
      t = nx;
    }
    return -1;
  }

  /// adjacency walk from hint; falls back to full scan (and counts it)
  std::int64_t locate(const Vec3 &p, std::int64_t &hint, Vec3 &lam) const {
    std::int64_t t = (hint >= 0 && hint < (std::int64_t)tets.size()) ? hint : 0;
    for (std::size_t step = 0; step < tets.size(); step++) {
      barycentric(t, p, lam);
      if (inside(lam)) {
        hint = t;
        return t;
      }
      // move towards the face whose barycentric coordinate is most negative;
      // local coords: lam0 belongs to vertex1... vertex0 coord = 1-sum
      const T l0 = 1 - lam[0] - lam[1] - lam[2];
      int worst = 0; // local vertex index opposite the exit face
      T wv = l0;
      if (lam[0] < wv) { wv = lam[0]; worst = 1; }
      if (lam[1] < wv) { wv = lam[1]; worst = 2; }
      if (lam[2] < wv) { wv = lam[2]; worst = 3; }
      const std::int64_t next = adjacency[t][worst];
      if (next < 0)
        break; // hit BZ boundary: p wrapped to other side -> scan
      t = next;
    }
    // the mesh's fundamental domain (e.g. Wigner-Seitz BZ) need not contain
    // the round-reduced point: try bucket-seeded walks on all G-images
    for (int gx = -1; gx <= 1; gx++)
      for (int gy = -1; gy <= 1; gy++)
        for (int gz = -1; gz <= 1; gz++) {
          const Vec3 q{p[0] + gx, p[1] + gy, p[2] + gz};
          const std::int64_t bkt = bucketOf(q);
          if (bkt < 0 || buckets[bkt] < 0)
            continue;
          const std::int64_t t2 = walkFrom(buckets[bkt], q, lam, 256);
          if (t2 >= 0) {
            hint = t2;
            return t2;
          }
        }
    // fallback: full scan over all tets and all images (rare). Legacy meshes
    // (float32-era coordinates) leave thin uncovered shells at the BZ boundary;
    // accept the best tetrahedron within a small extrapolation tolerance.
    fallbackScans.fetch_add(1, std::memory_order_relaxed);
    std::int64_t bestTet = -1;
    T bestWorst = -std::numeric_limits<T>::max();
    Vec3 bestLam{};
    for (int gx = -1; gx <= 1; gx++)
      for (int gy = -1; gy <= 1; gy++)
        for (int gz = -1; gz <= 1; gz++) {
          const Vec3 q{p[0] + gx, p[1] + gy, p[2] + gz};
          for (std::int64_t s = 0; s < (std::int64_t)tets.size(); s++) {
            barycentric(s, q, lam);
            const T l0 = 1 - lam[0] - lam[1] - lam[2];
            const T worst = std::min(std::min(lam[0], lam[1]),
                                     std::min(lam[2], l0));
            if (worst > bestWorst) {
              bestWorst = worst;
              bestTet = s;
              bestLam = lam;
            }
            if (worst >= -eps) {
              hint = s;
              return s;
            }
          }
        }
    // snap to nearest tet: legacy meshes have boundary shells of varying
    // width; linear extrapolation error is bounded by |bestWorst| * tet
    // energy range. A large deficit indicates a real bug, not a shell.
    if (bestWorst < static_cast<T>(-0.05))
      throw std::runtime_error(
          "emcNumericBandStructure: point far outside mesh (deficit " +
          std::to_string(static_cast<double>(-bestWorst)) + ")");
    for (T cur = worstShellDeficit.load(std::memory_order_relaxed);
         bestWorst < cur &&
         !worstShellDeficit.compare_exchange_weak(cur, bestWorst,
                                                  std::memory_order_relaxed);)
      ;   // keep the running minimum without a lock
    lam = bestLam;
    hint = bestTet;
    return bestTet;
  }

  Vec3 tetVelocity(std::int64_t t, std::size_t band) const {
    const auto &tet = tets[t];
    const auto &E = energies[band];
    const auto &M = tetInv[t];
    const T e0 = E[tet[0]];
    const T dE1 = E[tet[1]] - e0, dE2 = E[tet[2]] - e0, dE3 = E[tet[3]] - e0;
    // grad_frac E = M^T * dE ; grad_cart E = recBasisInv^T applied per chain
    Vec3 gf{M[0] * dE1 + M[3] * dE2 + M[6] * dE3,
            M[1] * dE1 + M[4] * dE2 + M[7] * dE3,
            M[2] * dE1 + M[5] * dE2 + M[8] * dE3};
    // dE/dk_cart = sum_j (dE/df_j) (df_j/dk_cart) = recBasisInv^T row combo
    Vec3 v{};
    for (int i = 0; i < 3; i++)
      // df_j/dk_i = ((B^-1)^T)_{ji} = (B^-1)_{ij}
      v[i] = (gf[0] * recBasisInv[i][0] + gf[1] * recBasisInv[i][1] +
              gf[2] * recBasisInv[i][2]) *
             qe / hbar; // eV*m -> J*m, /hbar -> m/s
    return v;
  }

  void precomputeTetGeometry() {
    tetInv.resize(tets.size());
    for (std::size_t t = 0; t < tets.size(); t++) {
      const auto &v0 = points[tets[t][0]];
      std::array<std::array<T, 3>, 3> A;
      for (int c = 1; c < 4; c++)
        for (int r = 0; r < 3; r++)
          A[r][c - 1] = points[tets[t][c]][r] - v0[r];
      std::array<std::array<T, 3>, 3> Ai;
      invert3x3(A, Ai);
      for (int r = 0; r < 3; r++)
        for (int c = 0; c < 3; c++)
          tetInv[t][r * 3 + c] = Ai[r][c];
    }
  }

  static void invert3x3(const std::array<std::array<T, 3>, 3> &A,
                        std::array<std::array<T, 3>, 3> &out) {
    const T det = A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) -
                  A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0]) +
                  A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);
    if (std::abs(det) < 1e-300)
      throw std::runtime_error("matpkg: degenerate tetrahedron/basis");
    const T id = 1 / det;
    out[0][0] = (A[1][1] * A[2][2] - A[1][2] * A[2][1]) * id;
    out[0][1] = (A[0][2] * A[2][1] - A[0][1] * A[2][2]) * id;
    out[0][2] = (A[0][1] * A[1][2] - A[0][2] * A[1][1]) * id;
    out[1][0] = (A[1][2] * A[2][0] - A[1][0] * A[2][2]) * id;
    out[1][1] = (A[0][0] * A[2][2] - A[0][2] * A[2][0]) * id;
    out[1][2] = (A[0][2] * A[1][0] - A[0][0] * A[1][2]) * id;
    out[2][0] = (A[1][0] * A[2][1] - A[1][1] * A[2][0]) * id;
    out[2][1] = (A[0][1] * A[2][0] - A[0][0] * A[2][1]) * id;
    out[2][2] = (A[0][0] * A[1][1] - A[0][1] * A[1][0]) * id;
  }

  // ------------------------- HDF5 helpers -------------------------
  static hid_t openDS(hid_t file, const char *path) {
    hid_t ds = H5Dopen2(file, path, H5P_DEFAULT);
    if (ds < 0)
      throw std::runtime_error(std::string("matpkg: missing dataset ") + path);
    return ds;
  }

  void readMatrix3(hid_t file, const char *path,
                   std::array<std::array<T, 3>, 3> &out) {
    hid_t ds = openDS(file, path);
    double buf[9];
    H5Dread(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
    H5Dclose(ds);
    for (int r = 0; r < 3; r++)
      for (int c = 0; c < 3; c++)
        out[r][c] = static_cast<T>(buf[r * 3 + c]);
  }

  template <class U, std::size_t N>
  void read2D(hid_t file, const char *path,
              std::vector<std::array<U, N>> &out, std::size_t expectCols) {
    hid_t ds = openDS(file, path);
    hid_t sp = H5Dget_space(ds);
    hsize_t dims[2];
    H5Sget_simple_extent_dims(sp, dims, nullptr);
    if (dims[1] != expectCols)
      throw std::runtime_error(std::string("matpkg: bad shape for ") + path);
    out.resize(dims[0]);
    if constexpr (std::is_integral<U>::value) {
      std::vector<std::int64_t> buf(dims[0] * N);
      H5Dread(ds, H5T_NATIVE_INT64, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf.data());
      for (hsize_t i = 0; i < dims[0]; i++)
        for (std::size_t j = 0; j < N; j++)
          out[i][j] = static_cast<U>(buf[i * N + j]);
    } else {
      std::vector<double> buf(dims[0] * N);
      H5Dread(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf.data());
      for (hsize_t i = 0; i < dims[0]; i++)
        for (std::size_t j = 0; j < N; j++)
          out[i][j] = static_cast<U>(buf[i * N + j]);
    }
    H5Sclose(sp);
    H5Dclose(ds);
  }

  void readBands(hid_t file) {
    hid_t ds = openDS(file, "/bands/electron/energies");
    hid_t sp = H5Dget_space(ds);
    hsize_t dims[2];
    H5Sget_simple_extent_dims(sp, dims, nullptr);
    std::vector<double> buf(dims[0] * dims[1]);
    H5Dread(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf.data());
    H5Sclose(sp);
    H5Dclose(ds);
    energies.resize(dims[0]);
    bandMin.resize(dims[0]);
    for (hsize_t b = 0; b < dims[0]; b++) {
      energies[b].assign(buf.begin() + b * dims[1],
                         buf.begin() + (b + 1) * dims[1]);
      bandMin[b] = *std::min_element(energies[b].begin(), energies[b].end());
    }
    if (!energies.empty() && energies[0].size() != points.size())
      throw std::runtime_error("matpkg: energies/points size mismatch");
    readVertexVelocities(file);
  }

  /// optional /bands/electron/velocities (Nb, Npt, 3) [m/s]
  void readVertexVelocities(hid_t file) {
    hid_t ds = H5Dopen2(file, "/bands/electron/velocities", H5P_DEFAULT);
    if (ds < 0)
      return;                       // absent: fall back to per-tet gradients
    hid_t sp = H5Dget_space(ds);
    hsize_t d[3];
    if (H5Sget_simple_extent_dims(sp, d, nullptr) != 3 ||
        d[1] != points.size() || d[2] != 3) {
      H5Sclose(sp);
      H5Dclose(ds);
      return;
    }
    std::vector<double> buf(d[0] * d[1] * 3);
    H5Dread(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf.data());
    H5Sclose(sp);
    H5Dclose(ds);
    vertexVel.assign(d[0], {});
    for (hsize_t b = 0; b < d[0]; b++) {
      vertexVel[b].resize(d[1]);
      for (hsize_t i = 0; i < d[1]; i++)
        for (int c = 0; c < 3; c++)
          vertexVel[b][i][c] = static_cast<T>(buf[(b * d[1] + i) * 3 + c]);
    }
  }

  std::vector<std::int64_t> ibzMap;
  std::vector<std::int64_t> ibzOpMap;
  std::vector<std::array<std::array<T, 3>, 3>> symOps;

  void readIbzMap(hid_t file) {
    hid_t ds = H5Dopen2(file, "/kmesh/ibz/map_full_to_ibz", H5P_DEFAULT);
    if (ds < 0)
      return; // optional for the engine (required by the validator)
    hid_t sp = H5Dget_space(ds);
    hsize_t n;
    H5Sget_simple_extent_dims(sp, &n, nullptr);
    ibzMap.resize(n);
    H5Dread(ds, H5T_NATIVE_INT64, H5S_ALL, H5S_ALL, H5P_DEFAULT, ibzMap.data());
    H5Sclose(sp);
    H5Dclose(ds);
    ds = H5Dopen2(file, "/kmesh/ibz/map_sym_op", H5P_DEFAULT);
    if (ds >= 0) {
      sp = H5Dget_space(ds);
      H5Sget_simple_extent_dims(sp, &n, nullptr);
      ibzOpMap.resize(n);
      H5Dread(ds, H5T_NATIVE_INT64, H5S_ALL, H5S_ALL, H5P_DEFAULT,
              ibzOpMap.data());
      H5Sclose(sp);
      H5Dclose(ds);
    }
    ds = H5Dopen2(file, "/lattice/symmetry/rotations", H5P_DEFAULT);
    if (ds >= 0) {
      sp = H5Dget_space(ds);
      hsize_t dims[3];
      H5Sget_simple_extent_dims(sp, dims, nullptr);
      std::vector<std::int64_t> buf(dims[0] * 9);
      H5Dread(ds, H5T_NATIVE_INT64, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf.data());
      symOps.resize(dims[0]);
      for (hsize_t s2 = 0; s2 < dims[0]; s2++)
        for (int r = 0; r < 3; r++)
          for (int c = 0; c < 3; c++)
            symOps[s2][r][c] = static_cast<T>(buf[(s2 * 3 + r) * 3 + c]);
      H5Sclose(sp);
      H5Dclose(ds);
    }
  }

  void readDOS(hid_t file) {
    for (auto [path, vec] :
         {std::pair<const char *, std::vector<T> *>{
              "/bands/electron/dos/energy_grid", &dosGrid},
          std::pair<const char *, std::vector<T> *>{
              "/bands/electron/dos/total", &dosTotal}}) {
      hid_t ds = openDS(file, path);
      hid_t sp = H5Dget_space(ds);
      hsize_t n;
      H5Sget_simple_extent_dims(sp, &n, nullptr);
      std::vector<double> buf(n);
      H5Dread(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf.data());
      vec->assign(buf.begin(), buf.end());
      H5Sclose(sp);
      H5Dclose(ds);
    }
  }
};

#endif // EMC_NUMERIC_BAND_STRUCTURE_HPP
