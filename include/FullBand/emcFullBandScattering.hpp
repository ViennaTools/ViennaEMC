#ifndef EMC_FULL_BAND_SCATTERING_HPP
#define EMC_FULL_BAND_SCATTERING_HPP

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include <hdf5.h>

#include <FullBand/emcNumericBandStructure.hpp>

/*! \brief Table-driven full-band scattering from a .matpkg package
 * (spec v0.1, /scattering/phaseA).
 *
 * Loads energy-resolved rate tables Gamma_nu(T, band, E) and provides:
 *  - total scatter rate Gamma(E) (linear interpolation, clamped - the spec
 *    forbids extrapolation),
 *  - self-scattering constant Gamma0 = safety * max(Gamma),
 *  - mechanism selection proportional to partial rates,
 *  - DOS-weighted final-state selection: k' is sampled on the constant-energy
 *    surface E' of the tetrahedral mesh (per-tet weights ~ tet DOS, uniform
 *    sampling on the iso-surface polygon inside the chosen tetrahedron).
 *
 * Phase-A tables are isotropic in the coupling; anisotropy enters only through
 * the band structure itself (spec 8.1). Phase-B (k-resolved) replaces the
 * final-state sampler, not this interface.
 */
template <class T> class emcFullBandScattering {
public:
  using Vec3 = std::array<T, 3>;

  struct Mechanism {
    std::string name;
    T deltaE;               // energy change per event [eV]
    std::vector<T> grid;    // energy grid [eV] (phaseA)
    std::vector<T> rates;   // [1/s] on grid (phaseA)
    // phaseB anisotropy factor per full-mesh point: a(p) = Gamma_B(p) /
    // Gamma_A(E(p)), dimensionless. The sharp energy dependence stays in the
    // phaseA tables (evaluated at the particle's exact energy); barycentric
    // interpolation of raw point RATES would smear thresholds over the tet's
    // vertex-energy spread and break detailed balance (+10 meV measured).
    std::vector<T> ptRates;
    // Stage-2 |g|^2 final-state rows (spec v0.3): CSR over final tets
    std::vector<std::int64_t> g2Offsets;  // (Nibz+1)
    std::vector<std::int64_t> g2Tets;
    std::vector<T> g2Cum;                 // per-row cumulative weights
  };

  bool hasG2Tables() const { return g2Loaded; }

  /// mode: 0 = auto (best available), 1 = force phaseA energy tables,
  ///       2 = phaseB k-rates without g2 final states
  emcFullBandScattering(const std::string &packageFile,
                        const emcNumericBandStructure<T> &bandStructure,
                        std::size_t band, T temperature, T maxEnergy,
                        T binWidth = static_cast<T>(2e-3), int mode = 0)
      : bs(bandStructure), bandIdx(band) {
    loadMechanisms(packageFile, temperature);
    if (mode != 1)
      kResolved = loadPhaseB(packageFile, temperature);
    buildEnergyBins(maxEnergy, binWidth);
    if (kResolved && mode == 0)
      g2Loaded = loadG2Tables(packageFile);
  }

  /// Stage-2 final-state sampling: draw k' from the |g|^2-weighted row of
  /// the particle's IBZ class for this mechanism; falls back to the
  /// isotropic (phaseA-style) sampler when no row entry straddles E'.
  template <class RNG>
  bool sampleFinalStateG2(std::size_t mech, const Vec3 &kCart, T energyFinal,
                          RNG &rng, Vec3 &kOut, std::int64_t &tetHint) const {
    if (!g2Loaded)
      return sampleFinalState(energyFinal, rng, kOut, tetHint);
    // current IBZ class via dominant vertex of the containing tet
    Vec3 lam;
    std::int64_t th = tetHint;
    const std::int64_t t0 = bs.locateTet(kCart, th, lam);
    const auto &tet = bs.getTetrahedra()[t0];
    const T l0 = 1 - lam[0] - lam[1] - lam[2];
    int v = 0;
    T best = l0;
    for (int i = 0; i < 3; i++)
      if (lam[i] > best) { best = lam[i]; v = i + 1; }
    const std::int64_t ibz = bs.getIbzMap()[tet[v]];
    const std::int64_t iop = bs.getIbzOpMap()[tet[v]];
    const auto &m = mechanisms[mech];
    const std::int64_t r0 = m.g2Offsets[ibz], r1 = m.g2Offsets[ibz + 1];
    if (r1 <= r0)
      return sampleFinalState(energyFinal, rng, kOut, tetHint);
    // draw from cumulative weights, then find a straddling tet
    std::uniform_real_distribution<T> U(0, m.g2Cum[r1 - 1]);
    const T r = U(rng);
    std::int64_t pick = r0;
    while (pick < r1 - 1 && m.g2Cum[pick] < r)
      pick++;
    Vec3 fOut;
    for (std::int64_t off = 0; off < r1 - r0; off++) {
      const std::int64_t cand = r0 + (pick - r0 + off) % (r1 - r0);
      if (samplePointOnIsoFrac(m.g2Tets[cand], energyFinal, rng, fOut)) {
        // tables live in the IBZ representative's frame; rotate the sampled
        // final state into this particle's frame (R : k_ibz -> k_full)
        const auto &R = bs.getSymOps()[iop];
        Vec3 fr{};
        for (int i = 0; i < 3; i++)
          fr[i] = R[i][0] * fOut[0] + R[i][1] * fOut[1] + R[i][2] * fOut[2];
        kOut = bs.fracToCart(fr);
        tetHint = -1; // rotated frame: let the locator re-seed
        g2Hits.fetch_add(1, std::memory_order_relaxed);
        return true;
      }
    }
    g2Fallbacks.fetch_add(1, std::memory_order_relaxed);
    return sampleFinalState(energyFinal, rng, kOut, tetHint);
  }

  std::size_t getG2Hits() const { return g2Hits.load(); }
  std::size_t getG2Fallbacks() const { return g2Fallbacks.load(); }

  bool isKResolved() const { return kResolved; }

  /// total rate at Cartesian k [1/m]: phaseB = exact-energy phaseA lookup
  /// modulated by the barycentric anisotropy factor, else plain phaseA
  T getTotalRate(const Vec3 &kCart, T energy, std::int64_t &tetHint) const {
    if (!kResolved)
      return getTotalRate(energy);
    Vec3 lam;
    const std::int64_t t = bs.locateTet(kCart, tetHint, lam);
    T sum = 0;
    for (const auto &m : mechanisms)
      sum += interp(m, energy) * interpPt(m, t, lam);
    return sum;
  }

  template <class RNG>
  std::size_t selectMechanism(const Vec3 &kCart, T energy, std::int64_t &tetHint,
                              RNG &rng) const {
    if (!kResolved)
      return selectMechanism(energy, rng);
    Vec3 lam;
    const std::int64_t t = bs.locateTet(kCart, tetHint, lam);
    T tot = 0;
    for (const auto &m : mechanisms)
      tot += interp(m, energy) * interpPt(m, t, lam);
    std::uniform_real_distribution<T> U(0, tot);
    T r = U(rng), acc = 0;
    for (std::size_t i = 0; i < mechanisms.size(); i++) {
      acc += interp(mechanisms[i], energy) * interpPt(mechanisms[i], t, lam);
      if (r <= acc)
        return i;
    }
    return mechanisms.size() - 1;
  }

  /// total scatter rate at energy [eV] -> [1/s]; clamped at table ends
  T getTotalRate(T energy) const {
    T sum = 0;
    for (const auto &m : mechanisms)
      sum += interp(m, energy);
    return sum;
  }

  /// self-scattering constant (>= max total rate on the tables)
  T getGamma0(T safety = static_cast<T>(1.2)) const {
    T mx = 0;
    if (kResolved) {
      // at mesh points the product a(p) * Gamma_A(E(p)) equals the raw
      // package rate, so this recovers the true pointwise maximum
      const auto &Ept = bs.getBandEnergies(bandIdx);
      const std::size_t npt = mechanisms[0].ptRates.size();
      for (std::size_t p = 0; p < npt; p++) {
        T s = 0;
        for (const auto &m : mechanisms)
          s += m.ptRates[p] * interp(m, Ept[p]);
        mx = std::max(mx, s);
      }
    } else {
      for (std::size_t i = 0; i < mechanisms[0].grid.size(); i++) {
        T s = 0;
        for (const auto &m : mechanisms)
          s += m.rates[std::min(i, m.rates.size() - 1)];
        mx = std::max(mx, s);
      }
    }
    return safety * mx;
  }

  /// selects mechanism index proportional to partial rates at given energy
  template <class RNG> std::size_t selectMechanism(T energy, RNG &rng) const {
    std::uniform_real_distribution<T> U(0, getTotalRate(energy));
    T r = U(rng), acc = 0;
    for (std::size_t i = 0; i < mechanisms.size(); i++) {
      acc += interp(mechanisms[i], energy);
      if (r <= acc)
        return i;
    }
    return mechanisms.size() - 1;
  }

  T getDeltaE(std::size_t mech) const { return mechanisms[mech].deltaE; }
  std::size_t getNrMechanisms() const { return mechanisms.size(); }
  const std::string &getMechanismName(std::size_t i) const {
    return mechanisms[i].name;
  }

  /// samples k' (Cartesian [1/m]) on the constant-energy surface at energy
  /// [eV]; returns false if the energy is outside the binned band range
  template <class RNG>
  bool sampleFinalState(T energy, RNG &rng, Vec3 &kOut,
                        std::int64_t &tetHint) const {
    const std::int64_t bin = binOf(energy);
    if (bin < 0 || bins[bin].tets.empty())
      return false;
    const auto &B = bins[bin];
    // pick tet ~ weight
    std::uniform_real_distribution<T> U(0, B.totalWeight);
    T r = U(rng), acc = 0;
    std::size_t pick = B.tets.size() - 1;
    for (std::size_t i = 0; i < B.tets.size(); i++) {
      acc += B.weights[i];
      if (r <= acc) {
        pick = i;
        break;
      }
    }
    const std::int64_t t = B.tets[pick];
    if (samplePointOnIso(t, energy, rng, kOut)) {
      tetHint = t;
      return true;
    }
    // degenerate polygon: centroid fallback
    {
      const auto &tetc = bs.getTetrahedra()[t];
      const auto &ptsc = bs.getPoints();
      Vec3 c{};
      for (int i = 0; i < 4; i++)
        for (int j = 0; j < 3; j++)
          c[j] += ptsc[tetc[i]][j] / 4;
      kOut = bs.fracToCart(c);
      tetHint = t;
      return true;
    }
  }

  /// uniform point on the iso-surface polygon of `energy` inside tet t;
  /// false if the tet does not straddle that energy
  template <class RNG>
  bool samplePointOnIso(std::int64_t t, T energy, RNG &rng, Vec3 &kOut) const {
    Vec3 f;
    if (!samplePointOnIsoFrac(t, energy, rng, f))
      return false;
    kOut = bs.fracToCart(f);
    return true;
  }

  /// as samplePointOnIso but returns FRACTIONAL coordinates
  template <class RNG>
  bool samplePointOnIsoFrac(std::int64_t t, T energy, RNG &rng,
                            Vec3 &fOut) const {
    // iso-surface polygon of level `energy` inside tet t (in frac coords)
    const auto &tet = bs.getTetrahedra()[t];
    const auto &pts = bs.getPoints();
    const auto &E = bs.getBandEnergies(bandIdx);
    std::array<Vec3, 4> v;
    std::array<T, 4> e;
    for (int i = 0; i < 4; i++) {
      v[i] = pts[tet[i]];
      e[i] = E[tet[i]];
    }
    Vec3 poly[4];
    int np = 0;
    static const int edges[6][2] = {{0, 1}, {0, 2}, {0, 3},
                                    {1, 2}, {1, 3}, {2, 3}};
    for (const auto &ed : edges) {
      T f;
      // solved consistently with the engine's energy interpolation
      if (bs.edgeCrossing(t, bandIdx, ed[0], ed[1], energy, f)) {
        for (int c = 0; c < 3; c++)
          poly[np][c] = v[ed[0]][c] + f * (v[ed[1]][c] - v[ed[0]][c]);
        if (++np == 4)
          break;
      }
    }
    if (np < 3)
      return false; // tet does not straddle this energy (or grazes a vertex)
    if (np == 4) // order quad vertices to avoid bow-tie (swap if needed)
      orderQuad(poly);
    // triangulate fan (1 or 2 triangles), pick by area, sample uniformly
    const int ntri = np - 2;
    T areas[2] = {0, 0};
    for (int i = 0; i < ntri; i++)
      areas[i] = triArea(poly[0], poly[i + 1], poly[i + 2]);
    std::uniform_real_distribution<T> U01(0, 1);
    const int tri =
        (ntri == 2 && U01(rng) * (areas[0] + areas[1]) > areas[0]) ? 1 : 0;
    T a = U01(rng), b2 = U01(rng);
    if (a + b2 > 1) {
      a = 1 - a;
      b2 = 1 - b2;
    }
    for (int c = 0; c < 3; c++)
      fOut[c] = poly[0][c] + a * (poly[tri + 1][c] - poly[0][c]) +
                b2 * (poly[tri + 2][c] - poly[0][c]);
    return true;
  }

private:
  const emcNumericBandStructure<T> &bs;
  std::size_t bandIdx;
  std::vector<Mechanism> mechanisms;
  bool kResolved = false;
  bool g2Loaded = false;
  // atomic: sampling is const and may run from several ensemble threads
  mutable std::atomic<std::size_t> g2Hits{0}, g2Fallbacks{0};

  /// loads /scattering/g2bins (spec v0.3): CSR rows per mechanism by name
  bool loadG2Tables(const std::string &file) {
    hid_t f = H5Fopen(file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (f < 0)
      return false;
    hid_t grp = H5Gopen2(f, "/scattering/g2bins", H5P_DEFAULT);
    if (grp < 0) {
      H5Fclose(f);
      return false;
    }
    bool any = false;
    for (auto &m : mechanisms) {
      hid_t mg = H5Gopen2(grp, m.name.c_str(), H5P_DEFAULT);
      if (mg < 0)
        continue;
      auto read1 = [&](const char *n, auto &vec, hid_t memType) {
        hid_t ds = H5Dopen2(mg, n, H5P_DEFAULT);
        hid_t sp = H5Dget_space(ds);
        hsize_t cnt;
        H5Sget_simple_extent_dims(sp, &cnt, nullptr);
        vec.resize(cnt);
        H5Dread(ds, memType, H5S_ALL, H5S_ALL, H5P_DEFAULT, vec.data());
        H5Sclose(sp);
        H5Dclose(ds);
      };
      read1("row_offsets", m.g2Offsets, H5T_NATIVE_INT64);
      read1("tet_ids", m.g2Tets, H5T_NATIVE_INT64);
      std::vector<double> w;
      read1("weights", w, H5T_NATIVE_DOUBLE);
      m.g2Cum.resize(w.size());
      // per-row cumulative sums
      for (std::size_t r = 0; r + 1 < m.g2Offsets.size(); r++) {
        T acc = 0;
        for (std::int64_t i = m.g2Offsets[r]; i < m.g2Offsets[r + 1]; i++) {
          acc += static_cast<T>(w[i]);
          m.g2Cum[i] = acc;
        }
      }
      H5Gclose(mg);
      any = true;
    }
    H5Gclose(grp);
    H5Fclose(f);
    return any;
  }

  /// barycentric interpolation of per-point rates within tet t
  T interpPt(const Mechanism &m, std::int64_t t, const Vec3 &lam) const {
    const auto &tet = bs.getTetrahedra()[t];
    const T l0 = 1 - lam[0] - lam[1] - lam[2];
    return l0 * m.ptRates[tet[0]] + lam[0] * m.ptRates[tet[1]] +
           lam[1] * m.ptRates[tet[2]] + lam[2] * m.ptRates[tet[3]];
  }

  /// loads /scattering/phaseB (spec v0.2) as per-point anisotropy factors
  /// a(p) = Gamma_B(p) / Gamma_A(E(p)) attached to the phaseA mechanisms
  /// (matched by group name). Returns false if the group is absent.
  bool loadPhaseB(const std::string &file, T /*temperature*/) {
    const auto &ibz = bs.getIbzMap();
    if (ibz.empty())
      return false;
    hid_t f = H5Fopen(file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (f < 0)
      throw std::runtime_error("matpkg: cannot open " + file);
    hid_t grp = H5Gopen2(f, "/scattering/phaseB", H5P_DEFAULT);
    if (grp < 0) {
      H5Fclose(f);
      return false;
    }
    const auto &Ept = bs.getBandEnergies(bandIdx);
    bool any = false;
    H5G_info_t info;
    H5Gget_info(grp, &info);
    for (hsize_t i = 0; i < info.nlinks; i++) {
      char name[256];
      H5Lget_name_by_idx(grp, ".", H5_INDEX_NAME, H5_ITER_NATIVE, i, name,
                         sizeof(name), H5P_DEFAULT);
      Mechanism *mech = nullptr;
      for (auto &m : mechanisms)
        if (m.name == name) {
          mech = &m;
          break;
        }
      if (!mech)
        continue; // phaseB group with no phaseA counterpart: ignore
      hid_t mg = H5Gopen2(grp, name, H5P_DEFAULT);
      hid_t ds = H5Dopen2(mg, "k_rates", H5P_DEFAULT);
      hid_t sp = H5Dget_space(ds);
      hsize_t dims[3];
      H5Sget_simple_extent_dims(sp, dims, nullptr);
      std::vector<double> buf(dims[0] * dims[1] * dims[2]);
      H5Dread(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf.data());
      H5Sclose(sp);
      H5Dclose(ds);
      H5Gclose(mg);
      // temperature index 0 (single-T v0.2), band slice = bandIdx
      const std::size_t b = std::min<std::size_t>(bandIdx, dims[1] - 1);
      mech->ptRates.resize(ibz.size());
      for (std::size_t p = 0; p < ibz.size(); p++) {
        const T raw = static_cast<T>(buf[b * dims[2] + ibz[p]]);
        const T ref = interp(*mech, Ept[p]);
        // ref == 0 only where the phaseA rate vanishes (e.g. emission below
        // threshold); the runtime product is 0 there whatever a is
        mech->ptRates[p] = ref > 0 ? raw / ref : static_cast<T>(1);
      }
      any = true;
    }
    H5Gclose(grp);
    H5Fclose(f);
    if (any) // mechanisms without k-data act isotropically (a = 1)
      for (auto &m : mechanisms)
        if (m.ptRates.empty())
          m.ptRates.assign(ibz.size(), static_cast<T>(1));
    return any;
  }

  struct Bin {
    std::vector<std::int64_t> tets;
    std::vector<T> weights;
    T totalWeight = 0;
  };
  T binLo = 0, binW = 0;
  std::vector<Bin> bins;

  std::int64_t binOf(T energy) const {
    const std::int64_t b = static_cast<std::int64_t>((energy - binLo) / binW);
    return (b >= 0 && b < (std::int64_t)bins.size()) ? b : -1;
  }

  static T interp(const Mechanism &m, T energy) {
    if (energy <= m.grid.front())
      return m.rates.front(); // clamp (spec: no extrapolation)
    if (energy >= m.grid.back())
      return m.rates.back();
    const auto it = std::upper_bound(m.grid.begin(), m.grid.end(), energy) - 1;
    const std::size_t i = it - m.grid.begin();
    const T f = (energy - m.grid[i]) / (m.grid[i + 1] - m.grid[i]);
    return m.rates[i] * (1 - f) + m.rates[i + 1] * f;
  }

  static T triArea(const Vec3 &a, const Vec3 &b, const Vec3 &c) {
    const T u[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
    const T w[3] = {c[0] - a[0], c[1] - a[1], c[2] - a[2]};
    const T x = u[1] * w[2] - u[2] * w[1], y = u[2] * w[0] - u[0] * w[2],
            z = u[0] * w[1] - u[1] * w[0];
    return std::sqrt(x * x + y * y + z * z) / 2;
  }

  /// order 4 coplanar points into a simple quad (angle sort around centroid)
  static void orderQuad(Vec3 *p) {
    Vec3 c{};
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < 3; j++)
        c[j] += p[i][j] / 4;
    // plane basis from first two edges
    Vec3 u{p[1][0] - p[0][0], p[1][1] - p[0][1], p[1][2] - p[0][2]};
    Vec3 n{}; // normal via cross with another edge
    const Vec3 w{p[2][0] - p[0][0], p[2][1] - p[0][1], p[2][2] - p[0][2]};
    n = {u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2],
         u[0] * w[1] - u[1] * w[0]};
    Vec3 v2{n[1] * u[2] - n[2] * u[1], n[2] * u[0] - n[0] * u[2],
            n[0] * u[1] - n[1] * u[0]};
    std::array<std::pair<T, int>, 4> ang;
    for (int i = 0; i < 4; i++) {
      const Vec3 d{p[i][0] - c[0], p[i][1] - c[1], p[i][2] - c[2]};
      const T x = d[0] * u[0] + d[1] * u[1] + d[2] * u[2];
      const T y = d[0] * v2[0] + d[1] * v2[1] + d[2] * v2[2];
      ang[i] = {std::atan2(y, x), i};
    }
    std::sort(ang.begin(), ang.end());
    const Vec3 tmp[4] = {p[ang[0].second], p[ang[1].second], p[ang[2].second],
                         p[ang[3].second]};
    for (int i = 0; i < 4; i++)
      p[i] = tmp[i];
  }

  void loadMechanisms(const std::string &file, T temperature) {
    hid_t f = H5Fopen(file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (f < 0)
      throw std::runtime_error("matpkg: cannot open " + file);
    hid_t grp = H5Gopen2(f, "/scattering/phaseA", H5P_DEFAULT);
    if (grp < 0)
      throw std::runtime_error("matpkg: no /scattering/phaseA group");
    H5G_info_t info;
    H5Gget_info(grp, &info);
    for (hsize_t i = 0; i < info.nlinks; i++) {
      char name[256];
      H5Lget_name_by_idx(grp, ".", H5_INDEX_NAME, H5_ITER_NATIVE, i, name,
                         sizeof(name), H5P_DEFAULT);
      Mechanism m;
      m.name = name;
      hid_t mg = H5Gopen2(grp, name, H5P_DEFAULT);
      hid_t at = H5Aopen(mg, "delta_E", H5P_DEFAULT);
      double dE = 0;
      H5Aread(at, H5T_NATIVE_DOUBLE, &dE);
      H5Aclose(at);
      m.deltaE = static_cast<T>(dE);
      // temperature index: nearest entry on the temperature grid
      std::vector<double> tg = read1D(mg, "temperature_grid");
      std::size_t ti = 0;
      for (std::size_t j = 1; j < tg.size(); j++)
        if (std::abs(tg[j] - temperature) < std::abs(tg[ti] - temperature))
          ti = j;
      std::vector<double> eg = read1D(mg, "energy_grid");
      m.grid.assign(eg.begin(), eg.end());
      // rates shape (NT, Nb, NE) - v0.1 prototype consumes band slice 0
      hid_t ds = H5Dopen2(mg, "rates", H5P_DEFAULT);
      hid_t sp = H5Dget_space(ds);
      hsize_t dims[3];
      H5Sget_simple_extent_dims(sp, dims, nullptr);
      std::vector<double> buf(dims[0] * dims[1] * dims[2]);
      H5Dread(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
              buf.data());
      m.rates.resize(dims[2]);
      for (hsize_t e = 0; e < dims[2]; e++)
        m.rates[e] = static_cast<T>(buf[(ti * dims[1] + 0) * dims[2] + e]);
      H5Sclose(sp);
      H5Dclose(ds);
      H5Gclose(mg);
      mechanisms.push_back(std::move(m));
    }
    H5Gclose(grp);
    H5Fclose(f);
    if (mechanisms.empty())
      throw std::runtime_error("matpkg: no phaseA mechanisms found");
  }

  static std::vector<double> read1D(hid_t loc, const char *name) {
    hid_t ds = H5Dopen2(loc, name, H5P_DEFAULT);
    hid_t sp = H5Dget_space(ds);
    hsize_t n;
    H5Sget_simple_extent_dims(sp, &n, nullptr);
    std::vector<double> out(n);
    H5Dread(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, out.data());
    H5Sclose(sp);
    H5Dclose(ds);
    return out;
  }

  /// per-tet DOS-proportional weights per energy bin:
  /// weight = vol(tet) * overlap(bin, [emin,emax]) / (emax - emin)
  void buildEnergyBins(T maxEnergy, T width) {
    const auto &tets = bs.getTetrahedra();
    const auto &pts = bs.getPoints();
    const auto &E = bs.getBandEnergies(bandIdx);
    binLo = bs.getBandMinimum(bandIdx);
    binW = width;
    const std::size_t nb =
        static_cast<std::size_t>((maxEnergy - binLo) / width) + 1;
    bins.assign(nb, Bin{});
    // energy range of each tet, CONSISTENT with the engine's interpolation:
    // with the quadratic term the band can dip below the vertex minimum
    // (a band minimum inside a tet), so vertex-only ranges would exclude
    // exactly the tets that matter near the band edge.
    static const T probe[11][3] = {
        {0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1},        // vertices
        {.5, 0, 0}, {0, .5, 0}, {0, 0, .5},                // edges from v0
        {.5, .5, 0}, {.5, 0, .5}, {0, .5, .5},             // remaining edges
        {.25, .25, .25}};                                  // centroid
    for (std::int64_t t = 0; t < (std::int64_t)tets.size(); t++) {
      T lo = E[tets[t][0]], hi = lo;
      for (int i = 1; i < 4; i++) {
        lo = std::min(lo, E[tets[t][i]]);
        hi = std::max(hi, E[tets[t][i]]);
      }
      if (bs.hasQuadraticEnergy()) {
        for (const auto &pr : probe) {
          const Vec3 lam{pr[0], pr[1], pr[2]};
          const T ev = bs.tetEnergy(t, bandIdx, lam);
          lo = std::min(lo, ev);
          hi = std::max(hi, ev);
        }
      }
      if (hi <= binLo || lo >= maxEnergy || hi - lo < 1e-12)
        continue;
      // frac-space volume as DOS proxy (constant Jacobian to cartesian)
      const auto &a = pts[tets[t][0]];
      T e1[3], e2[3], e3[3];
      for (int c = 0; c < 3; c++) {
        e1[c] = pts[tets[t][1]][c] - a[c];
        e2[c] = pts[tets[t][2]][c] - a[c];
        e3[c] = pts[tets[t][3]][c] - a[c];
      }
      const T vol = std::abs(
          e1[0] * (e2[1] * e3[2] - e2[2] * e3[1]) -
          e1[1] * (e2[0] * e3[2] - e2[2] * e3[0]) +
          e1[2] * (e2[0] * e3[1] - e2[1] * e3[0])) / 6;
      const std::int64_t b0 = std::max<std::int64_t>(0, binOf(lo));
      const std::int64_t b1 =
          std::min<std::int64_t>(bins.size() - 1, binOf(std::min(hi, maxEnergy)) < 0
                                                      ? bins.size() - 1
                                                      : binOf(std::min(hi, maxEnergy)));
      for (std::int64_t b = b0; b <= b1; b++) {
        const T bl = binLo + b * binW, bh = bl + binW;
        const T ov = std::min(hi, bh) - std::max(lo, bl);
        if (ov <= 0)
          continue;
        const T w = vol * ov / (hi - lo);
        bins[b].tets.push_back(t);
        bins[b].weights.push_back(w);
        bins[b].totalWeight += w;
      }
    }
  }
};

#endif // EMC_FULL_BAND_SCATTERING_HPP
