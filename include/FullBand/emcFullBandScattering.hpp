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

#include <cstdio>
#include <limits>
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
    // WP3-INTERBAND: destination BAND per row entry, parallel to g2Tets.
    // EMPTY means "every final state is in the source band", which is exactly
    // what every package written before spec v0.4 contains - so an old package
    // behaves bit-identically. The converter fills this once it searches other
    // bands' iso-surfaces for final states.
    std::vector<std::int64_t> g2Band;
    std::vector<T> g2Cum;                 // per-row cumulative weights
    // WP3d: screened-Coulomb angular kernel, beta^2 [1/m^2]. Zero for every
    // mechanism except ionized impurities (read from the phaseA attribute
    // `beta_1_per_m` that add_impurity_mechanism.py writes). When set, the
    // final state is placed WITHIN the drawn tet with density 1/(q^2+beta^2)^2
    // rather than uniformly on the iso-polygon - see sampleFinalStateG2.
    T coulombBeta2 = 0;
    // WP3d: destination band of a Coulomb mechanism (phaseA attr `dest_band`);
    // SIZE_MAX means the source band. With the Bloch overlap in the tables,
    // one mechanism per destination band carries the intra/interband split.
    std::size_t coulombDestBand = SIZE_MAX;
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
    if (kResolved) {
      // rejection bound for the a(k')-weighted final-state sampler
      aMaxMech.assign(mechanisms.size(), T(0));
      for (std::size_t j = 0; j < mechanisms.size(); j++)
        for (T v : mechanisms[j].ptRates)
          aMaxMech[j] = std::max(aMaxMech[j], v);
    }
    loadDos(packageFile, band);
    buildEnergyBins(maxEnergy, binWidth);
    buildValleyFrames();
    if (kResolved && mode == 0)
      g2Loaded = loadG2Tables(packageFile);
  }

  /// Stage-2 final-state sampling: draw k' from the |g|^2-weighted row of
  /// the particle's IBZ class for this mechanism; falls back to the
  /// isotropic (phaseA-style) sampler when no row entry straddles E'.
  /*! \brief Final state for a screened-Coulomb mechanism: KERNEL-CORRECTED
   * rejection sampling around the particle's actual k, landed exactly on the
   * numeric band. Nothing about the band is analytic here - only the PROPOSAL
   * distribution is, and rejection sampling returns the target distribution
   * |V(k,k')|^2 x (surface density of states) for any proposal whose weight is
   * bounded; the proposal only sets the acceptance rate.
   *
   * WHY NOT ROWS. The g2 rows belong to the IBZ representative of the tet's
   * dominant vertex, so for a kernel peaked below the mesh cell the forward
   * core is selected around k_rep, up to one cell from the particle - at
   * N_I = 1e17 that is ~3 screening lengths. Measured on the parabolic
   * package against the exact BH+phonon mobility 595.3: rows with uniform
   * placement 512-530, with the within-tet kernel correction 547-549, and
   * the event-averaged <cos theta> 0.606 against the exact 0.711. No
   * mesh-side fix reaches that; this is the field's full-band Coulomb
   * treatment (Fischetti-Laux lineage).
   *
   * HOW. (1) nearest stamped valley minimum k_min, r = |k - k_min|;
   * (2) theta from the EXACT Brooks-Herring inverse CDF on the sphere of
   * radius r, azimuth uniform; rotate (k - k_min); (3) rescale along that
   * ray by a bracket-and-bisect on the numeric band so E(k') = E exactly;
   * (4) accept with weight w = (|k'-k_min|/r)^3 * f(q_actual)/f(q_proposed)
   * against a cap: the r^3 factor is the angular density of states of a
   * non-spherical surface relative to the sphere, the kernel ratio corrects
   * the proposal's q. For a parabolic isotropic band w == 1, so the
   * parabolic test is exact by construction; for silicon's ellipsoids the
   * weight is correct but the cap costs efficiency (Herring-Vogt proposal is
   * the follow-up). Weight above the cap is COUNTED - a bias indicator.
   */
  template <class RNG>
  bool sampleFinalStateCoulomb(std::size_t mech, const Vec3 &kCart, T energyFinal,
                               RNG &rng, Vec3 &kOut, std::int64_t &tetHint,
                               int *bandOut = nullptr) const {
    const std::size_t nd = mechanisms[mech].coulombDestBand == SIZE_MAX
                               ? bandIdx : mechanisms[mech].coulombDestBand;
    if (nd >= valleysByBand.size() || valleysByBand[nd].empty() || nd > maxDestBand)
      return false;
    const auto &valleys = valleysByBand[nd];
    const bool valleysIsotropic = isoByBand[nd];
    const T b2 = mechanisms[mech].coulombBeta2;
    // Bloch overlap needs the mesh point nearest the SOURCE state
    std::int64_t hSrc = tetHint;
    const std::int64_t vSrc = bs.hasGauge ? bs.nearestVertex(kCart, hSrc) : -1;
    // nearest valley image (frames hold cartesian minima; minimum image in frac)
    const Vec3 kf = bs.cartToFrac(kCart);
    std::size_t vb = 0; Vec3 dbest{}; T r2 = std::numeric_limits<T>::max();
    for (std::size_t v = 0; v < valleys.size(); v++) {
      const Vec3 mf = bs.cartToFrac(valleys[v].kmin);
      Vec3 d{};
      for (int i = 0; i < 3; i++) { d[i] = kf[i] - mf[i]; d[i] -= std::round(d[i]); }
      const Vec3 dc = bs.fracToCart(d);
      const T q = dc[0]*dc[0] + dc[1]*dc[1] + dc[2]*dc[2];
      if (q < r2) { r2 = q; dbest = dc; vb = v; }
    }
    if (!(r2 > 0) || r2 == std::numeric_limits<T>::max()) return false;
    const auto &V = valleys[vb];
    const T rEff = std::sqrt(r2);                    // real radius, for the BH draw
    Vec3 kmin{};
    for (int i = 0; i < 3; i++) kmin[i] = kCart[i] - dbest[i];
    // HERRING-VOGT: in k* = M^{-1/2}(k - k_min) the stamped ellipsoid is a
    // sphere and the angular density of states is UNIFORM, so a direction
    // proposed there needs no r^3 Jacobian; only the kernel mismatch between
    // the proposal's q and the real q is corrected. Reduces to the plain
    // sphere for an isotropic band.
    Vec3 ks{};
    for (int r = 0; r < 3; r++) ks[r] = V.mih[r][0]*dbest[0] + V.mih[r][1]*dbest[1] + V.mih[r][2]*dbest[2];
    const T rho = std::sqrt(ks[0]*ks[0] + ks[1]*ks[1] + ks[2]*ks[2]);
    if (!(rho > 0)) return false;
    Vec3 u0{}; for (int i = 0; i < 3; i++) u0[i] = ks[i] / rho;
    Vec3 a{}; a[std::abs(u0[0]) < 0.6 ? 0 : (std::abs(u0[1]) < 0.6 ? 1 : 2)] = 1;
    Vec3 e1{u0[1]*a[2]-u0[2]*a[1], u0[2]*a[0]-u0[0]*a[2], u0[0]*a[1]-u0[1]*a[0]};
    { const T n = std::sqrt(e1[0]*e1[0]+e1[1]*e1[1]+e1[2]*e1[2]); for (int i=0;i<3;i++) e1[i]/=n; }
    const Vec3 e2{u0[1]*e1[2]-u0[2]*e1[1], u0[2]*e1[0]-u0[0]*e1[2], u0[0]*e1[1]-u0[1]*e1[0]};
    // COULOMB_CAP overrides the importance cap (A/B for the truncation bias)
    static const T capEnv = std::getenv("COULOMB_CAP") ? static_cast<T>(std::atof(std::getenv("COULOMB_CAP"))) : T(0);
    // Default cap 8 for anisotropic valleys. Measured on Si 1e17 with a
    // 512-proposal budget (no fallbacks): cap 4 -> 881 +/- 28, cap 8 -> 861
    // +/- 13 (0.38% of proposals capped), cap 16 -> 836 +/- 29 (0.00%). The
    // truncation at cap 4 biases mu high by ~3-4%; 8 is within the SEs of 16
    // at half the cost. The tail comes from the (sm/s0)^3 Jacobian where the
    // real surface bulges beyond the stamped ellipsoid (non-parabolicity); a
    // Kane-corrected proposal radius would shrink it - not done.
    const T wcap = capEnv > 0 ? capEnv : (valleysIsotropic ? T(2) : T(8));
    std::uniform_real_distribution<T> U01(0, 1);
    // proposal budget: exhausting it falls back to the ROWS, which are the
    // smeared, over-relaxing path - measured on Si 1e17: cap 16 with a budget
    // of 64 sent 4% of events to the rows and pulled mu from 840 to 793.
    // Spending more proposals is always cheaper than that.
    static const int maxProp = std::getenv("COULOMB_MAXPROP") ? std::atoi(std::getenv("COULOMB_MAXPROP")) : 512;
    for (int prop = 0; prop < maxProp; prop++) {
      coulombProposals.fetch_add(1, std::memory_order_relaxed);
      // Brooks-Herring theta on a sphere of radius rEff: x = 1 - cos(theta)
      const T u = U01(rng);
      const T A = 1 / b2 - u * (1 / b2 - 1 / (4 * r2 + b2));
      const T x = std::min<T>(2, std::max<T>(0, (1 / A - b2) / (2 * r2)));
      const T ct = 1 - x, st = std::sqrt(std::max<T>(0, 1 - ct * ct));
      const T phi = 2 * static_cast<T>(M_PI) * U01(rng);
      Vec3 os{};                                     // direction in k*
      for (int i = 0; i < 3; i++)
        os[i] = ct * u0[i] + st * (std::cos(phi) * e1[i] + std::sin(phi) * e2[i]);
      Vec3 d0{};                                     // k'_0 - k_min on the ellipsoid
      for (int r = 0; r < 3; r++)
        d0[r] = rho * (V.mh[r][0]*os[0] + V.mh[r][1]*os[1] + V.mh[r][2]*os[2]);
      const T s0 = std::sqrt(d0[0]*d0[0] + d0[1]*d0[1] + d0[2]*d0[2]);
      if (!(s0 > 0)) continue;
      Vec3 ud{}; for (int i = 0; i < 3; i++) ud[i] = d0[i] / s0;
      // land on the numeric band along the ray from k_min
      std::int64_t h = tetHint;
      auto Eat = [&](T sc) {
        Vec3 k{}; for (int i = 0; i < 3; i++) k[i] = kmin[i] + sc * ud[i];
        return bs.getEnergy(k, nd, h);
      };
      T sLo = s0, sHi = s0, eLo = Eat(s0), eHi = eLo;
      bool ok = true;
      if (eLo > energyFinal) {
        int it = 0;
        while (eLo > energyFinal && it++ < 14) { sLo *= T(0.7); eLo = Eat(sLo); }
        ok = eLo <= energyFinal;
      } else {
        int it = 0;
        while (eHi < energyFinal && it++ < 14) { sHi *= T(1.4); eHi = Eat(sHi); }
        ok = eHi >= energyFinal;
      }
      if (!ok) { coulombBracketFail.fetch_add(1, std::memory_order_relaxed); continue; }
      T sm = s0;
      for (int it = 0; it < 40; it++) {
        sm = T(0.5) * (sLo + sHi);
        const T em = Eat(sm);
        if (std::abs(em - energyFinal) < T(1e-8)) break;
        if (em < energyFinal) sLo = sm; else sHi = sm;
      }
      Vec3 kp{};
      for (int i = 0; i < 3; i++) kp[i] = kmin[i] + sm * ud[i];
      // importance weight = target / proposal:
      //   kernel ratio   - proposal q on the rEff sphere vs the real q;
      //   (sm/s0)^3      - the angular density of states of the REAL surface
      //                    relative to the stamped ellipsoid, i.e. the exact
      //                    Jacobian of projecting the ellipsoid onto the numeric
      //                    surface along the ray from k_min. With it the fast
      //                    path is exact for ANY surface star-shaped about the
      //                    valley minimum, not only for a parabolic ellipsoid;
      //                    on one it is 1 up to mesh noise.
      const T qp2 = 2 * r2 * x;
      T qa2 = 0;
      for (int i = 0; i < 3; i++) qa2 += (kp[i] - kCart[i]) * (kp[i] - kCart[i]);
      const T ratio = (qp2 + b2) / (qa2 + b2);
      const T jac = (sm / s0) * (sm / s0) * (sm / s0);
      // Bloch overlap |<u_nd,k'|u_n,k>|^2 (<= 1): the rate tables already
      // carry it, so it enters the acceptance as a plain factor
      T I = 1;
      if (vSrc >= 0) {
        std::int64_t hd = h;
        const std::int64_t vDst = bs.nearestVertex(kp, hd);
        if (vDst >= 0) I = bs.overlapGroup(vSrc, bandIdx, vDst, nd);
      }
      const T w = ratio * ratio * jac * I;
      if (w > wcap) coulombCapped.fetch_add(1, std::memory_order_relaxed);
      if (U01(rng) * wcap <= w) {
        kOut = kp;
        tetHint = h;
        if (bandOut) *bandOut = static_cast<int>(nd);
        coulombAnalytic.fetch_add(1, std::memory_order_relaxed);
        return true;
      }
    }
    return false;
  }

  template <class RNG>
  bool sampleFinalStateG2(std::size_t mech, const Vec3 &kCart, T energyFinal,
                          RNG &rng, Vec3 &kOut, std::int64_t &tetHint,
                          int *bandOut = nullptr) const {
    // bandOut, when given, receives the destination BAND of the sampled final
    // state. It stays at the source band unless the package carries dest_band
    // rows, so this is inert on pre-v0.4 packages.
    if (bandOut) *bandOut = static_cast<int>(bandIdx);
    // WP3d: closed-form kernels are placed by kernel-corrected rejection
    // sampling from the particle's actual k - see sampleFinalStateCoulomb. Rows for such a mechanism are
    // only the fallback.
    if (mech < mechanisms.size() && mechanisms[mech].coulombBeta2 > 0 &&
        bs.edgePar.hasStar) {
      if (sampleFinalStateCoulomb(mech, kCart, energyFinal, rng, kOut, tetHint, bandOut))
        return true;
      coulombFallbacks.fetch_add(1, std::memory_order_relaxed);
    }
    if (!g2Loaded)
      return sampleFinalState(energyFinal, rng, kOut, tetHint, mech);
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
    // A mechanism WITHOUT g2 rows in a package where others have them (an
    // ionized-impurity mechanism added to a phaseA-only package, say) leaves
    // g2Offsets EMPTY - loadG2Tables skips absent groups. Indexing it was
    // undefined behaviour; every v0.4 package so far carried rows for every
    // mechanism, so it never fired. Such a mechanism samples DOS-weighted.
    if (m.g2Offsets.size() < static_cast<std::size_t>(ibz) + 2)
      return sampleFinalState(energyFinal, rng, kOut, tetHint, mech);
    const std::int64_t r0 = m.g2Offsets[ibz], r1 = m.g2Offsets[ibz + 1];
    if (r1 <= r0)
      return sampleFinalState(energyFinal, rng, kOut, tetHint, mech);
    // REJECTION SAMPLING over the |g|^2 row.
    //
    // The previous version drew ONE |g|^2-weighted tet and, if it did not
    // straddle E', walked FORWARD through the row taking whatever came next in
    // storage order. That is a deterministic substitution: tets that happen to
    // follow a non-straddling one get systematically over-selected, so the
    // final states are not |g|^2-distributed. Measured on the production
    // package it fired for 17.6% of accepted draws with a mean walk of 3.8
    // tets - i.e. nearly one final state in five came from the wrong measure,
    // while the reported "wrong measure" figure counted only the 0.75% that
    // fell through to the isotropic sampler.
    //
    // Redrawing instead yields exactly the correct conditional distribution:
    // |g|^2 restricted to the tets that straddle E'. Rejection sampling is
    // unbiased by construction - each accepted sample is distributed as the
    // proposal conditioned on acceptance.
    std::uniform_real_distribution<T> U(0, m.g2Cum[r1 - 1]);
    Vec3 fOut;
    const int nsub = bs.hasSubdivision() ? 8 : 1;
    std::uniform_int_distribution<int> Usub(0, nsub - 1);
    constexpr int MAXTRY = 64;
    for (int trial = 0; trial < MAXTRY; trial++) {
      const T r = U(rng);
      // binary search the cumulative weights (rows reach ROWCAP=640)
      std::int64_t pick =
          std::lower_bound(m.g2Cum.begin() + r0, m.g2Cum.begin() + r1, r) -
          m.g2Cum.begin();
      if (pick >= r1) pick = r1 - 1;
      const std::int64_t id = bs.hasSubdivision()
                                  ? ((m.g2Tets[pick] << 3) | Usub(rng))
                                  : m.g2Tets[pick];
      // v0.4: the row entry carries its own destination band, and the
      // iso-surface must be resolved against THAT band - see
      // samplePointOnIsoFrac. Without dest_band (pre-v0.4) g2Band is empty
      // and this stays the source band, so old packages are unaffected.
      const std::size_t db =
          pick < static_cast<std::int64_t>(m.g2Band.size())
              ? static_cast<std::size_t>(m.g2Band[pick])
              : bandIdx;
      // The engine may not have a scattering instance for this destination
      // (ENGINE_NBANDS=1 on a v0.4 package). REJECT and redraw - do NOT accept
      // the draw and then move the carrier back to the source band: k' was
      // placed on the DESTINATION band's iso-surface at E', so re-labelling
      // its band leaves it off the energy shell by the inter-band separation.
      // Rejecting instead samples |g|^2 restricted to reachable destinations,
      // renormalised - which is exactly the pre-v0.4 row content, and so makes
      // ENGINE_NBANDS=1 a faithful no-interband CONTROL.
      if (db > maxDestBand) {
        g2Rejected.fetch_add(1, std::memory_order_relaxed);
        continue;
      }
      if (samplePointOnIsoFrac(id, energyFinal, rng, fOut, db)) {
        const auto &R = bs.getSymOps()[iop];
        Vec3 fr{};
        for (int i = 0; i < 3; i++)
          fr[i] = R[i][0] * fOut[0] + R[i][1] * fOut[1] + R[i][2] * fOut[2];
        // KERNEL-AWARE PLACEMENT WITHIN THE TET (ionized impurities). The row
        // weight is the kernel INTEGRATED over the tet's iso-polygon, so the
        // tet marginal is right; but placing k' uniformly on that polygon is
        // not, because 1/(q^2+beta^2)^2 is peaked BELOW the tet scale whenever
        // the screening length is comparable to a mesh cell (N_I = 1e17:
        // 12.9 nm vs 7.7 nm on N48). Uniform placement in the tets around k
        // inflates <q^2> per event ~2.3x and over-relaxes momentum; measured
        // on the parabolic package against the exact BH+phonon mobility
        // (595.3): Einstein 512-530, i.e. -12%. Rejection against the kernel
        // with the TET-LOCAL bound f_max = 1/(q_min^2+beta^2)^2, re-proposing
        // in the SAME tet on rejection, makes the within-tet density
        // proportional to the kernel and leaves the tet marginal untouched.
        // q is the minimum-image distance k' - k in the particle's frame.
        if (m.coulombBeta2 > 0) {
          const T b2 = m.coulombBeta2;
          const Vec3 kf = bs.cartToFrac(kCart);
          const std::int64_t par = bs.hasSubdivision() ? (id >> 3) : id;
          const auto &tv = bs.getTetrahedra()[par];
          const auto &P = bs.getPoints();
          T qmin2 = std::numeric_limits<T>::max();
          for (int v = 0; v < 4; v++) {
            Vec3 pv{}, dv{};
            for (int i = 0; i < 3; i++)
              pv[i] = R[i][0] * P[tv[v]][0] + R[i][1] * P[tv[v]][1] + R[i][2] * P[tv[v]][2];
            for (int i = 0; i < 3; i++) {
              dv[i] = pv[i] - kf[i];
              dv[i] -= std::round(dv[i]);
            }
            const Vec3 dc = bs.fracToCart(dv);
            qmin2 = std::min(qmin2, dc[0]*dc[0] + dc[1]*dc[1] + dc[2]*dc[2]);
          }
          const T fmax = T(1) / ((qmin2 + b2) * (qmin2 + b2));
          std::uniform_real_distribution<T> U01c(0, 1);
          bool ok = false;
          for (int rep = 0; rep < 256; rep++) {
            Vec3 dq{};
            for (int i = 0; i < 3; i++) {
              dq[i] = fr[i] - kf[i];
              dq[i] -= std::round(dq[i]);
            }
            const Vec3 dqc = bs.fracToCart(dq);
            const T q2 = dqc[0]*dqc[0] + dqc[1]*dqc[1] + dqc[2]*dqc[2];
            const T f = T(1) / ((q2 + b2) * (q2 + b2));
            coulombTries.fetch_add(1, std::memory_order_relaxed);
            if (U01c(rng) * fmax <= f) { ok = true; break; }
            // re-propose in the SAME tet (keeps the row-weighted marginal)
            if (!samplePointOnIsoFrac(id, energyFinal, rng, fOut, db)) break;
            for (int i = 0; i < 3; i++)
              fr[i] = R[i][0] * fOut[0] + R[i][1] * fOut[1] + R[i][2] * fOut[2];
          }
          if (!ok) coulombExhausted.fetch_add(1, std::memory_order_relaxed);
        }
        if (bandOut) *bandOut = static_cast<int>(db);
        g2Walk.fetch_add(trial > 0 ? 1 : 0, std::memory_order_relaxed);
        g2WalkLen.fetch_add(static_cast<std::size_t>(trial),
                            std::memory_order_relaxed);
        kOut = bs.fracToCart(fr);
        tetHint = -1; // rotated frame: let the locator re-seed
        g2Hits.fetch_add(1, std::memory_order_relaxed);
        return true;
      }
    }
    g2Fallbacks.fetch_add(1, std::memory_order_relaxed);
    return sampleFinalState(energyFinal, rng, kOut, tetHint, mech);
  }

  /// Draws an energy from the package's OWN equilibrium distribution
  /// D(E) exp(-(E-Emin)/kT). The energy bins already carry the tetrahedral
  /// DOS weight, so this is material-agnostic - no parabolic assumption.
  ///
  /// The driver used to draw E ~ exp(-E/kT) with no density-of-states
  /// factor, which gives <E> = kT instead of 1.5kT for a parabolic band.
  /// With a constant lifetime mu = e*tau/m* is distribution-independent, so
  /// the error stayed invisible on the constant-rate regression package;
  /// with an energy-dependent lifetime - or any elastic package, which can
  /// never relax the energy distribution it was handed - it biases every
  /// transport number the engine reports.
  /// Build the thermal CDF up front, from ONE thread, before any worker runs.
  /// MUST be called before the replica threads start - see the race note on
  /// sampleThermalEnergy below.
  void prepareThermal(T kT) const {
    if (thermKT != kT)
      buildThermalCdf(kT);
  }

  template <class RNG> T sampleThermalEnergy(T kT, RNG &rng) const {
    // DATA RACE if this is the first call and several replica threads reach it
    // together: `thermCum` is a shared mutable vector, so two threads can both
    // see the cache stale and assign() it concurrently, or one can read
    // thermCum.back() while another is mid-assign. That was undefined
    // behaviour in every threaded run, and it surfaced as RUN-TO-RUN VARIATION
    // in mu with identical inputs (1442.2 vs 1408.8) once a second scattering
    // instance shifted thread timing enough to expose it. The driver now calls
    // prepareThermal() before spawning workers; this lazy path remains only as
    // a fallback for single-threaded callers.
    if (thermKT != kT)
      buildThermalCdf(kT);
    std::uniform_real_distribution<T> U(0, thermCum.back());
    const T r = U(rng);
    std::size_t b = static_cast<std::size_t>(
        std::lower_bound(thermCum.begin(), thermCum.end(), r) -
        thermCum.begin());
    if (b >= bins.size())
      b = bins.size() - 1;
    std::uniform_real_distribution<T> Ub(0, binW);
    return binLo + static_cast<T>(b) * binW + Ub(rng);
  }

  /// Enforce detailed balance against the bin DOS the sampler actually
  /// draws from. For a phonon pair (absorption +hw, emission -hw) the
  /// stationary distribution is Boltzmann on the sampler's own measure only if
  ///     Gamma_abs(E) D(E) e^-E/kT == Gamma_em(E+hw) D(E+hw) e^-(E+hw)/kT.
  /// The tables are built against a DIFFERENT density - analytic sqrt(E) for
  /// the test packages, the true ab-initio DOS for silicon - so the condition
  /// is violated and the ensemble relaxes to the wrong temperature (measured:
  /// -2.8% on the parabolic package, +7.0% on Si 64^3, both against the
  /// equilibrium the bins themselves imply).
  ///
  /// The repair is symmetric: both partners are moved to the geometric mean of
  /// the two sides, which fixes the ratio while preserving the pair's overall
  /// scattering strength. It does NOT fix the separate error that the binned
  /// DOS differs from the true band DOS (~+10% for Si at 64^3); that is
  /// discretisation and needs a finer mesh.
  std::size_t enforceDetailedBalance(T kT) {
    const T dosFloor = 1e-30;
    // Density for the balance condition. Rate grids and the package DOS grid
    // both carry ABSOLUTE energies (Si CBM = 6.67 eV); binOf() subtracts
    // binLo itself, so adding binLo here once double-offset every lookup out
    // of range and made the whole repair a silent no-op on any package with a
    // non-zero CBM.
    //
    // Prefer the package's linear-tetrahedron DOS: that IS the density the
    // DOS-weighted iso-surface sampler realises at an exact energy. The bin
    // average is only a 2 meV-smeared stand-in for it, and smearing is worst
    // exactly where the DOS varies fastest - the band edge, which dominates
    // transport. (Measured: bin-averaged repair left +4.8% at 32^3.)
    // MEASURED: the bin average is the RIGHT density here, and the package's
    // exact linear-tetrahedron DOS is worse (Si <E> 0.0427 -> 0.0448 at 64^3,
    // 0.0462 -> 0.0560 at 32^3). The sampler draws a tetrahedron FROM A BIN,
    // so the density it realises is the bin-averaged one; the tet DOS is the
    // better description of the band but not of the sampler. Balance must be
    // enforced against the measure actually sampled. DOSSRC=pkg to A/B it.
    // Use the engine's own BIN DOS. The package's stored
    // /bands/electron/dos/total is defective: on the analytic parabolic
    // packages, where D(E) ~ sqrt(E) exactly and <E> = 1.5kT = 0.0388 eV, it
    // implies 0.0419-0.0427 at EVERY mesh from 32^3 to 80^3 - a flat +8-10%
    // bias that does not converge. The bin DOS implies 0.0392 there (+1%).
    //
    // Do not be fooled by self-consistency: enforcing balance against ANY
    // density D makes the stationary distribution D(E)exp(-E/kT), so the MC
    // will always reproduce the equilibrium of whichever density was used.
    // Each candidate must be judged against an EXTERNAL exact answer, never
    // against its own implied equilibrium.
    // DOSSRC=pkg to A/B (expect it to be ~8% hot until matforge is fixed).
    const char *dsrc = std::getenv("DOSSRC");
    const bool useDos = !dosGrid.empty() && dsrc && std::string(dsrc) == "pkg";
    auto dosOf = [&](T eAbs) -> T {
      if (useDos) {
        if (eAbs <= dosGrid.front() || eAbs >= dosGrid.back())
          return T(0);
        const std::size_t k =
            std::lower_bound(dosGrid.begin(), dosGrid.end(), eAbs) -
            dosGrid.begin();
        const T f = (eAbs - dosGrid[k - 1]) / (dosGrid[k] - dosGrid[k - 1]);
        return dosVals[k - 1] + f * (dosVals[k] - dosVals[k - 1]);
      }
      const std::int64_t b = binOf(eAbs);
      return b < 0 ? T(0) : bins[b].dos / binW;
    };
    auto binDos = dosOf;

    // A-WEIGHTED MEASURE FOR THE K-RESOLVED PATH.
    // The balance condition on the ENERGY distribution is
    //     f(E) Gamma_A,abs(E) Z(E) = f(E') Gamma_A,em(E') Z(E'),
    // where Z(E) = sum_{k at E} a(k) w(k) is the measure the flux actually
    // carries. With a k-resolved out-rate a(k)*Gamma_A(E), Z = <a>*D, NOT D.
    // Repairing against D (as this function used to) is therefore correct only
    // for mode 1, where a == 1. In modes 2/0 it repairs against the wrong
    // measure - which is exactly why BALANCE=1 improved mode 1 slightly
    // (1.067 -> 1.041) while DEGRADING mode 2 (1.005 -> 0.911), and why the
    // modes relaxed to different zero-field equilibria with it off
    // (<E> 0.0433 vs 0.0446 - an equilibrium may not depend on the rate model).
    // a_abs is evaluated at the initial energy E, a_em at its initial E'.
    std::vector<std::vector<T>> zbin;
    if (kResolved) {
      zbin.assign(mechanisms.size(), std::vector<T>(bins.size(), T(0)));
      for (std::size_t m = 0; m < mechanisms.size(); m++) {
        if (mechanisms[m].ptRates.empty()) continue;
        const auto &pr = mechanisms[m].ptRates;
        for (std::size_t b = 0; b < bins.size(); b++) {
          T z = 0;
          const auto &B = bins[b];
          for (std::size_t i = 0; i < B.tets.size(); i++) {
            const std::int64_t par =
                bs.hasSubdivision() ? (B.tets[i] >> 3) : B.tets[i];
            const auto &tv = bs.getTetrahedra()[par];
            z += static_cast<T>(0.25) *
                 (pr[tv[0]] + pr[tv[1]] + pr[tv[2]] + pr[tv[3]]) * B.weights[i];
          }
          zbin[m][b] = z / binW;
        }
      }
    }
    auto zOf = [&](const Mechanism &m, T eAbs) -> T {
      if (!kResolved || m.ptRates.empty() || useDos) return binDos(eAbs);
      const std::int64_t b = binOf(eAbs);
      if (b < 0 || b >= static_cast<std::int64_t>(bins.size())) return T(0);
      const std::size_t mi = static_cast<std::size_t>(&m - &mechanisms[0]);
      return mi < zbin.size() ? zbin[mi][b] : binDos(eAbs);
    };

    std::size_t npair = 0;
    for (auto &ma : mechanisms) {
      if (ma.deltaE <= 0)
        continue;                       // absorption drives the pairing
      const T hw = ma.deltaE;
      Mechanism *me = nullptr;
      for (auto &m : mechanisms)
        if (m.deltaE < 0 && std::fabs(m.deltaE + hw) < 1e-6 * hw &&
            partnerName(m.name) == partnerName(ma.name))
          me = &m;
      if (!me)
        continue;
      npair++;
      nrewrite_ = 0;
      // PASS 1: g(E) = sqrt(lhs*rhs) from the ORIGINAL tables only. Repairing
      // in place instead lets later grid points balance against already-
      // repaired values and lets several source points collapse onto one
      // nearest destination point - which drove <E> further off, not closer.
      std::vector<T> g(ma.grid.size(), T(-1));
      for (std::size_t i = 0; i < ma.grid.size(); i++) {
        const T E = ma.grid[i];
        const T Da = zOf(ma, E), Db = zOf(*me, E + hw);
        if (Da <= dosFloor || Db <= dosFloor)
          continue;
        const T lhs = ma.rates[i] * Da;
        const T rhs = interp(*me, E + hw) * Db * std::exp(-hw / kT);
        if (lhs > 0 && rhs > 0) {
          g[i] = std::sqrt(lhs * rhs);
          nrewrite_++;
        }
      }
      // PASS 2: write both partners from g. Emission is written on its OWN
      // grid, evaluating g at E'-hw, so every destination point is set once.
      std::vector<T> aNew = ma.rates, eNew = me->rates;
      for (std::size_t i = 0; i < ma.grid.size(); i++) {
        const T Da = zOf(ma, ma.grid[i]);
        if (g[i] > 0 && Da > dosFloor)
          aNew[i] = g[i] / Da;
      }
      for (std::size_t j = 0; j < me->grid.size(); j++) {
        const T Ep = me->grid[j];
        const T Db = zOf(*me, Ep);
        if (Db <= dosFloor || Ep - hw < ma.grid.front() ||
            Ep - hw > ma.grid.back())
          continue;
        const T gi = interpArr(ma.grid, g, Ep - hw);
        if (gi > 0)
          eNew[j] = gi * std::exp(hw / kT) / Db;
      }
      // How far the repair moves the AB-INITIO numbers. This is the honest
      // measure of how much "fitting" it does: a legitimate discretisation
      // consistency correction is O(mesh error) and shrinks on refinement.
      // A large value means the rates are being bent to compensate for a bad
      // mesh DOS, which has no ab-initio meaning.
      for (std::size_t i = 0; i < ma.rates.size(); i++)
        if (ma.rates[i] > 0) {
          const T rel = std::fabs(aNew[i] - ma.rates[i]) / ma.rates[i];
          relSum_ += rel; relN_++;
          if (rel > relMax_) relMax_ = rel;
        }
      for (std::size_t i = 0; i < me->rates.size(); i++)
        if (me->rates[i] > 0) {
          const T rel = std::fabs(eNew[i] - me->rates[i]) / me->rates[i];
          relSum_ += rel; relN_++;
          if (rel > relMax_) relMax_ = rel;
        }
      ma.rates = aNew;
      me->rates = eNew;
    }
    return npair;
  }

  std::size_t getRewrittenPoints() const { return nrewrite_; }
  T getRateChangeMean() const { return relN_ ? relSum_ / relN_ : T(0); }
  T getRateChangeMax() const { return relMax_; }

  /// Equilibrium <E>-Emin implied by the engine's OWN energy bins:
  /// sum_b w_b E_b e^-E_b/kT / sum_b w_b e^-E_b/kT, with w_b the same bin
  /// weights the final-state sampler draws from. This is the temperature
  /// the ensemble MUST relax to if detailed balance holds - the only
  /// self-consistent reference for the measured <E>. Comparing against a
  /// DOS quadrature computed some other way (tet means, vertex integrals,
  /// point sums) compares two different discretisations, not the engine
  /// against the truth.
  /// Boltzmann weight of this band's occupiable states, sum_bins dos x exp(-E/kT):
  /// the relative population a THERMAL ensemble puts on this band. Used by the
  /// driver to draw each particle's initial band - every particle used to start
  /// in band 0, which for a degenerate valence manifold is the heavy branch alone.
  T thermalWeight(T kT) const {
    T z = 0;
    for (std::size_t i = 0; i < bins.size(); i++) {
      const T Ec = binLo + (static_cast<T>(i) + T(0.5)) * binW;   // absolute
      z += bins[i].dos * std::exp(-(Ec - bs.getBandMinimum(0)) / kT);
    }
    return z;
  }

  T equilibriumEnergy(T kT) const {
    T num = 0, den = 0;
    for (std::size_t i = 0; i < bins.size(); i++) {
      const T Ec = (static_cast<T>(i) + T(0.5)) * binW;
      const T w = bins[i].dos * std::exp(-Ec / kT);
      num += w * Ec;
      den += w;
    }
    return den > 0 ? num / den : T(0);
  }

  std::size_t getCentroidFallbacks() const {
    return centroidFallbacks.load();
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
  /// Piecewise-constant self-scattering bound: Gamma0 for the slab holding
  /// `energy`, instead of one global maximum over the whole band.
  ///
  /// A single global Gamma0 must bound the rate EVERYWHERE, so on silicon it
  /// sits ~45x above the thermal rate and 97.8% of flights are self-
  /// scattering - 45 wasted flights per real event. Slabbing the energy axis
  /// bounds each region separately and brings that to tens of percent.
  ///
  /// Self-scattering is a variance-reduction device, not physics: any bound
  /// >= the true rate gives identical results. The ONLY correctness
  /// requirement is that the bound is never exceeded, so each slab takes the
  /// max over itself PLUS its neighbours - a flight moves k by ~2% of a cell,
  /// far less than a slab - and violations are counted so a breach cannot
  /// pass silently.
  T gamma0At(T energy) const {
    if (slabG0.empty())
      return gamma0Global;
    const std::int64_t i =
        static_cast<std::int64_t>((energy - slabLo) / slabW);
    if (i < 0 || i >= static_cast<std::int64_t>(slabG0.size()))
      return gamma0Global;
    return slabG0[i];
  }

  /// slab width [eV]. The driver needs it to cap a free flight so the
  /// particle cannot leave the slab whose Gamma0 it drew its flight time
  /// from - without that cap the bound is field-blind and breaks at high F.
  T getSlabWidth() const { return slabG0.empty() ? T(0) : slabW; }

  std::size_t getAWeightTries() const { return aWeightTries.load(); }
  std::size_t getAWeightRej() const { return aWeightRej.load(); }

  std::size_t getG2Walks() const { return g2Walk.load(); }
  /// v0.4 packages carry per-source-band g2 rows; v0.3 do not, and then every
  /// band reads band 0's. Only meaningful once loadG2Tables has run.
  bool hasPerSourceBandG2() const { return g2PerSourceBand; }
  /// restrict final states to destination bands <= b (the bands the driver
  /// actually built instances for). Unset = unrestricted.
  void setMaxDestBand(std::size_t b) { maxDestBand = b; }
  std::size_t getG2Rejected() const { return g2Rejected.load(); }
  std::size_t getG2WalkLen() const { return g2WalkLen.load(); }

  std::size_t getG0Violations() const { return g0Viol.load(); }
  void noteG0Violation() const { g0Viol.fetch_add(1, std::memory_order_relaxed); }

  /// builds the slab table; call once after the rate tables are loaded
  /// 10 meV default: measured on Si, self-scattering 97.8% -> 40.7% with
  /// ZERO bound violations. Narrower slabs cut self-scattering further
  /// (28.5% at 4 meV, 24.0% at 2 meV) but START BREACHING THE BOUND (81 and
  /// 621 violations) because the one-slab neighbour margin no longer covers
  /// how far a flight can move in energy - and a breached bound is wrong
  /// physics, not merely slow. Wall time plateaus at ~6 s either way, since
  /// real events dominate once the wasted flights are gone, so there is
  /// nothing to buy below 10 meV. Widen the margin before going narrower.
  /// exact maximum of the piecewise-linear phaseA table over [e0, e1].
  /// interp() is linear between grid nodes, so the max is attained at an
  /// endpoint or at a grid node strictly inside the interval - no scanning
  /// approximation, no assumption that Gamma_A is monotone in E.
  static T maxRateOver(const Mechanism &m, T e0, T e1) {
    T mx = std::max(interp(m, e0), interp(m, e1));
    for (auto it = std::upper_bound(m.grid.begin(), m.grid.end(), e0);
         it != m.grid.end() && *it < e1; ++it)
      mx = std::max(mx, m.rates[it - m.grid.begin()]);
    return mx;
  }

  /// Builds a PROVABLE upper bound on the scattering rate, per energy slab.
  ///
  /// The old construction sampled the rate at MESH POINTS and padded with a
  /// one-slab neighbour margin. That is not a bound: at runtime the rate is
  ///     sum_m  interp(m, E) * interpPt(m, t, lam)
  /// i.e. a PRODUCT OF TWO INTERPOLANTS, and a product can exceed the max of
  /// its vertex products inside a tetrahedron. On the production package that
  /// leaked 274 breaches at 10 meV slabs, and the engine printed "results are
  /// NOT valid". Widening slabs only raised the bound by luck; it left the
  /// guarantee absent and the width a per-material tuning knob - fatal for a
  /// material-agnostic pipeline, since a new material would silently need a
  /// different constant and g0viol is only a detector, not a guarantee.
  ///
  /// The bound below is exact by construction. Inside tetrahedron t:
  ///   - lam are barycentric, non-negative, summing to 1, so interpPt is a
  ///     CONVEX COMBINATION of the four vertex ptRates and cannot exceed
  ///     their maximum;
  ///   - E is the linear interpolant of the vertex energies, so it lies in
  ///     [Emin(t), Emax(t)] and interp(m, E) <= maxRateOver(m, Emin, Emax).
  /// Hence  rate(k in t) <= sum_m max_v ptRates_m[v] * maxRateOver(m, ...).
  /// Summing per-mechanism maxima over-estimates the max of the sum, which is
  /// still a valid upper bound. Every slab the tet's energy span touches then
  /// takes at least that value, so any particle at energy E inside t reads a
  /// slab value >= its own rate. gamma0At(E) is unchanged.
  ///
  /// slabWidth now controls RESOLUTION ONLY - correctness no longer depends on
  /// it. `safety` (>1) absorbs locateTet's small extrapolation tolerance,
  /// which can push lam marginally outside the simplex.
  void buildGamma0Slabs(T slabWidth = static_cast<T>(0.010),
                        T safety = static_cast<T>(1.25),
                        T eGuard = static_cast<T>(0.005)) {
    gamma0Global = getGamma0(safety);
    const auto &Ept = bs.getBandEnergies(bandIdx);
    T lo = Ept[0], hi = Ept[0];
    for (T e : Ept) { lo = std::min(lo, e); hi = std::max(hi, e); }
    slabLo = lo;
    slabW = slabWidth;
    const std::size_t n =
        static_cast<std::size_t>((hi - lo) / slabW) + 2;
    slabG0.assign(n, T(0));

    const auto &tets = bs.getTetrahedra();
    std::vector<T> mrScratch(mechanisms.size(), T(0)); // hoisted: per-tet alloc
    for (std::size_t t = 0; t < tets.size(); t++) {
      const auto &tv = tets[t];
      T e0 = Ept[tv[0]], e1 = Ept[tv[0]];
      for (int v = 1; v < 4; v++) {
        e0 = std::min(e0, Ept[tv[v]]);
        e1 = std::max(e1, Ept[tv[v]]);
      }
      std::size_t i0 = static_cast<std::size_t>((e0 - slabLo) / slabW);
      std::size_t i1 = static_cast<std::size_t>((e1 - slabLo) / slabW);
      if (i1 >= n) i1 = n - 1;
      for (std::size_t i = i0; i <= i1; i++) {
        // Bound slab i using only the part of the tet that can actually be
        // seen from slab i: a particle counted in slab i has E in slab i AND
        // in [e0,e1], so the energy max is taken over the INTERSECTION. Using
        // the whole tet span instead lets one wide tetrahedron push its
        // high-energy rate onto every low-energy slab it touches - which cost
        // mode 1 a jump from 38% to 70% self-scattering for no physics.
        // eGuard is an ABSOLUTE numerical tolerance (eV), not a per-material
        // tuning knob: locateTet accepts a tetrahedron within a small
        // extrapolation tolerance, so lam can fall marginally outside the
        // simplex and the particle's energy marginally outside [e0,e1].
        // Without it the bound is exact only down to slabW ~ 10 meV; at 2 and
        // 4 meV slabs it leaked 19 and 13 breaches - the bound's logic was
        // sound, the locate slop was not covered.
        const T sLo = std::max(e0, slabLo + static_cast<T>(i) * slabW) - eGuard;
        const T sHi = std::min(e1, slabLo + static_cast<T>(i + 1) * slabW) + eGuard;
        for (std::size_t j = 0; j < mechanisms.size(); j++)
          mrScratch[j] = maxRateOver(mechanisms[j], sLo, sHi);
        // sum_m mr[m] * interpPt(m,t,lam) is LINEAR in lam - the barycentric
        // interpolation of the vertex values sum_m mr[m]*ptRates_m[v] - and a
        // linear function on a simplex is maximised at a vertex. Taking the
        // max of the SUM is strictly tighter than the sum of per-mechanism
        // maxima, which would let every mechanism peak at a different vertex.
        T bound = 0;
        if (kResolved) {
          for (int v = 0; v < 4; v++) {
            T s = 0;
            for (std::size_t j = 0; j < mechanisms.size(); j++)
              s += mrScratch[j] * mechanisms[j].ptRates[tv[v]];
            bound = std::max(bound, s);
          }
        } else {
          for (std::size_t j = 0; j < mechanisms.size(); j++)
            bound += mrScratch[j];
        }
        bound *= safety;
        slabG0[i] = std::max(slabG0[i], bound);
      }
    }
    // FLIGHT-EXCURSION MARGIN - do not remove.
    // The driver picks Gamma0 = gamma0At(E) at the START of a free flight and
    // then accelerates the particle before testing it against the rate at the
    // END energy, so the slab value must bound the rate over the whole energy
    // range the flight traverses, not just the starting slab. The excursion is
    // NOT negligible: near the band edge the rates are small, so tau ~ 1/Gamma0
    // reaches ~1 ps and dE = qF*v*tau ~ 2 meV at 200 V/cm - the same order as a
    // slab. (An earlier estimate of ~1 ueV assumed the peak Gamma0 and was
    // wrong; that mistake is why this margin was briefly deleted, which cost
    // mode 1 ten breaches at 10 meV slabs while mode 0 stayed clean, mode 0
    // having larger Gamma0 and therefore shorter flights.)
    // One slab of headroom on each side, on top of eGuard.
    {
      std::vector<T> pad(slabG0);
      for (std::size_t i = 0; i < n; i++) {
        if (i)         pad[i] = std::max(pad[i], slabG0[i - 1]);
        if (i + 1 < n) pad[i] = std::max(pad[i], slabG0[i + 1]);
      }
      slabG0.swap(pad);
    }
    // a slab no tetrahedron touches can never be sampled; keep it finite
    for (std::size_t i = 0; i < n; i++)
      if (slabG0[i] <= 0) slabG0[i] = gamma0Global;
  }

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
                        std::int64_t &tetHint,
                        std::size_t mech = static_cast<std::size_t>(-1)) const {
    const std::int64_t bin = binOf(energy);
    if (bin < 0 || bins[bin].tets.empty())
      return false;
    const auto &B = bins[bin];
    // Weighted draw of a tet, REPEATED until one actually straddles this
    // exact energy. A bin is 2 meV wide and its list holds every tet that
    // OVERLAPS the bin, so a tet drawn for energy E often does not straddle
    // E itself. The previous code then placed the carrier at the tet
    // CENTROID - off the iso-surface, at the wrong energy - and returned
    // success, silently breaking per-event energy conservation on several
    // percent of all events (33.5k on the elastic test package, 59.6k on
    // the inelastic one). Redrawing keeps the DOS weighting exact
    // conditional on straddling; only a pathological bin now fails, and it
    // fails visibly instead of fabricating a state.
    // DETAILED BALANCE IN k. The out-rate in the k-resolved modes is
    // a(k)*Gamma_A(E), but this sampler weights final states by the DOS alone.
    // W(k->k')p(k) = W(k'->k)p(k') then requires a(k)w(k') = a(k')w(k), which
    // is false in general - so the ensemble relaxes to the WRONG equilibrium,
    // by an amount that depends on the rate model. Measured at zero field:
    // mode 1 <E> = 0.0433 eV vs modes 2/0 0.0446 eV, a 3% split that cannot
    // exist if balance holds (the equilibrium may not depend on the rates).
    // Weighting the final state by a(k')w(k') makes the k-parts symmetric:
    // a(k)a(k')w(k') vs a(k')a(k)w(k). Implemented as rejection ON TOP of the
    // existing DOS draw - accept with probability a_tet/aMax - which needs no
    // per-mechanism tables and leaves the DOS weighting exact.
    const bool aWeight = kResolved && mech < mechanisms.size() &&
                         !mechanisms[mech].ptRates.empty() &&
                         mech < aMaxMech.size() && aMaxMech[mech] > 0;
    std::uniform_real_distribution<T> U(0, B.totalWeight);
    std::uniform_real_distribution<T> U01(0, 1);
    const int maxAttempt = aWeight ? 512 : 64;
    for (int attempt = 0; attempt < maxAttempt; attempt++) {
      T r = U(rng), acc = 0;
      std::size_t pick = B.tets.size() - 1;
      for (std::size_t i = 0; i < B.tets.size(); i++) {
        acc += B.weights[i];
        if (r <= acc) {
          pick = i;
          break;
        }
      }
      // bin entries are PACKED ids (tet << 3 | sub-tet) when subdivision is
      // on; everything that indexes mesh arrays needs the parent id
      const std::int64_t id = B.tets[pick];
      if (aWeight) {
        const std::int64_t par = bs.hasSubdivision() ? (id >> 3) : id;
        const auto &tv = bs.getTetrahedra()[par];
        const auto &pr = mechanisms[mech].ptRates;
        const T aTet = static_cast<T>(0.25) *
                       (pr[tv[0]] + pr[tv[1]] + pr[tv[2]] + pr[tv[3]]);
        aWeightTries.fetch_add(1, std::memory_order_relaxed);
        if (U01(rng) * aMaxMech[mech] > aTet) {
          aWeightRej.fetch_add(1, std::memory_order_relaxed);
          continue;   // rejected: yields the a(k')*w(k') measure
        }
      }
      if (samplePointOnIso(id, energy, rng, kOut)) {
        tetHint = bs.hasSubdivision() ? (id >> 3) : id;
        return true;
      }
    }
    centroidFallbacks.fetch_add(1, std::memory_order_relaxed);
    return false;
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
  ///
  /// `band` selects which band's energies define the iso-surface; SIZE_MAX
  /// means this instance's own band. INTERBAND final states must pass the
  /// DESTINATION band - resolving E' against the source band would test the
  /// wrong tets for straddling and place k' on the wrong surface, so a
  /// cross-band transition would land at an energy that violates
  /// E' = E +/- hbar*omega by the inter-band separation.
  template <class RNG>
  bool samplePointOnIsoFrac(std::int64_t packed, T energy, RNG &rng,
                            Vec3 &fOut,
                            std::size_t band = SIZE_MAX) const {
    if (band == SIZE_MAX) band = bandIdx;
    // iso-surface polygon of level `energy` inside a tet (frac coords).
    // With subdivision the id packs (tet << 3 | sub-tet); the energy is
    // linear within a sub-tet by construction, so the planar polygon and
    // the Blochl weights remain exact.
    const bool sd = bs.hasSubdivision();
    const std::int64_t t = sd ? (packed >> 3) : packed;
    const int sub = sd ? static_cast<int>(packed & 7) : -1;
    const auto &tet = bs.getTetrahedra()[t];
    const auto &pts = bs.getPoints();
    const auto &E = bs.getBandEnergies(band);
    std::array<Vec3, 4> v;
    std::array<T, 4> e;
    for (int i = 0; i < 4; i++) {
      if (sd) {
        const int node = emcNumericBandStructure<T>::SUBTET[sub][i];
        v[i] = bs.nodePos(t, node);
        e[i] = bs.nodeEnergy(t, band, node);
      } else {
        v[i] = pts[tet[i]];
        e[i] = E[tet[i]];
      }
    }
    Vec3 poly[4];
    int np = 0;
    static const int edges[6][2] = {{0, 1}, {0, 2}, {0, 3},
                                    {1, 2}, {1, 3}, {2, 3}};
    for (const auto &ed : edges) {
      T f;
      // sub-tets carry a linear field; whole tets use the consistent solver
      const bool cross =
          sd ? linearCrossing(e[ed[0]], e[ed[1]], energy, f)
             : bs.edgeCrossing(t, band, ed[0], ed[1], energy, f);
      if (cross) {
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
  /// highest destination band this engine can place a carrier in; SIZE_MAX
  /// means unrestricted. Set by the driver from ENGINE_NBANDS.
  std::size_t maxDestBand = SIZE_MAX;
  mutable std::atomic<std::size_t> coulombTries{0}, coulombExhausted{0};
  /// WP3d: one frame per stamped valley image of THIS band - its minimum
  /// (cartesian) and M^{+1/2}, M^{-1/2} built from the stamped principal
  /// masses and axes rotated to that image. Used by the Herring-Vogt
  /// proposal in sampleFinalStateCoulomb.
  struct ValleyFrame {
    Vec3 kmin;
    std::array<std::array<T, 3>, 3> mh, mih;
  };
  std::vector<std::vector<ValleyFrame>> valleysByBand;   // [package band]
  std::vector<bool> isoByBand;

  void buildValleyFrames() {
    valleysByBand.assign(bs.getNrBands(), {});
    isoByBand.assign(bs.getNrBands(), true);
    for (std::size_t b = 0; b < static_cast<std::size_t>(bs.getNrBands()); b++)
      buildValleyFramesFor(b);
  }

  void buildValleyFramesFor(std::size_t band) {
    auto &valleys = valleysByBand[band];
    valleys.clear();
    const auto &ep = bs.edgePar;
    if (!ep.hasStar) return;
    const std::ptrdiff_t ip = ep.find(static_cast<std::int64_t>(band));
    if (ip < 0) return;
    const std::size_t I = static_cast<std::size_t>(ip);
    const double *m = &ep.masses[I * 3];
    if (!(m[0] > 0)) return;
    isoByBand[band] = (m[2] / m[0] < 1.05);
    std::array<Vec3, 3> ax0;
    for (int i = 0; i < 3; i++)
      for (int c = 0; c < 3; c++)
        ax0[i][c] = static_cast<T>(ep.axes[I * 9 + i * 3 + c]);
    Vec3 s0{};
    for (int c = 0; c < 3; c++) s0[c] = static_cast<T>(ep.star[(I * ep.nvMax) * 3 + c]);
    const auto &ops = bs.getSymOps();
    for (std::size_t v = 0; v < ep.nvMax; v++) {
      const double *sv = &ep.star[(I * ep.nvMax + v) * 3];
      if (std::isnan(sv[0])) continue;
      // the frac rotation carrying image 0 to image v (mod G); identity if none
      std::array<std::array<T, 3>, 3> R{{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};
      for (const auto &op : ops) {
        bool hit = true;
        for (int i = 0; i < 3 && hit; i++) {
          T x = op[i][0] * s0[0] + op[i][1] * s0[1] + op[i][2] * s0[2] - static_cast<T>(sv[i]);
          x -= std::round(x);
          if (std::abs(x) > 1e-5) hit = false;
        }
        if (hit) { R = op; break; }
      }
      ValleyFrame f{};
      Vec3 fv{}; for (int c = 0; c < 3; c++) fv[c] = static_cast<T>(sv[c]);
      f.kmin = bs.fracToCart(fv);
      std::array<Vec3, 3> ax;
      for (int i = 0; i < 3; i++) {
        // a cartesian direction transforms as fracToCart(R * cartToFrac(e))
        const Vec3 ef = bs.cartToFrac(ax0[i]);
        Vec3 rf{};
        for (int r = 0; r < 3; r++) rf[r] = R[r][0]*ef[0] + R[r][1]*ef[1] + R[r][2]*ef[2];
        ax[i] = bs.fracToCart(rf);
        const T n = std::sqrt(ax[i][0]*ax[i][0] + ax[i][1]*ax[i][1] + ax[i][2]*ax[i][2]);
        for (int c = 0; c < 3; c++) ax[i][c] /= n;
      }
      for (int r = 0; r < 3; r++)
        for (int c = 0; c < 3; c++) {
          T h = 0, ih = 0;
          for (int i = 0; i < 3; i++) {
            const T sq = std::sqrt(static_cast<T>(m[i]));
            h += sq * ax[i][r] * ax[i][c];
            ih += ax[i][r] * ax[i][c] / sq;
          }
          f.mh[r][c] = h; f.mih[r][c] = ih;
        }
      valleys.push_back(f);
    }
  }

  mutable std::atomic<std::size_t> coulombAnalytic{0}, coulombProposals{0},
      coulombBracketFail{0}, coulombCapped{0}, coulombFallbacks{0};
public:
  std::size_t getCoulombTries() const { return coulombTries.load(); }
  std::size_t getCoulombExhausted() const { return coulombExhausted.load(); }
  std::size_t getCoulombAnalytic() const { return coulombAnalytic.load(); }
  std::size_t getCoulombProposals() const { return coulombProposals.load(); }
  std::size_t getCoulombBracketFail() const { return coulombBracketFail.load(); }
  std::size_t getCoulombCapped() const { return coulombCapped.load(); }
  std::size_t getCoulombFallbacks() const { return coulombFallbacks.load(); }
private:
  mutable std::atomic<std::size_t> g2Rejected{0};
  /// true when the package carried per-source-band g2 rows (v0.4). False
  /// means every band is reading the band-0 rows, which is only correct for
  /// a single-band simulation - the driver reports it so an interband run
  /// cannot silently proceed on v0.3 tables.
  bool g2PerSourceBand = false;
  std::vector<Mechanism> mechanisms;
  bool kResolved = false;
  bool g2Loaded = false;
  // atomic: sampling is const and may run from several ensemble threads
  std::vector<T> aMaxMech;   // per-mechanism max of a(k), for the
                             // rejection bound in sampleFinalState
  mutable std::atomic<std::size_t> aWeightTries{0}, aWeightRej{0};
  mutable std::atomic<std::size_t> g2Hits{0}, g2Fallbacks{0};
  mutable std::atomic<std::size_t> g2Walk{0}, g2WalkLen{0};
  mutable std::atomic<std::size_t> centroidFallbacks{0};
  mutable std::atomic<std::size_t> g0Viol{0};
  std::vector<T> slabG0;
  T slabLo = 0, slabW = 0.025, gamma0Global = 0;
  std::size_t nrewrite_ = 0;
  T relSum_ = 0, relMax_ = 0;
  std::size_t relN_ = 0;
  std::vector<T> dosGrid, dosVals;   // package linear-tetrahedron DOS
  mutable std::vector<T> thermCum;   // cumulative DOS*Boltzmann over bins
  mutable T thermKT = -1;
  void buildThermalCdf(T kT) const {
    thermCum.assign(bins.size(), T(0));
    T acc = 0;
    for (std::size_t i = 0; i < bins.size(); i++) {
      const T Ec = binLo + (static_cast<T>(i) + T(0.5)) * binW;
      acc += bins[i].dos * std::exp(-(Ec - binLo) / kT);
      thermCum[i] = acc;
    }
    thermKT = kT;
  }

  /// loads the package's own DOS (linear tetrahedron, absolute energies)
  void loadDos(const std::string &file, std::size_t band) {
    hid_t f = H5Fopen(file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (f < 0)
      return;
    auto rd = [&](const char *name, std::vector<T> &out) {
      hid_t d = H5Dopen2(f, name, H5P_DEFAULT);
      if (d < 0)
        return false;
      hid_t sp = H5Dget_space(d);
      hsize_t dims[2] = {0, 0};
      const int nd = H5Sget_simple_extent_dims(sp, dims, nullptr);
      const hsize_t n = nd == 2 ? dims[0] * dims[1] : dims[0];
      std::vector<double> buf(n);
      H5Dread(d, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf.data());
      H5Sclose(sp);
      H5Dclose(d);
      if (nd == 2 && dims[0] > band)       // per_band: take this band's row
        out.assign(buf.begin() + band * dims[1],
                   buf.begin() + (band + 1) * dims[1]);
      else
        out.assign(buf.begin(), buf.end());
      return true;
    };
    std::vector<T> g, v;
    if (rd("/bands/electron/dos/energy_grid", g) &&
        rd("/bands/electron/dos/total", v) && g.size() == v.size()) {
      dosGrid = g;
      dosVals = v;
    }
    H5Fclose(f);
  }

  /// loads /scattering/g2bins: CSR rows per mechanism by name
  ///
  /// v0.4 groups the rows by SOURCE band under `b{n}/`, so each per-band
  /// instance reads its own final states. v0.3 stored one flat set built from
  /// band-0 sources; when `b{n}/` is absent we read that, which reproduces
  /// the old behaviour exactly (every band shared the band-0 rows).
  bool loadG2Tables(const std::string &file) {
    hid_t f = H5Fopen(file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (f < 0)
      return false;
    hid_t root = H5Gopen2(f, "/scattering/g2bins", H5P_DEFAULT);
    if (root < 0) {
      H5Fclose(f);
      return false;
    }
    const std::string sub = "b" + std::to_string(bandIdx);
    hid_t grp = root;
    bool perBand = false;
    if (H5Lexists(root, sub.c_str(), H5P_DEFAULT) > 0) {
      hid_t bg = H5Gopen2(root, sub.c_str(), H5P_DEFAULT);
      if (bg >= 0) { grp = bg; perBand = true; }
    }
    g2PerSourceBand = perBand;
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
      // optional (spec v0.4); absent on every earlier package
      if (H5Lexists(mg, "dest_band", H5P_DEFAULT) > 0)
        read1("dest_band", m.g2Band, H5T_NATIVE_INT64);
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
    if (perBand)
      H5Gclose(grp);
    H5Gclose(root);
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
    T totalWeight = 0;   // sum of `weights`: the tet-SAMPLING normaliser
    T dos = 0;           // the density the DOS consumers read - equals
                         // totalWeight unless the band-edge patch replaced it
  };
  T binLo = 0, binW = 0;
  std::vector<Bin> bins;

  std::int64_t binOf(T energy) const {
    const std::int64_t b = static_cast<std::int64_t>((energy - binLo) / binW);
    return (b >= 0 && b < (std::int64_t)bins.size()) ? b : -1;
  }

  /// "eph_m3_abs" / "eph_m3_em" -> "eph_m3": pairs a phonon mode's two channels
  static std::string partnerName(const std::string &n) {
    const std::size_t p = n.rfind('_');
    return p == std::string::npos ? n : n.substr(0, p);
  }

  /// linear interpolation of `v` on `grid`, skipping the -1 holes left where
  /// a grid point had no usable DOS
  static T interpArr(const std::vector<T> &grid, const std::vector<T> &v, T e) {
    if (e <= grid.front() || e >= grid.back())
      return T(-1);
    const std::size_t k =
        std::lower_bound(grid.begin(), grid.end(), e) - grid.begin();
    if (k == 0 || k >= grid.size())
      return T(-1);
    if (v[k] < 0 || v[k - 1] < 0)
      return v[k] >= 0 ? v[k] : v[k - 1];
    const T f = (e - grid[k - 1]) / (grid[k] - grid[k - 1]);
    return v[k - 1] + f * (v[k] - v[k - 1]);
  }

  static std::size_t nearestGrid(const Mechanism &m, T e) {
    std::size_t best = 0;
    T bd = std::fabs(m.grid[0] - e);
    for (std::size_t i = 1; i < m.grid.size(); i++) {
      const T d = std::fabs(m.grid[i] - e);
      if (d < bd) { bd = d; best = i; }
    }
    return best;
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

  static bool linearCrossing(T ea, T eb, T energy, T &f) {
    if ((ea - energy) * (eb - energy) >= 0)
      return false;
    f = (energy - ea) / (eb - ea);
    return true;
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
      if (H5Aexists(mg, "beta_1_per_m") > 0) {
        hid_t ab = H5Aopen(mg, "beta_1_per_m", H5P_DEFAULT);
        double beta = 0;
        H5Aread(ab, H5T_NATIVE_DOUBLE, &beta);
        H5Aclose(ab);
        m.coulombBeta2 = static_cast<T>(beta * beta);
        if (H5Aexists(mg, "dest_band") > 0) {
          hid_t ad = H5Aopen(mg, "dest_band", H5P_DEFAULT);
          long long db = -1;
          H5Aread(ad, H5T_NATIVE_LLONG, &db);
          H5Aclose(ad);
          if (db >= 0) m.coulombDestBand = static_cast<std::size_t>(db);
        }
      }
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
  /// fractional-space volume of tet t (constant Jacobian to cartesian)
  T vfrac(std::int64_t t) const {
    const auto &pts = bs.getPoints();
    const auto &tet = bs.getTetrahedra()[t];
    const auto &a = pts[tet[0]];
    T e1[3], e2[3], e3[3];
    for (int c = 0; c < 3; c++) {
      e1[c] = pts[tet[1]][c] - a[c];
      e2[c] = pts[tet[2]][c] - a[c];
      e3[c] = pts[tet[3]][c] - a[c];
    }
    return std::abs(e1[0] * (e2[1] * e3[2] - e2[2] * e3[1]) -
                    e1[1] * (e2[0] * e3[2] - e2[2] * e3[0]) +
                    e1[2] * (e2[0] * e3[1] - e2[1] * e3[0])) / 6;
  }

  /// distribute a (sub-)tet of volume `vol` spanning [lo,hi] over the bins
  /// Fraction of a linear tetrahedron's VOLUME with energy below E.
  /// Exact cubic CDF for vertex energies sorted e0<=e1<=e2<=e3.
  static T volFrac(T e0, T e1, T e2, T e3, T E) {
    const T tiny = 1e-12;
    if (E <= e0)
      return 0;
    if (E >= e3)
      return 1;
    if (E < e1) {
      const T d = std::max((e1 - e0) * (e2 - e0) * (e3 - e0), tiny);
      return (E - e0) * (E - e0) * (E - e0) / d;
    }
    if (E < e2) {
      const T d20 = std::max(e2 - e0, tiny), d30 = std::max(e3 - e0, tiny);
      const T d21 = std::max(e2 - e1, tiny), d31 = std::max(e3 - e1, tiny);
      const T x = E - e1, e10 = e1 - e0;
      // Lehmann-Taut middle interval, integrated
      return (e10 * e10 + 3 * e10 * x + 3 * x * x
              - (d20 + d31) * x * x * x / (d21 * d31)) / (d20 * d30);
    }
    const T d = std::max((e3 - e0) * (e3 - e1) * (e3 - e2), tiny);
    return 1 - (e3 - E) * (e3 - E) * (e3 - E) / d;
  }

  /// Bin weights from the EXACT energy distribution of the tet's volume.
  ///
  /// The previous form spread `vol` UNIFORMLY over [lo,hi]
  /// (w = vol*overlap/(hi-lo)). A linear tetrahedron's volume is not
  /// uniformly distributed in energy - it follows the linear-tetrahedron
  /// law - so the sampler's measure was biased, and being a bias in the
  /// SHAPE it did not improve as the bins narrowed: the implied equilibrium
  /// sat at 0.0431 eV for every bin width from 2 meV down to 0.1 meV, while
  /// the mesh's own tetrahedron DOS gives 0.0396 eV.
  void addToBinsExact(std::int64_t id, T e0, T e1, T e2, T e3, T vol) {
    const T maxE = binLo + binW * (bins.size() - 1) + binW;
    const std::int64_t b0 = std::max<std::int64_t>(0, binOf(e0));
    const std::int64_t bh = binOf(std::min(e3, maxE));
    const std::int64_t b1 =
        std::min<std::int64_t>(bins.size() - 1, bh < 0 ? bins.size() - 1 : bh);
    for (std::int64_t b = b0; b <= b1; b++) {
      const T bl = binLo + b * binW, bhh = bl + binW;
      const T f = volFrac(e0, e1, e2, e3, std::min(bhh, e3)) -
                  volFrac(e0, e1, e2, e3, std::max(bl, e0));
      if (f <= 0)
        continue;
      const T w = vol * f;
      bins[b].tets.push_back(id);
      bins[b].weights.push_back(w);
      bins[b].totalWeight += w;
    }
  }

  void addToBins(std::int64_t id, T lo, T hi, T vol) {
    const T maxE = binLo + binW * (bins.size() - 1) + binW;
    const std::int64_t b0 = std::max<std::int64_t>(0, binOf(lo));
    const std::int64_t bh = binOf(std::min(hi, maxE));
    const std::int64_t b1 =
        std::min<std::int64_t>(bins.size() - 1, bh < 0 ? bins.size() - 1 : bh);
    for (std::int64_t b = b0; b <= b1; b++) {
      const T bl = binLo + b * binW, bhh = bl + binW;
      const T ov = std::min(hi, bhh) - std::max(lo, bl);
      if (ov <= 0)
        continue;
      const T w = vol * ov / (hi - lo);
      bins[b].tets.push_back(id);
      bins[b].weights.push_back(w);
      bins[b].totalWeight += w;
    }
  }

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
      if (bs.hasSubdivision()) {           // bin each sub-tet separately
        // A sub-tet's energy field is LINEAR by construction (same reason
        // samplePointOnIsoFrac can use a planar polygon there), so the exact
        // linear-tetrahedron law applies and `addToBins` - which spreads the
        // volume UNIFORMLY over [lo,hi] - is simply the wrong measure. That
        // uniform spread is the bias `addToBinsExact` was written to remove,
        // and it was still in force on every subdivided run.
        for (int sb = 0; sb < 8; sb++) {
          T slo, shi;
          bs.subRange(t, bandIdx, sb, slo, shi);
          if (shi <= binLo || slo >= maxEnergy || shi - slo < 1e-12)
            continue;
          T sv[4];
          for (int i = 0; i < 4; i++)
            sv[i] = bs.nodeEnergy(
                t, bandIdx, emcNumericBandStructure<T>::SUBTET[sb][i]);
          std::sort(sv, sv + 4);
          addToBinsExact((t << 3) | sb, sv[0], sv[1], sv[2], sv[3],
                         vfrac(t) / 8);
        }
        continue;
      }
      if (bs.hasQuadraticEnergy()) {
        for (const auto &pr : probe) {
          const Vec3 lam{pr[0], pr[1], pr[2]};
          const T ev = bs.tetEnergy(t, bandIdx, lam);
          lo = std::min(lo, ev);
          hi = std::max(hi, ev);
        }
      }
      // Without subdivision the quadratic interpolant is neither linear nor
      // uniform in energy, so NEITHER measure below is exact for it and it
      // falls through to the uniform ramp. SUBDIV=1 is what makes the
      // quadratic path exact - see the branch above. Measured cost of the
      // uniform ramp on this path: modes 2 and 0, which must agree at
      // equilibrium, land at <E> = 0.0431/0.0433 instead of 0.0445/0.0446.
      if (!bs.hasQuadraticEnergy() && hi > binLo && lo < maxEnergy
          && hi - lo >= 1e-12) {
        T sv[4] = {E[tets[t][0]], E[tets[t][1]], E[tets[t][2]], E[tets[t][3]]};
        std::sort(sv, sv + 4);
        addToBinsExact(t, sv[0], sv[1], sv[2], sv[3], vfrac(t));
        continue;
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
    for (auto &bn : bins)
      bn.dos = bn.totalWeight;
    patchEdgeDos();
  }

  /*! \brief Replace the linear-tetrahedron DOS below E_sw by an analytic
   * Kane form fitted to the mesh's own DOS above E_sw.
   *
   * WHY. Bands are interpolated LINEARLY inside each tetrahedron, and a linear
   * function has no interior minimum: the tet containing the true band
   * minimum has all four vertices ABOVE it, so the interpolated DOS is
   * deficient - zero below the lowest vertex, then low by O((dk)^2 x
   * curvature). Measured on Si REFINE=6: 0.59 of parabolic at 10 meV, 0.77 at
   * 20; REFINE=12 reaches 0.87 / 0.90. The consequence is not cosmetic: the
   * equilibrium sits HOT at any field (<E> = 0.0446 eV at 25 V/cm against
   * 1.5kT = 0.0388), every ensemble average inherits it, and the fix by mesh
   * alone costs a 10-33 h table rebuild per attempt because the scattering
   * tables are indexed by mesh state.
   *
   * This is the field's standard remedy - analytic band edge, full band above
   * (Fischetti-Laux lineage; Jacoboni-Reggiani) - applied to the DOS only. It
   * is the region where a mass tensor plus Kane non-parabolicity is exact to
   * a few percent and a linear mesh is worst, so it is the right region to be
   * analytic in. The form is A * sqrt(x(1+ax)) * (1+2ax), the Kane DOS with x
   * measured from the band minimum, and BOTH A and alpha are fitted to this
   * instance's own bins on [E_sw, E_hi]: A absorbs the mesh's volume units
   * and the valley count, alpha the curvature. Nothing is imported - the gate
   * docstring records that silicon's alpha is NOT 1/Eg, so a formula would be
   * a fit in disguise.
   *
   * WHAT IT DOES NOT DO. Final states are still placed on the linear-tet
   * iso-surface; bins with no straddling tet stay empty (the sampler's
   * tets.empty() guard makes that safe) and the analytic weight that lands in
   * them is reported as `unsampleable`. Placing k' analytically below E_sw is
   * the next stage.
   *
   * VALIDATION IS EXTERNAL OR IT IS NOTHING. Enforcing balance against ANY
   * density makes the stationary distribution D(E)exp(-E/kT), so the
   * bin-implied equilibrium landing on 1.5kT proves nothing - see
   * enforceDetailedBalance. The checks that count: the parabolic test
   * packages, whose exact DOS is sqrt(E), and mesh-independence.
   *
   * EDGE_PATCH=0 disables (A/B). EDGE_SW [eV] switch energy (default 0.05),
   * EDGE_FIT_HI [eV] top of the fit window (default 0.15).
   */
  /*! \brief Carry the corrected edge DOS into the RATES - the only place mode 0
   * reads it.
   *
   * The equilibrium is set by the dynamics, and the dynamics read
   * Gamma_m(E) x a(k) (phaseA grid x phaseB anisotropy). Both were integrated
   * by the converter with the same linear-tetrahedron DOS at the FINAL state,
   * so transitions INTO the band edge are under-counted by exactly the deficit
   * the DOS patch measures: emission into an edge bin is too rare, the edge is
   * under-populated, the ensemble sits hot. Rescaling Gamma_m(E) by
   * r(E + dE_m) = dos/totalWeight of the FINAL bin fixes that at its source:
   * emission into the edge gets the analytic DOS, absorption out of it (final
   * states higher up, on a trusted mesh) is untouched, and the em/abs ratio at
   * an edge pair becomes D_an(E - hw)/D_mesh(E) - the correct one. The
   * anisotropy factor a(k) = Gamma_B/Gamma_A(E) is a ratio and needs no
   * change. The g2 rows need none either: E' fixes the bin, and all of a row's
   * tets share it, so r cancels within the row.
   *
   * Measured before this existed (2026-09-02): the DOS-only patch moved the
   * bin-implied equilibrium 0.0439 -> 0.0408 on Si REFINE=6 and the ENSEMBLE
   * not at all (0.0450 -> 0.0450), with or without BALANCE=1 - which
   * rebalances the phaseA tables mode 0 does not use, and stays off (PLAN).
   *
   * Interband rows: the final state of a cross-band mechanism sits on ANOTHER
   * band's edge; this uses the source instance's bins, i.e. it assumes the
   * intraband part dominates Gamma_m(E). EDGE_RATES=0 disables (A/B against
   * DOS-only).
   */
  void rescaleRatesToEdgeDos() {
    const char *er = std::getenv("EDGE_RATES");
    if (er && std::atoi(er) == 0) {
      std::printf("# band %zu: edge-DOS rate rescale OFF (EDGE_RATES=0)\n", bandIdx);
      return;
    }
    std::size_t nTouched = 0;
    double rMin = 1e300, rMax = 0;
    for (auto &m : mechanisms) {
      for (std::size_t i = 0; i < m.grid.size() && i < m.rates.size(); i++) {
        const std::int64_t b = binOf(m.grid[i] + m.deltaE);
        if (b < 0) continue;
        const auto &bn = bins[b];
        if (bn.totalWeight <= 0 || bn.dos <= 0 || bn.dos == bn.totalWeight)
          continue;
        const double r = static_cast<double>(bn.dos / bn.totalWeight);
        m.rates[i] *= static_cast<T>(r);
        nTouched++;
        rMin = std::min(rMin, r); rMax = std::max(rMax, r);
      }
    }
    if (nTouched)
      std::printf("# band %zu: edge-DOS rate rescale: %zu grid entries across %zu "
                  "mechanism(s) scaled by r(E_final) in [%.3f, %.3f]\n",
                  bandIdx, nTouched, mechanisms.size(), rMin, rMax);
    else
      std::printf("# band %zu: edge-DOS rate rescale: nothing to do\n", bandIdx);
  }

  void patchEdgeDos() {
    // OFF BY DEFAULT. What is actually measured (2026-09-02):
    //   pv_A1 N72, 20000 x 8, exact mu = e tau/m* = 674.9:
    //     unpatched 665.0 +/- 19.5, DOS-only 683.6 +/- 16.0 - both within 1
    //     sigma of exact; <E> 0.0397 -> 0.0390 toward the exact 0.0388. The
    //     DOS-only patch is correct and harmless there. (Earlier 8000 x 4 runs
    //     read 767/817 and were replica noise - the 4-replica SE is ~2x
    //     optimistic, see the driver.)
    //   Si REFINE=6, 8000 x 4: DOS-only leaves the ensemble untouched (0.0450,
    //     mode 0 never reads bins[].dos); DOS+RATES cools it 0.0450 -> 0.0407
    //     as designed but moves mu 1348 -> 1586 with rate factors of 9-35x in
    //     the near-empty lowest bins - multiplying the rate INTO a bin whose
    //     one or two tets place k' on the chord's iso-surface (the linear
    //     interpolant overestimates E by 1.5-1.7x below 10 meV on this mesh)
    //     amplifies the placement error. The rate rescale must not run
    //     without analytic placement.
    // Vertex velocities are EXACT (0.00% vs the tensor on pv); the interpolant
    // errors are in E(k) and hence in placement and binning, not in v.
    // EDGE_PATCH=1 opts in (DOS-only unless EDGE_RATES=1) for development of
    // the full hybrid (analytic E, placement, DOS below E_sw).
    const char *ep = std::getenv("EDGE_PATCH");
    if (!ep || std::atoi(ep) == 0) {
      if (ep)
        std::printf("# band %zu: band-edge DOS patch OFF (EDGE_PATCH=0)\n",
                    bandIdx);
      return;
    }
    const double eSw = std::getenv("EDGE_SW") ? std::atof(std::getenv("EDGE_SW"))
                                              : 0.05;
    const double eHi = std::getenv("EDGE_FIT_HI")
                           ? std::atof(std::getenv("EDGE_FIT_HI")) : 0.15;
    auto kane = [](double x, double a) {
      return std::sqrt(x * (1 + a * x)) * (1 + 2 * a * x);
    };
    // ---- MASS-ANCHORED FORM (preferred) ---------------------------------
    // The converter stamps m_d, nv, e0 per band from a quadratic fit to E(k)
    // at MESH POINTS - exact Wannier energies, so immune to the very
    // linear-tetrahedron deficit being corrected. That pins the AMPLITUDE
    // absolutely: the fractional BZ volume below x is
    //     nv * (4 pi / 3) * k(x)^3 / |det B|,   k(x) = sqrt(2 m_d m_e x e)/hbar
    // and the bins are in fractional volume (vfrac), so no unit fit is
    // needed. Only the Kane alpha is fitted, on [E_sw, E_hi]; a one-parameter
    // fit with a pinned amplitude survives a +/-15% wiggly window that a
    // two-parameter shape fit did not (it extrapolated alpha = 2.0 on Si
    // REFINE=6 and REDUCED the edge DOS - see the shape-fit fallback below).
    const std::ptrdiff_t ip = bs.edgePar.has
                                  ? bs.edgePar.find(static_cast<std::int64_t>(bandIdx))
                                  : -1;
    if (ip >= 0) {
      const double mD = bs.edgePar.mD[ip];
      const double nv = static_cast<double>(bs.edgePar.nv[ip]);
      // |det B| from the band structure's own frac->cart map
      const Vec3 b1 = bs.fracToCart(Vec3{1, 0, 0}), b2 = bs.fracToCart(Vec3{0, 1, 0}),
                 b3 = bs.fracToCart(Vec3{0, 0, 1});
      const double detB = std::abs(
          b1[0] * (b2[1] * b3[2] - b2[2] * b3[1]) -
          b1[1] * (b2[0] * b3[2] - b2[2] * b3[0]) +
          b1[2] * (b2[0] * b3[1] - b2[1] * b3[0]));
      constexpr double ME = 9.1093837015e-31, QE = 1.602176634e-19,
                       HBAR = 1.054571817e-34;
      // fractional volume of the Kane band below x (from THIS band's minimum)
      auto volBelow = [&](double x, double al) {
        if (x <= 0) return 0.0;
        const double g = x * (1 + al * x);                 // Kane: E(1+aE)
        const double k = std::sqrt(2 * mD * ME * g * QE) / HBAR;
        return nv * (4.0 * M_PI / 3.0) * k * k * k / detB;
      };
      // alpha only: least-squares on the trusted window with A pinned
      std::vector<std::size_t> win;
      for (std::size_t i = 0; i < bins.size(); i++) {
        const double x = (static_cast<double>(i) + 0.5) * binW;
        if (x >= eSw && x <= eHi && bins[i].totalWeight > 0) win.push_back(i);
      }
      double bestAl = 0, bestRes = std::numeric_limits<double>::max();
      double ymean = 0;
      for (auto i : win) ymean += bins[i].totalWeight;
      ymean = win.empty() ? 1 : ymean / win.size();
      if (win.size() >= 8) {
        for (int ia = 0; ia <= 80; ia++) {
          const double al = 0.05 * ia;
          double res = 0;
          for (auto i : win) {
            const double xl = i * binW, xh = xl + binW;
            const double d = bins[i].totalWeight - (volBelow(xh, al) - volBelow(xl, al));
            res += d * d;
          }
          if (res < bestRes) { bestRes = res; bestAl = al; }
        }
      }
      const double rms = win.size() >= 8 ? std::sqrt(bestRes / win.size()) / ymean : -1;
      // HEALTH CHECK, not a fit: mesh/analytic on the trusted window. Far
      // from 1 means a wrong mass stamp, a wrong valley count, or a units
      // slip - all of which this line makes visible.
      double meshWin = 0, anWin = 0;
      for (auto i : win) {
        meshWin += bins[i].totalWeight;
        anWin += volBelow((i + 1) * binW, bestAl) - volBelow(i * binW, bestAl);
      }
      // GUARD. The analytic form is only usable if it agrees with the mesh
      // where the mesh is trusted. A ratio far from 1 means the stamp does
      // not describe this edge: Si conduction band 1 (a 137 meV feature that
      // is not a parabolic valley) read 0.446, and the sorted branches of the
      // warped degenerate valence manifold read 0.447 / 0.570 - a per-branch
      // parabola over-counts the heavy branch 2x. Patching on those numbers
      // over-corrected the hole equilibrium to 0.0314 eV. Leave such a band
      // on the mesh DOS and say so.
      const double health = anWin > 0 ? meshWin / anWin : 0.0;
      if (health < 0.8 || health > 1.25) {
        std::printf("# band %zu: band-edge DOS patch NOT APPLIED - the stamped "
                    "parabola (m_d = %.4f m_e, nv = %.0f) disagrees with the "
                    "mesh on [%.0f, %.0f] meV: mesh/analytic = %.3f. This "
                    "edge is not a parabolic valley (degenerate or warped "
                    "manifold, or not a minimum); it keeps the mesh DOS.\n",
                    bandIdx, mD, nv, eSw * 1e3, eHi * 1e3, health);
        return;
      }
      double meshBelow = 0, anBelow = 0, unsampleable = 0;
      std::size_t nPatched = 0;
      for (std::size_t i = 0; i < bins.size(); i++) {
        const double xl = i * binW, xh = xl + binW;
        if (xl + 0.5 * binW >= eSw) break;
        const double an = volBelow(xh, bestAl) - volBelow(xl, bestAl);
        meshBelow += bins[i].totalWeight;
        anBelow += an;
        if (bins[i].tets.empty()) { unsampleable += an; continue; }
        bins[i].dos = static_cast<T>(an);
        nPatched++;
      }
      std::printf("# band %zu: band-edge DOS patch (MASS-ANCHORED, m_d = %.4f m_e, "
                  "nv = %.0f) below %.0f meV: Kane alpha = %.2f 1/eV (rms %.1f%% on "
                  "[%.0f, %.0f] meV), mesh/analytic on that window = %.3f; %zu "
                  "bins patched, mesh held %.1f%% of the analytic weight below, "
                  "%.2f%% unsampleable\n",
                  bandIdx, mD, nv, eSw * 1e3, bestAl, 100 * rms, eSw * 1e3,
                  eHi * 1e3, anWin > 0 ? meshWin / anWin : 0.0, nPatched,
                  anBelow > 0 ? 100 * meshBelow / anBelow : 0.0,
                  anBelow > 0 ? 100 * unsampleable / anBelow : 0.0);
      rescaleRatesToEdgeDos();
      return;
    }
    // ---- SHAPE-FIT FALLBACK (no stamps) ---------------------------------
    // Fits BOTH A and alpha to the mesh DOS on [E_sw, E_hi]. Proven on the
    // smooth single-valley parabolic packages (N24: <E> 0.0442 -> 0.0394 vs
    // exact 0.0388); NOT trustworthy on a real material's coarse mesh.
    std::vector<double> xs, ys;
    for (std::size_t i = 0; i < bins.size(); i++) {
      const double x = (static_cast<double>(i) + 0.5) * binW;
      if (x >= eSw && x <= eHi && bins[i].totalWeight > 0) {
        xs.push_back(x);
        ys.push_back(static_cast<double>(bins[i].totalWeight));
      }
    }
    if (xs.size() < 8) {
      std::printf("# band %zu: band-edge DOS patch SKIPPED - only %zu bins in "
                  "[%.3f, %.3f] eV\n", bandIdx, xs.size(), eSw, eHi);
      return;
    }
    // alpha by scan, A by linear least squares at each alpha
    double bestA = 0, bestAl = 0, bestRes = std::numeric_limits<double>::max();
    for (int ia = 0; ia <= 80; ia++) {
      const double al = 0.05 * ia;                      // 0 .. 4 1/eV
      double sgg = 0, syg = 0;
      for (std::size_t i = 0; i < xs.size(); i++) {
        const double g = kane(xs[i], al);
        sgg += g * g; syg += ys[i] * g;
      }
      const double A = sgg > 0 ? syg / sgg : 0;
      double res = 0;
      for (std::size_t i = 0; i < xs.size(); i++) {
        const double d = ys[i] - A * kane(xs[i], al);
        res += d * d;
      }
      if (res < bestRes) { bestRes = res; bestA = A; bestAl = al; }
    }
    double ymean = 0;
    for (double y : ys) ymean += y;
    ymean /= ys.size();
    const double rms = std::sqrt(bestRes / xs.size()) / ymean;
    // apply below E_sw
    double meshBelow = 0, anBelow = 0, unsampleable = 0;
    std::size_t nPatched = 0;
    for (std::size_t i = 0; i < bins.size(); i++) {
      const double x = (static_cast<double>(i) + 0.5) * binW;
      if (x >= eSw) break;
      const double an = bestA * kane(x, bestAl);
      meshBelow += bins[i].totalWeight;
      anBelow += an;
      if (bins[i].tets.empty()) { unsampleable += an; continue; }
      bins[i].dos = static_cast<T>(an);
      nPatched++;
    }
    std::printf("# band %zu: band-edge DOS patch (SHAPE-FIT fallback, no stamps) below %.0f meV: Kane alpha = "
                "%.2f 1/eV, fit rms %.1f%% on [%.0f, %.0f] meV; %zu bins "
                "patched, mesh held %.1f%% of the analytic weight there, "
                "%.2f%% unsampleable (bins with no tet)\n",
                bandIdx, eSw * 1e3, bestAl, 100 * rms, eSw * 1e3, eHi * 1e3,
                nPatched, anBelow > 0 ? 100 * meshBelow / anBelow : 0.0,
                anBelow > 0 ? 100 * unsampleable / anBelow : 0.0);
    rescaleRatesToEdgeDos();
  }
};

#endif // EMC_FULL_BAND_SCATTERING_HPP
