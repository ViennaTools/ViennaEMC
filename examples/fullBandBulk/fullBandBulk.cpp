/*! fullBandBulk - full-band ensemble MC for bulk transport (tasks 2.2-2.4).
 *
 * Reads a .matpkg material package (bands on a tetrahedral BZ mesh + Phase-A
 * scattering tables) and runs an ensemble of electrons under a homogeneous
 * field: event-driven free flight with self-scattering (constant Gamma0),
 * k-update dk/dt = -eF/hbar, table-driven scattering with DOS-weighted
 * final-state selection on constant-energy surfaces.
 *
 * Output: time-averaged drift velocity, mobility and mean energy per field.
 */
#include <FullBand/emcFullBandScattering.hpp>
#include <FullBand/emcNumericBandStructure.hpp>

#include <atomic>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <memory>
#include <random>
#include <string>
#include <thread>
#include <vector>

/// true when the package declares carrier="hole" on /bands/electron.
static bool packageIsHole(const std::string &pkg) {
  hid_t f = H5Fopen(pkg.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if (f < 0) return false;
  bool hole = false;
  hid_t g = H5Gopen2(f, "/bands/electron", H5P_DEFAULT);
  if (g >= 0) {
    if (H5Aexists(g, "carrier") > 0) {
      hid_t a = H5Aopen(g, "carrier", H5P_DEFAULT);
      hid_t t = H5Aget_type(a);
      char buf[64] = {0};
      if (H5Tis_variable_str(t)) {
        char *p = nullptr;
        if (H5Aread(a, t, &p) >= 0 && p) { std::snprintf(buf, sizeof buf, "%s", p); free(p); }
      } else {
        H5Aread(a, t, buf);
      }
      hole = (std::string(buf).find("hole") != std::string::npos);
      H5Tclose(t); H5Aclose(a);
    }
    H5Gclose(g);
  }
  H5Fclose(f);
  return hole;
}

int main(int argc, char **argv) {
  using NBS = emcNumericBandStructure<double>;
  using FBS = emcFullBandScattering<double>;
  constexpr double HBAR = 1.054571817e-34, QE = 1.602176634e-19,
                   KB = 8.617333262e-5; // eV/K
  const std::string pkg =
      argc > 1 ? argv[1] : "/home/filipov/Software/FullBandMC/matforge/Si.matpkg";
  const double T = 300.0;               // K
  // argv[5] = particles per replica, argv[6] = independent replicas.
  // Replicas give an honest standard error: single-run mu noise is +-5-6%
  // at the default 2000 x 20 ps, which hides every few-percent effect.
  const int nPart = argc > 5 ? std::stoi(argv[5]) : 2000;
  const int nRep = argc > 6 ? std::stoi(argv[6]) : 1;
  const double tTotal = 20e-12, tTransient = 5e-12;

  NBS bs(pkg);
  // QUADE=0 disables the quadratic energy interpolation (A/B testing)
  if (const char *q = std::getenv("QUADE"))
    bs.setQuadraticEnergy(std::atoi(q) != 0);
  // SUBDIV=1: 1-to-8 tet subdivision of the ENERGY field with quadratic
  // midpoints (4x lower interpolation error, all machinery stays linear)
  if (const char *d = std::getenv("SUBDIV"))
    bs.setSubdivision(std::atoi(d) != 0);
  const double cbm = bs.getBandMinimum(0);
  // CARRIER: a hole package stores E = VBM - E_electron, so its carriers drift
  // WITH the field. Everything else in the engine is carrier-agnostic because
  // the flipped axis still has a MINIMUM; only the drift sign differs.
  const double carrierSign = packageIsHole(pkg) ? -1.0 : 1.0;
  if (carrierSign < 0)
    std::printf("# HOLE package: carriers drift WITH the field; mobility sign "
                "flipped accordingly\n");
  const int mode = argc > 3 ? std::stoi(argv[3]) : 0;  // 0 auto, 1 phaseA, 2 no-g2
  // BINW: final-state energy bin width [eV]. The inelastic step hbar*omega
  // must be resolved by these bins; when it is comparable to the bin width
  // the engine cannot place the final state at the right energy.
  const double binW = std::getenv("BINW") ? std::atof(std::getenv("BINW")) : 2e-3;
  // WP3-INTERBAND, engine half. One scattering instance PER BAND rather than
  // refactoring the class internals: bins, Gamma0 slabs and the DOS are all
  // per-band, and emcFullBandScattering already handles "a band" correctly -
  // so N of them is far less invasive than making one of them multi-band.
  // ENGINE_NBANDS=1 (default) builds exactly one and costs nothing, so normal
  // runs are unchanged; =2 builds the band-1 instance too, which is what lets
  // a particle in band 1 be propagated and sampled at all.
  // Whether this is PHYSICS or just machinery depends on the package. A v0.4
  // package carries g2 rows grouped by SOURCE band and tagged with a
  // DESTINATION band, and then each instance scatters with its own final
  // states. A v0.3 package has one flat set of band-0 rows, every instance
  // reads it, and a band-1 particle is propagated on band-0 physics - useful
  // for exercising the machinery, meaningless as a transport number. The
  // engine cannot tell these apart silently, so it says which it got.
  const int engineNBands =
      std::getenv("ENGINE_NBANDS") ? std::atoi(std::getenv("ENGINE_NBANDS")) : 1;
  // Interband bookkeeping. Time-weighted occupancy is the number that matters
  // - an event count says transfers happen, occupancy says whether the two
  // bands reach a steady split, which is what Gamma->L transfer in GaAs will
  // be judged on. Relaxed atomics: these are diagnostics, not physics.
  static std::atomic<long> nXband{0};
  // occupancy in femtoseconds: atomic<double> would need C++20 fetch_add, and
  // integer fs is exact and plenty - a flight is O(1e-13 s).
  static std::atomic<long> bandFs[8] = {};
  std::vector<std::unique_ptr<FBS>> scatBand;
  for (int b = 0; b < engineNBands; b++)
    scatBand.emplace_back(new FBS(pkg, bs, b, T, cbm + 1.0, binW, mode));
  FBS &scat = *scatBand[0];
  for (auto &sb : scatBand)
    sb->setMaxDestBand(static_cast<std::size_t>(engineNBands - 1));
  if (engineNBands > 1) {
    bool perBand = true;
    for (auto &sb : scatBand)
      perBand = perBand && sb->hasPerSourceBandG2();
    if (perBand)
      std::printf("# ENGINE_NBANDS=%d: per-source-band g2 rows (pkg v0.4) "
                  "- interband transport is live\n", engineNBands);
    else
      std::printf("# ENGINE_NBANDS=%d: MACHINERY ONLY - this package has no "
                  "per-source-band g2 rows (v0.3), so every band is running "
                  "band-0 final states. Do NOT read the mobility as "
                  "interband physics.\n", engineNBands);
  }
  // BALANCE=1: repair detailed balance against the sampler's own bin DOS
  if (std::getenv("BALANCE") && std::atoi(std::getenv("BALANCE")) != 0) {
    const std::size_t np = scat.enforceDetailedBalance(KB * T);
    std::printf("# detailed balance enforced against the bin DOS: %zu pair(s), "
                "%zu grid points rewritten\n", np, scat.getRewrittenPoints());
    std::printf("# ab-initio rates moved by %.1f%% on average, %.1f%% max "
                "<-- this is the size of the correction, not physics\n",
                100.0 * scat.getRateChangeMean(),
                100.0 * scat.getRateChangeMax());
    if (scat.getRewrittenPoints() == 0)
      std::printf("# WARNING: balance repair was a NO-OP - check energy units\n");
  }
  // GAMMA0_SLAB=0 reverts to one global bound (A/B). Slabbing costs nothing
  // in accuracy - self-scattering is variance reduction, not physics - and
  // removes most of the wasted flights.
  //
  // GAMMA0_W is RESOLUTION ONLY. It used to be a correctness knob: the old
  // point-sampled bound leaked 274 breaches at 10 meV on this package and the
  // run printed "results are NOT valid", so the width had to be tuned upward
  // per package - fatal for a material-agnostic pipeline. buildGamma0Slabs now
  // bounds the rate per tetrahedron per slab, plus an explicit margin for the
  // energy a free flight traverses (Gamma0 is chosen at the flight's START
  // energy). Verified g0viol=0 in BOTH modes at 10 and 20 meV. 0.010 is the
  // cheaper of the two: mode 1 48.0% vs 59.5%, mode 0 71.1% vs 75.4%
  // self-scattering. That cost - against 43.8% for the old UNSAFE bound - is
  // the honest price of a bound that holds instead of one that happened to.
  const bool slab = !(std::getenv("GAMMA0_SLAB") &&
                      std::atoi(std::getenv("GAMMA0_SLAB")) == 0);
  // EVERY band instance needs its own slab table. This used to run on
  // `scat` alone (= band 0); a band-1 instance then had slabG0 empty,
  // getSlabWidth() returned 0, its flights were never capped, and it fell
  // back to the GLOBAL Gamma0. That is only invisible while nothing ever
  // reaches band 1 - the moment interband transport went live it showed up as
  // g_real/g_tab = 0.81.
  if (slab) {
    const double gw = std::getenv("GAMMA0_W")
                          ? std::atof(std::getenv("GAMMA0_W")) : 0.010;
    for (auto &sb : scatBand)
      sb->buildGamma0Slabs(gw);
  }
  const double gamma0 = scat.getGamma0();
  std::printf("# %s: CBM %.4f eV, Gamma0 = %.3e 1/s, %zu mechanism(s), %s\n",
              pkg.c_str(), cbm, gamma0, scat.getNrMechanisms(),
              scat.hasG2Tables() ? "k-RESOLVED + g2 FINAL STATES (Stage 2)"
              : (scat.isKResolved() ? "k-RESOLVED (phaseB)" : "energy tables (phaseA)"));
  std::printf("# bin-implied equilibrium <E>-CBM = %.4f eV  (physical 1.5kT = "
              "%.4f eV)\n", scat.equilibriumEnergy(KB * T), 1.5 * KB * T);
  std::printf("# ensemble: %d particles x %d replica(s) x %.0f ps\n",
              nPart, nRep, tTotal * 1e12);
  std::printf("# %8s %14s %14s %10s %14s %12s %10s\n", "F[V/cm]", "vd[m/s]",
              "mu[cm2/Vs]", "SE(mu)", "<E>-CBM[eV]", "selfscat[%]",
              "<cos_th>", "fail[%]", "accel", "<v^2>");

  // Build the thermal CDF ONCE, on this thread, before any worker starts.
  // Without this the first parallel replicas race on the lazily-built cache.
  for (auto &sb : scatBand)
    sb->prepareThermal(KB * T);

  std::vector<double> fields{500.0, 1000.0, 2000.0, 5000.0};
  if (argc > 2) {                      // comma-separated F list [V/cm]
    fields.clear();
    std::string s(argv[2]);
    size_t p = 0;
    while (p < s.size()) {
      size_t c = s.find(',', p);
      if (c == std::string::npos) c = s.size();
      fields.push_back(std::stod(s.substr(p, c - p)));
      p = c + 1;
    }
  }
  // field direction (argv[4]): "100" (default), "110", "111"
  double dir[3] = {1.0, 0.0, 0.0};
  if (argc > 4) {
    const std::string d(argv[4]);
    for (int i = 0; i < 3 && i < (int)d.size(); i++)
      dir[i] = d[i] - '0';
    const double n = std::sqrt(dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2]);
    for (int i = 0; i < 3; i++) dir[i] /= n;
    std::printf("# field along <%s>\n", d.c_str());
  }

  // one independent replica -> (mu, <E>, selfscat%)
  static constexpr int EHIST_N = 150;
  static constexpr double EHIST_DE = 0.002;   // 2 meV bins to 0.30 eV
  struct Res { double mu, energy, self, cosTheta, retention, speedRatio, failFrac, accel, v2, muEns;
               double gTab, gReal, dEmean, dEmax; std::vector<double> hist; };
  // DTCAP_SCALE: test knob. The flight cap must be physically INERT - it
  // only chops the trajectory more finely, it does not change the trajectory
  // or the scattering physics. Shrinking it by 10x must leave mu unchanged
  // within error bars; if it does not, the cap is doing something it should
  // not be. Default 0.5 (half a slab per flight).
  const double dtCapScale =
      std::getenv("DTCAP_SCALE") ? std::atof(std::getenv("DTCAP_SCALE")) : 0.5;
  // PER BAND. The flight cap keeps a flight's energy excursion inside half a
  // Gamma0 slab; with the slab-wise bound that is what makes the accept /
  // self-scatter bookkeeping consistent. Taking vMax and the slab width from
  // BAND 0 and applying them to a band-1 carrier under-restricts its
  // excursion, it crosses more slabs than the cap intends, and the identity
  // g_real = <Gamma> breaks - measured as ratio 0.810 instead of ~1.00 on the
  // first package where carriers actually reached band 1.
  std::vector<double> vMaxB, slabWB;
  for (int b = 0; b < engineNBands; b++) {
    vMaxB.push_back(bs.getMaxSpeed(b));
    slabWB.push_back(scatBand[b]->getSlabWidth());
  }
  const double vMax = vMaxB[0];
  const double slabW_ = slabWB[0];
  for (int b = 0; b < engineNBands; b++)
    std::printf("# band %d: vMax = %.3e m/s, slab = %.4f eV\n", b, vMaxB[b],
                slabWB[b]);
  auto runReplica = [&](int rep, double F) {
    std::mt19937_64 rng(1234ull + 7919ull * static_cast<unsigned>(rep));
    std::uniform_real_distribution<double> U01(0.0, 1.0);

    // thermal initialization: E ~ Boltzmann above CBM, k via iso-surface
    // sampler (also exercises the final-state machinery)
    // WP3-INTERBAND, engine half: the particle now CARRIES its band index
    // instead of every lookup hardcoding 0. Nothing can change it yet - the
    // package has no interband channel (`g2bins` rows carry no band dimension
    // and phaseA rates are band-0 only), so `band` stays 0 for the whole
    // trajectory and every number is bit-identical. This is the plumbing that
    // has to exist before an interband final state can be sampled at all.
    // Sized 2026-08-27: the fraction of the distribution above silicon's
    // band-1 minimum (136.7 meV) is 2.4% at 200 V/cm but 20.4% at 20 kV/cm.
    struct P {
      int band = 0;
      NBS::Vec3 k;
      std::int64_t hint = -1;
    };
    std::vector<P> ps(nPart);
    for (auto &p : ps) {
      double E;
      do
        E = scat.sampleThermalEnergy(KB * T, rng);   // DOS-weighted, not bare
      while (!scat.sampleFinalState(E, rng, p.k, p.hint));
    }

    double sumVdt = 0, sumEdt = 0, sumT = 0, sumV2dt = 0, sumV2G = 0;
    // Realized vs tabulated scattering rate. The self-scattering scheme makes
    // real events occur at rate Gamma(k); if the rate the ensemble ACTUALLY
    // experiences differs from the time-average of the table, the dynamics is
    // inconsistent with the rates the SERTA integral is reading.
    double sumGdt = 0;
    long nRealAcc = 0;
    // Time-weighted histogram of E-CBM. <E> alone cannot distinguish "thermal
    // on the wrong DOS" from "non-thermal": both can give the same mean. The
    // shape does.
    std::vector<double> hist(EHIST_N, 0.0);
    // PER-EVENT ENERGY CONSERVATION. The final state is sampled on the
    // iso-surface of the SAME band field the engine walks on, so
    // E(k_new) must equal the requested E' exactly. Any residual here is
    // the energy mismatch at its source, not inferred from <E>.
    double dEsum = 0, dEmax = 0;
    long dEn = 0;
    long nSelf = 0, nReal = 0;
    double cosAcc = 0;      // <cos(theta)> momentum-randomization diagnostic
    long cosN = 0;
    // MOMENTUM RETENTION, the quantity mobility actually depends on:
    //     R = <v_post . v_pre> / <|v_pre|^2>
    // <cos(theta)> normalises by |v_post||v_pre| and therefore DISCARDS any
    // correlation between the final SPEED and the scattering angle. If fast
    // final states sit preferentially backward and slow ones forward, <cos>
    // can be ~0 while R is negative - momentum is over-relaxed and the drift
    // mobility falls below the SERTA estimate even though the angles look
    // isotropic. That is exactly the signature seen in mode 0
    // (drift/mu_ens = 0.875 with <cos> = 0.002).
    double projAcc = 0, vpre2Acc = 0, spdAcc = 0;
    long nFail = 0;         // real events whose final state could not be sampled
    // effective acceleration between real collisions: must equal eF/m*.
    // The ballistic test only verifies this over long flights (many cells);
    // here the field moves k by ~2% of a cell between collisions, which is
    // the regime that actually sets the mobility.
    double dvAcc = 0, dtAcc = 0;
    for (auto &p : ps) {
      double t = 0;
      double tLast = 0, vxLast = 0;
      bool haveLast = false;
      // trapezoidal drift accumulation: sampling only the pre-flight velocity
      // loses the in-flight acceleration, a bias of order Gamma/Gamma0 (it
      // destroyed a constant-rate test package where Gamma0 = 1.2 Gamma, and
      // costs ~1% where self-scattering dominates). The end-of-flight velocity
      // is the next flight's start velocity, so the trapezoid is free.
      auto vCur = bs.getVelocity(p.k, p.band, p.hint);
      double eCur = bs.getEnergy(p.k, p.band, p.hint);
      // 0.5 slab / (F * vMax), PER BAND - indexed by the particle's current
      // band inside the flight loop, since a real event can change it.
      std::vector<double> dtCapB(engineNBands,
                                 std::numeric_limits<double>::infinity());
      for (int b = 0; b < engineNBands; b++)
        if (slab && slabWB[b] > 0 && F * vMaxB[b] > 0)
          dtCapB[b] = dtCapScale * slabWB[b] / (F * vMaxB[b]);
      while (t < tTotal) {
        // dispatch scattering on the particle's CURRENT band
        FBS &sc = *scatBand[p.band < (int)scatBand.size() ? p.band : 0];
        const double g0 = slab ? sc.gamma0At(eCur) : gamma0;
        const double tau = -std::log(1 - U01(rng)) / g0;
        // SLAB-CROSSING CAP. gamma0At(eCur) bounds the rate only within the
        // slab holding eCur, but the field moves the particle in energy during
        // the flight at dE/dt = F*v_parallel [eV/s]. Cap dt so it cannot cross
        // more than half a slab, then redraw with the new slab's bound - exact
        // for a piecewise-constant rate, because at the moment the rate changes
        // the residual waiting time is again exponential with the NEW rate.
        // Without this the bound is FIELD-BLIND: it held at 200 and 1000 V/cm
        // but broke at 5000 (g0viol=2695, mu 836 and meaningless), and v-E
        // curves are a WP3a deliverable that runs to exactly those fields.
        // Cap with the band's MAXIMUM speed, not this flight's starting
        // velocity. v_parallel can be ~0 at the start - giving an infinite cap
        // - immediately before the field accelerates the particle across
        // several slabs; that mistake left g0viol=141 at 5 kV/cm and 20660 at
        // 20 kV/cm. vMax makes the cap an upper bound on the excursion for any
        // flight, and it costs nothing at low field where dtCap greatly
        // exceeds the mean free time.
        // SERTA integrand must pair v and Gamma at the SAME phase-space
        // point. The rate used to be evaluated from p.k AFTER the flight had
        // already advanced it, while v2c is the velocity at flight START -
        // so mu_ens mixed the two ends of a flight. The size of that error
        // scales with how far a flight moves in energy, which differs per
        // mode (self-scattering 48% in mode 1 vs 71% in mode 0), which is why
        // mu/mu_ens came out 1.067 / 1.005 / 0.905 across modes that should
        // all sit at 1 for isotropic scattering.
        std::int64_t hStart = p.hint;
        const double gStart = sc.getTotalRate(p.k, eCur, hStart);
        const double dtCap =
            dtCapB[p.band < (int)dtCapB.size() ? p.band : 0];
        const double dt = std::min(std::min(tau, dtCap), tTotal - t);
        const bool capped = (dtCap < tau) && (dtCap < tTotal - t);
        // free flight: dk = -e F dt / hbar (electron)
        for (int i = 0; i < 3; i++)
          p.k[i] -= QE * F * dir[i] * dt / HBAR;
        const auto vEnd = bs.getVelocity(p.k, p.band, p.hint);
        if (t + dt > tTransient) {
          const double w = std::min(dt, t + dt - tTransient);
          const double vMid = 0.5 * ((vCur[0] + vEnd[0]) * dir[0] +
                                     (vCur[1] + vEnd[1]) * dir[1] +
                                     (vCur[2] + vEnd[2]) * dir[2]);
          sumVdt += vMid * w;
          sumEdt += (bs.getEnergy(p.k, p.band, p.hint) - cbm) * w;
          // <v^2>: if the final-state sampler's measure were biased, the
          // ensemble would sit on states with the wrong velocity statistics
          // and mu would deviate from the Boltzmann/SERTA integral even with
          // identical rates and full momentum randomisation
          const double v2c = vCur[0]*vCur[0] + vCur[1]*vCur[1] + vCur[2]*vCur[2];
          sumV2dt += v2c * w;
          // SERTA integrand on the MC's own ensemble: mu = e<v^2/Gamma>/(3kT).
          // Same states, same rates, no drift estimator involved.
          {
            const double gg = gStart;   // paired with vCur, same point
            if (gg > 0) sumV2G += v2c / gg * w;
            sumGdt += gg * w;
          }
          {
            const int hb = static_cast<int>((bs.getEnergy(p.k, p.band, p.hint) - cbm)
                                            / EHIST_DE);
            if (hb >= 0 && hb < EHIST_N) hist[hb] += w;
          }
          sumT += w;
          if (p.band < 8)
            bandFs[p.band].fetch_add(static_cast<long>(w * 1e15),
                                     std::memory_order_relaxed);
        }
        vCur = vEnd;
        t += dt;
        if (t >= tTotal)
          break; // reached tTotal mid-flight
        if (capped) {
          // flight truncated at a slab boundary, NOT a scattering event:
          // refresh the energy and redraw against the new slab's bound
          eCur = bs.getEnergy(p.k, p.band, p.hint);
          continue;
        }
        const double E = bs.getEnergy(p.k, p.band, p.hint);
        const double rateNow = sc.getTotalRate(p.k, E, p.hint);
        // the bound must never be exceeded; count it loudly if it is
        if (rateNow > g0)
          sc.noteG0Violation();
        eCur = E;
        if (U01(rng) * g0 < rateNow) {
          nReal++;
          if (t > tTransient) nRealAcc++;
          const std::size_t m = sc.selectMechanism(p.k, E, p.hint, rng);
          const double Ef = E + sc.getDeltaE(m);
          NBS::Vec3 kNew;
          std::int64_t hNew = p.hint;
          // Destination band, drawn from the row entry's `dest_band` (v0.4).
          // The sampler resolves the iso-surface on THAT band, so the energy
          // check below is a real test of the transfer: a final state placed
          // on the source band's surface would miss E' by the inter-band
          // separation and show up immediately in dE_max.
          int bNew = p.band;   // filled by the sampler when dest_band exists
          const auto vPre = bs.getVelocity(p.k, p.band, p.hint);
          if (sc.sampleFinalStateG2(m, p.k, Ef, rng, kNew, hNew, &bNew)) {
            // momentum-randomization diagnostic: <cos(theta)> between the
            // velocity before and after a real scattering event. Isotropic
            // final states give 0; a positive value means the event only
            // partially destroys momentum, so the transport (momentum)
            // relaxation time is tau/(1-<cos>) - the quantity that sets
            // mobility, not the total scattering time.
            {
              std::int64_t hc = hNew;
              const double eGot = bs.getEnergy(kNew, bNew, hc);
              const double d = std::fabs(eGot - Ef);
              dEsum += d;
              if (d > dEmax) dEmax = d;
              dEn++;
            }
            const auto vPost = bs.getVelocity(kNew, bNew, hNew);
            const double n1 = std::sqrt(vPre[0]*vPre[0] + vPre[1]*vPre[1] +
                                        vPre[2]*vPre[2]);
            const double n2 = std::sqrt(vPost[0]*vPost[0] + vPost[1]*vPost[1] +
                                        vPost[2]*vPost[2]);
            if (n1 > 0 && n2 > 0) {
              const double dotpp = vPre[0]*vPost[0] + vPre[1]*vPost[1] +
                                   vPre[2]*vPost[2];
              cosAcc += dotpp / (n1 * n2);
              projAcc += dotpp;          // UNnormalised: keeps speed weighting
              vpre2Acc += n1 * n1;
              spdAcc += n2 / n1;
              cosN++;
            }
            if (haveLast && t > tLast) {
              const double vxPre = vPre[0]*dir[0] + vPre[1]*dir[1] +
                                   vPre[2]*dir[2];
              dvAcc += vxPre - vxLast;
              dtAcc += t - tLast;
            }
            tLast = t;
            vxLast = vPost[0]*dir[0] + vPost[1]*dir[1] + vPost[2]*dir[2];
            haveLast = true;
            p.k = kNew;
            p.hint = hNew;
            // No clamp here. Destinations the engine cannot represent are
            // REJECTED inside the sampler (setMaxDestBand), so bNew is always
            // a band we have an instance for. Clamping after the draw put the
            // carrier at a k' taken from the destination band's iso-surface
            // while labelling it the source band - off the energy shell, and
            // invisible to the dE check above, which runs on the pre-clamp
            // band. It cost ~30% of the mobility at 200 V/cm.
            if (bNew != p.band) nXband.fetch_add(1, std::memory_order_relaxed);
            p.band = bNew;
            eCur = Ef;         // sampler lands exactly on the requested E'
            vCur = vPost;      // discontinuous change at a real event
          } else {
            nFail++;   // no final state: momentum NOT randomized this event
          }
        } else {
          nSelf++;
        }
      }
    }
    const double vd = sumVdt / sumT;
    // Sign convention by CARRIER. Electrons drift AGAINST F, so -vd/F is
    // positive for them. A HOLE package stores a flipped axis (E = VBM - E_e)
    // and its carriers drift WITH F - measured vd = +3.19e3 m/s at 200 V/cm
    // once the stored velocities were corrected - so the same formula would
    // report a negative mobility. `carrierSign` is +1 for electrons, -1 for
    // holes, read from the package's /bands/electron carrier attribute.
    return Res{-carrierSign * vd / F * 1e4,   // cm^2/Vs
               sumEdt / sumT,
               100.0 * nSelf / std::max(1L, nSelf + nReal),
               cosN ? cosAcc / cosN : 0.0,
               vpre2Acc > 0 ? projAcc / vpre2Acc : 0.0,
               cosN ? spdAcc / cosN : 0.0,
               100.0 * nFail / std::max(1L, nReal),
               dtAcc > 0 ? dvAcc / dtAcc : 0.0,
               sumV2dt / sumT,
               sumV2G / sumT / (3.0 * KB * T) * 1e4,
               sumGdt / sumT,
               nRealAcc / (nPart * (tTotal - tTransient)),
               dEn ? dEsum / dEn : 0.0, dEmax, hist};
  };

  const unsigned nThreads = std::min<unsigned>(
      static_cast<unsigned>(nRep),
      std::max(1u, std::thread::hardware_concurrency()));

  for (const double Fcm : fields) {
    const double F = Fcm * 100.0; // V/m, along dir
    std::vector<Res> res(nRep);
    std::atomic<int> next{0};
    auto worker = [&]() {
      int r;
      while ((r = next.fetch_add(1)) < nRep)
        res[r] = runReplica(r, F);
    };
    std::vector<std::thread> pool;
    for (unsigned i = 0; i < nThreads; i++)
      pool.emplace_back(worker);
    for (auto &th : pool)
      th.join();

    double muM = 0, eM = 0, sM = 0, cM = 0, fM = 0, aM = 0, v2M = 0, meM = 0;

    double retM = 0, spdM = 0;
    double gtM = 0, grM = 0, deM = 0, dxM = 0;
    std::vector<double> hAcc(EHIST_N, 0.0);
    for (const auto &r : res) {
      muM += r.mu; eM += r.energy; sM += r.self; cM += r.cosTheta;
      retM += r.retention; spdM += r.speedRatio;
      fM += r.failFrac; aM += r.accel; v2M += r.v2; meM += r.muEns;
      gtM += r.gTab; grM += r.gReal; deM += r.dEmean;
      for (int i = 0; i < EHIST_N; i++) hAcc[i] += r.hist[i];
      if (r.dEmax > dxM) dxM = r.dEmax;
    }
    muM /= nRep; eM /= nRep; sM /= nRep; cM /= nRep; fM /= nRep; aM /= nRep;
    retM /= nRep; spdM /= nRep;
    v2M /= nRep; meM /= nRep; gtM /= nRep; grM /= nRep; deM /= nRep;
    double se = 0;
    if (nRep > 1) {
      double var = 0;
      for (const auto &r : res) var += (r.mu - muM) * (r.mu - muM);
      se = std::sqrt(var / (nRep - 1) / nRep);   // standard error of the mean
    }
    std::printf("  %8.0f %14.4e %14.1f %10.1f %14.4f %12.1f %10.3f %9.1f"
                " %11.3e %11.3e   cfb=%zu  mu_ens=%.1f"
                "  g_tab=%.4e g_real=%.4e ratio=%.3f"
                "  dE_mean=%.3e dE_max=%.3e eV  g0viol=%zu"
                "  retention=%+.4f  spd_ratio=%.4f\n",
                Fcm, -muM * F * 1e-4, muM, se, eM, sM, cM, fM, aM, v2M,
                scat.getCentroidFallbacks(), meM, gtM, grM,
                gtM > 0 ? grM / gtM : 0.0, deM, dxM,
                scat.getG0Violations(), retM, spdM);
    if (std::getenv("EHIST")) {
      double tot = 0;
      for (double v : hAcc) tot += v;
      std::FILE *fh = std::fopen(std::getenv("EHIST"), "w");
      for (int i = 0; i < EHIST_N; i++)
        std::fprintf(fh, "%.4f %.6e\n", (i + 0.5) * EHIST_DE, hAcc[i] / tot);
      std::fclose(fh);
      std::printf("# energy histogram -> %s\n", std::getenv("EHIST"));
    }
  }
  std::printf("# a-weighted final states: %zu draws, %zu rejected (%.1f%%)\n",
              scat.getAWeightTries(), scat.getAWeightRej(),
              scat.getAWeightTries() ? 100.0*scat.getAWeightRej()/scat.getAWeightTries() : 0.0);
  if (scat.getG0Violations())
    std::printf("# WARNING: self-scattering bound breached %zu times - "
                "results are NOT valid; widen GAMMA0_W\n",
                scat.getG0Violations());
  if (scat.hasG2Tables()) {
    if (scat.getG2Rejected() && engineNBands == 1)
      std::printf("# NOTE: %zu interband final-state draws were REJECTED and "
                  "redrawn (ENGINE_NBANDS=1 on a v0.4 package): |g|^2 "
                  "restricted to intraband destinations, renormalised. This "
                  "is the no-interband CONTROL - rerun with ENGINE_NBANDS=2 "
                  "for the interband case.\n", scat.getG2Rejected());
    if (engineNBands > 1) {
      double tot = 0;
      for (int b = 0; b < engineNBands && b < 8; b++) tot += bandFs[b].load();
      std::printf("# interband: %ld band-changing events; time-weighted "
                  "occupancy =", nXband.load());
      for (int b = 0; b < engineNBands && b < 8; b++)
        std::printf(" b%d=%.3f%%", b, tot > 0 ? 100.0 * bandFs[b].load() / tot : 0.0);
      std::printf("\n");
    }
    const double h = scat.getG2Hits(), fb = scat.getG2Fallbacks();
    std::printf("# g2 row walks: %zu of %.0f accepted draws (%.2f%%) took a "
                "NON-random substitute; mean walk %.2f tets\n",
                scat.getG2Walks(), h,
                h > 0 ? 100.0 * scat.getG2Walks() / h : 0.0,
                scat.getG2Walks() ? 1.0 * scat.getG2WalkLen() / scat.getG2Walks() : 0.0);
    // Events that miss every g2 row entry fall back to the ISOTROPIC bin-DOS
    // sampler - a different measure from the |g|^2 weights the rates were
    // built against, so each one breaks detailed balance on the g2 path.
    std::printf("# g2 sampling: %.0f hits, %.0f fallbacks (%.2f%% of g2 events "
                "sampled from the WRONG measure)\n", h, fb,
                100.0 * fb / std::max(1.0, h + fb));
  }
  return 0;
}
