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
#include <Averages/emcBulkAverages.hpp>
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
  // T_TOTAL / T_TRANSIENT [ps] override the run length: a degenerate valence
  // manifold relaxes its BRANCH populations only through interbranch
  // channels, and 20 ps may not reach the DOS-thermal split (Si valence at
  // 300 K: 82.4 / 14.8 / 2.7% on the package mesh).
  const double tTotal = (std::getenv("T_TOTAL") ? std::atof(std::getenv("T_TOTAL")) : 20.0) * 1e-12;
  const double tTransient = (std::getenv("T_TRANSIENT") ? std::atof(std::getenv("T_TRANSIENT")) : 5.0) * 1e-12;

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
  // NO SIGN FLIP FOR HOLES, and this is the representation's doing rather than
  // an oversight. The flipped axis (E_h = VBM - E_e) with NEGATED velocities
  // turns a hole into a formal ELECTRON on the inverted band: same charge sign
  // in the equation of motion the engine integrates, same v = grad(E)/hbar on
  // the stored field. Its k-trajectory is the mirror of the physical hole's,
  // so the drift it measures is MINUS the physical drift, and mu = -vd/F
  // recovers the correct POSITIVE mobility for both carriers.
  //
  // A `carrierSign = -1` flip was added on 2026-08-31 on the strength of
  // "measured vd = +3.19e3 m/s, so holes drift WITH the field". That number
  // came from the results column below, which printed `-mu * F * 1e-4` - the
  // mobility back-converted, NOT the drift velocity, and already sign-flipped
  // for a hole package. The flip therefore double-counted a sign the
  // representation had handled, and turned mu_h = +534 into -534. The column
  // now prints the true drift velocity so the mistake is not reachable again.
  const double carrierSign = 1.0;
  if (packageIsHole(pkg))
    std::printf("# HOLE package: flipped axis (E = VBM - E_e), velocities "
                "negated at build time; the carrier is a formal electron on "
                "the inverted band, so mu = -vd/F unchanged\n");
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
  // DEFAULTS TO THE PACKAGE'S BAND COUNT, not to 1.
  //
  // It used to default to 1, which silently made EVERY run a no-interband
  // control: with maxDestBand = 0 the sampler rejects every draw whose
  // dest_band is not 0 and redraws, so |g|^2 is restricted to intraband
  // destinations. On a two-band electron package that is nearly harmless -
  // band 1 sits 136.7 meV above the edge and is under 1% occupied at 500 V/cm.
  // On silicon's three-fold DEGENERATE valence manifold, where every band
  // touches the edge, it discards most of the available final states: a run on
  // 2026-09-02 rejected 147k draws and was read as a hole-transport defect.
  // The engine printed a NOTE saying so; a default nobody has to read is
  // better than a warning everybody must.
  //
  // ENGINE_NBANDS=1 remains the deliberate no-interband CONTROL.
  // ...but only when the package can actually SUPPLY per-source-band rows. On
  // a v0.3 package (flat g2bins) more than one instance runs every band on
  // band-0 final states - "MACHINERY ONLY", which the engine says not to read
  // as physics - so the format decides: v0.4 -> all bands, v0.3 -> 1. Band 0
  // is constructed first to ask.
  std::vector<std::unique_ptr<FBS>> scatBand;
  scatBand.emplace_back(new FBS(pkg, bs, 0, T, cbm + 1.0, binW, mode));
  const int engineNBands =
      std::getenv("ENGINE_NBANDS")
          ? std::atoi(std::getenv("ENGINE_NBANDS"))
          : (scatBand[0]->hasPerSourceBandG2()
                 ? static_cast<int>(bs.getNrBands()) : 1);
  // Interband bookkeeping. Time-weighted occupancy is the number that matters
  // - an event count says transfers happen, occupancy says whether the two
  // bands reach a steady split, which is what Gamma->L transfer in GaAs will
  // be judged on. Relaxed atomics: these are diagnostics, not physics.
  static std::atomic<long> nXband{0};
  // occupancy in femtoseconds: atomic<double> would need C++20 fetch_add, and
  // integer fs is exact and plenty - a flight is O(1e-13 s).
  static std::atomic<long> bandFs[8] = {};
  for (int b = 1; b < engineNBands; b++)
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
  struct Res { double mu, vd, energy, self, cosTheta, retention, speedRatio, failFrac, accel, v2, muEns, muEnsT;
               // ORDER MATCHES THE POSITIONAL INITIALIZER in runReplica: dif/muE/muET
               // come after dEmax there. Putting them after muEnsT once shifted
               // every later field by three slots and printed D as g_tab.
               double gTab, gReal, dEmean, dEmax, dif, muE, muET; std::vector<double> hist; };
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
      // real-space displacement, for the MSD / Einstein-mobility estimator.
      // Bulk is homogeneous, so r is a diagnostic only and never feeds the
      // dynamics; it is integrated with the same trapezoid as the drift.
      double r[3] = {0, 0, 0}, r0[3] = {0, 0, 0};
      int is = 0;          // next MSD checkpoint to record
    };
    std::vector<P> ps(nPart);
    // INITIAL BAND drawn from each band's thermal weight. Every particle used
    // to start in band 0 (P::band defaults to 0 and nothing set it). For
    // electrons that is the CBM band and harmless; for a degenerate valence
    // manifold it is the heavy branch alone, and the ensemble then relaxes
    // through the interbranch channels for the whole run - Si holes at 300 K
    // sat at 44/34/22% after 20 ps against a DOS-thermal 82/15/3%, i.e. the
    // light branches over-populated 2-8x, which inflates mu_h. INIT_BAND=0
    // restores the old behaviour for A/B.
    std::vector<double> zb(scatBand.size(), 0.0);
    double ztot = 0;
    for (std::size_t b = 0; b < scatBand.size(); b++) { zb[b] = scatBand[b]->thermalWeight(KB * T); ztot += zb[b]; }
    const bool initBandThermal = !(std::getenv("INIT_BAND") && std::atoi(std::getenv("INIT_BAND")) == 0);
    for (auto &p : ps) {
      if (initBandThermal && scatBand.size() > 1 && ztot > 0) {
        double u = U01(rng) * ztot, acc = 0;
        p.band = 0;
        for (std::size_t b = 0; b < scatBand.size(); b++) { acc += zb[b]; if (u <= acc) { p.band = static_cast<int>(b); break; } }
      }
      FBS &si = *scatBand[p.band < (int)scatBand.size() ? p.band : 0];
      double E;
      do
        E = si.sampleThermalEnergy(KB * T, rng);   // DOS-weighted, not bare
      while (!si.sampleFinalState(E, rng, p.k, p.hint));
    }
    if (rep == 0 && scatBand.size() > 1) {
      std::printf("# initial band split (thermal weights):");
      for (std::size_t b = 0; b < scatBand.size(); b++) std::printf(" b%zu=%.1f%%", b, ztot > 0 ? 100 * zb[b] / ztot : 0.0);
      std::printf("%s\n", initBandThermal ? "" : "   [INIT_BAND=0: all particles start in band 0]");
    }

    // The time-weighted moments - drift velocity, <E>, <v^2>, the SERTA
    // integrand <v^2/Gamma>, <Gamma> and the energy histogram - all live in
    // emcBulkAverages now rather than as loose accumulators here. They used to
    // be six private sums plus a hand-rolled histogram, which is why every new
    // bulk observable meant editing this example.
    //
    // The histogram is part of it: <E> alone cannot distinguish "thermal on
    // the wrong DOS" from "non-thermal", since both can give the same mean.
    // The shape does. <Gamma> is carried for the same reason - the
    // self-scattering scheme makes real events occur at rate Gamma(k), so if
    // the rate the ensemble ACTUALLY experiences differs from the time average
    // of the table, the dynamics is inconsistent with the rates the SERTA
    // integral reads.
    emcBulkAverages<double> avg(1, EHIST_N, EHIST_DE);
    // MSD checkpoints every 1 ps from the end of the transient; checkpoint 0
    // is each particle's reference position. D is fitted from checkpoint
    // MSD_FIT_FROM on, past the ballistic t^2 regime (tau ~ 1e-13 s).
    static constexpr int NT_MSD = 15, MSD_FIT_FROM = 5;
    static constexpr double DT_MSD = 1e-12;
    emcDisplacementStatistics<double> disp(NT_MSD, tTransient, DT_MSD);
    long nRealAcc = 0;
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
        {
          // displacement: trapezoid over the substep, and every MSD
          // checkpoint that falls inside (t, t+dt] is recorded at its exact
          // time by linear interpolation within the step. No RNG use, so
          // the trajectory stream is untouched.
          const double vm[3] = {0.5 * (vCur[0] + vEnd[0]),
                                0.5 * (vCur[1] + vEnd[1]),
                                0.5 * (vCur[2] + vEnd[2])};
          while (p.is < NT_MSD) {
            const double ts = tTransient + p.is * DT_MSD;
            if (ts > t + dt) break;
            const double f = std::max(0.0, ts - t);
            double rs[3];
            for (int c = 0; c < 3; c++) rs[c] = p.r[c] + vm[c] * f;
            if (p.is == 0) {
              for (int c = 0; c < 3; c++) p.r0[c] = rs[c];
            } else {
              double dr[3];
              for (int c = 0; c < 3; c++) dr[c] = rs[c] - p.r0[c];
              disp.addSample(static_cast<SizeType>(p.is), dr);
            }
            p.is++;
          }
          for (int c = 0; c < 3; c++) p.r[c] += vm[c] * dt;
        }
        if (t + dt > tTransient) {
          const double w = std::min(dt, t + dt - tTransient);
          const double vMid = 0.5 * ((vCur[0] + vEnd[0]) * dir[0] +
                                     (vCur[1] + vEnd[1]) * dir[1] +
                                     (vCur[2] + vEnd[2]) * dir[2]);
          // <v^2>: if the final-state sampler's measure were biased, the
          // ensemble would sit on states with the wrong velocity statistics
          // and mu would deviate from the Boltzmann/SERTA integral even with
          // identical rates and full momentum randomisation. It is evaluated
          // at vCur, the point gStart was paired with, while the energy is
          // taken AFTER the k-update - the pairing this loop has always used.
          const double v2c = vCur[0]*vCur[0] + vCur[1]*vCur[1] + vCur[2]*vCur[2];
          avg.addSample(vMid, bs.getEnergy(p.k, p.band, p.hint) - cbm,
                        v2c, gStart, w);
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
    // Diffusion from the MSD slope, averaged over the components TRANSVERSE to
    // the field (all three at F = 0): the parallel component also carries the
    // field-induced velocity dispersion. Einstein mobility on the lattice
    // temperature (the physical statement) and on the ensemble's own
    // temperature 2<E>/3 (the self-consistent one - they coincide only when
    // the ensemble is thermal at the lattice, which the band-edge mesh error
    // prevents on coarse packages).
    double dTr = 0;
    {
      int nc = 0;
      for (int c = 0; c < 3; c++) {
        if (F == 0.0 || std::fabs(dir[c]) < 1e-9) {
          dTr += disp.getDiffusion(c, MSD_FIT_FROM);
          nc++;
        }
      }
      dTr = nc ? dTr / nc : 0.0;
    }
    const double muE = emcDisplacementStatistics<double>::getEinsteinMobility(dTr, KB * T);
    const double muET = emcDisplacementStatistics<double>::getEinsteinMobility(
        dTr, 2.0 * avg.getMeanEnergy() / 3.0);
    // Sign, unit scaling and the SERTA integrand are emcBulkAverages' job now,
    // so a second caller cannot re-derive one of them slightly differently.
    // Sign convention by CARRIER. Electrons drift AGAINST F, so -vd/F is
    // positive for them. A HOLE package stores a flipped axis (E = VBM - E_e)
    // and its carriers drift WITH F - measured vd = +3.19e3 m/s at 200 V/cm
    // once the stored velocities were corrected - so the same formula would
    // report a negative mobility. `carrierSign` is +1 for electrons, -1 for
    // holes, read from the package's /bands/electron carrier attribute.
    return Res{avg.getMobility(F, carrierSign),          // cm^2/Vs
               avg.getDriftVelocity(),                  // m/s, AS MEASURED
               avg.getMeanEnergy(),
               100.0 * nSelf / std::max(1L, nSelf + nReal),
               cosN ? cosAcc / cosN : 0.0,
               vpre2Acc > 0 ? projAcc / vpre2Acc : 0.0,
               cosN ? spdAcc / cosN : 0.0,
               100.0 * nFail / std::max(1L, nReal),
               dtAcc > 0 ? dvAcc / dtAcc : 0.0,
               avg.getMeanSquaredVelocity(),
               avg.getMobilitySERTA(KB * T),
               avg.getMobilitySERTASelfConsistent(),   // 3kT -> 2<E>
               avg.getMeanScatterRate(),
               nRealAcc / (nPart * (tTotal - tTransient)),
               dEn ? dEsum / dEn : 0.0, dEmax, dTr, muE, muET,
               avg.getEnergyHistogram()};
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

    double muM = 0, eM = 0, sM = 0, cM = 0, fM = 0, aM = 0, v2M = 0, meM = 0, metM = 0;
    double vdM = 0;   // measured, not back-converted from mu
    double dM = 0, muEM = 0, muETM = 0;
    emcReplicaStatistics<double> repStat;   // SE of the Einstein mobility

    double retM = 0, spdM = 0;
    double gtM = 0, grM = 0, deM = 0, dxM = 0;
    std::vector<double> hAcc(EHIST_N, 0.0);
    for (const auto &r : res) {
      muM += r.mu; vdM += r.vd; eM += r.energy; sM += r.self; cM += r.cosTheta;
      dM += r.dif; muEM += r.muE; muETM += r.muET; repStat.add("muE", r.muE);
      retM += r.retention; spdM += r.speedRatio;
      fM += r.failFrac; aM += r.accel; v2M += r.v2; meM += r.muEns; metM += r.muEnsT;
      gtM += r.gTab; grM += r.gReal; deM += r.dEmean;
      for (int i = 0; i < EHIST_N; i++) hAcc[i] += r.hist[i];
      if (r.dEmax > dxM) dxM = r.dEmax;
    }
    muM /= nRep; vdM /= nRep; dM /= nRep; muEM /= nRep; muETM /= nRep; eM /= nRep; sM /= nRep; cM /= nRep; fM /= nRep; aM /= nRep;
    retM /= nRep; spdM /= nRep;
    v2M /= nRep; meM /= nRep; metM /= nRep; gtM /= nRep; grM /= nRep; deM /= nRep;
    double se = 0;
    if (nRep > 1) {
      double var = 0;
      for (const auto &r : res) var += (r.mu - muM) * (r.mu - muM);
      se = std::sqrt(var / (nRep - 1) / nRep);   // standard error of the mean
    }
    std::printf("  %8.0f %14.4e %14.1f %10.1f %14.4f %12.1f %10.3f %9.1f"
                " %11.3e %11.3e   cfb=%zu  mu_ens=%.1f  mu_ensT=%.1f  D=%.3fcm2/s  mu_E=%.1f+/-%.1f  mu_ET=%.1f"
                "  g_tab=%.4e g_real=%.4e ratio=%.3f"
                "  dE_mean=%.3e dE_max=%.3e eV  g0viol=%zu"
                "  retention=%+.4f  spd_ratio=%.4f\n",
                // vd is MEASURED now. It used to be `-muM * F * 1e-4`, the
                // mobility back-converted, which for a hole package printed
                // the NEGATIVE of the real drift velocity and was read as
                // evidence that holes drift with the field.
                Fcm, vdM, muM, se, eM, sM, cM, fM, aM, v2M,
                scat.getCentroidFallbacks(), meM, metM, dM * 1e4, muEM,
                repStat.getStandardError("muE"), muETM, gtM, grM,
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
  {
    std::size_t na = 0, np_ = 0, nb_ = 0, nc = 0, nf = 0;
    for (auto &sb : scatBand) {
      na += sb->getCoulombAnalytic(); np_ += sb->getCoulombProposals();
      nb_ += sb->getCoulombBracketFail(); nc += sb->getCoulombCapped();
      nf += sb->getCoulombFallbacks();
    }
    if (na || nf)
      std::printf("# coulomb kernel-corrected placement: %zu events, %.2f proposals/event, "
                  "%zu bracket failures, %zu weight-capped (%.2f%% of proposals), "
                  "%zu fell back to rows\n", na, na ? double(np_) / na : 0.0, nb_,
                  nc, np_ ? 100.0 * nc / np_ : 0.0, nf);
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
