#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "../bulkSimulation/basicBulkParticleHandler.hpp"
#include "Ga2O3Functions.hpp"

#include <ParticleType/emcElectron.hpp>
#include <emcPauliExclusion.hpp>
#include <emcDevice.hpp>

/**
 * Hot-phonon-limited high-field electron transport in bulk β-Ga₂O₃.
 *
 * Every published β-Ga₂O₃ Monte Carlo study — including the recent full-band
 * first-principles ones — evaluates the polar-optical scattering rate with the
 * phonon population held at its equilibrium Bose-Einstein value. β-Ga₂O₃ is
 * however the wide-bandgap material in which that assumption is least safe: it
 * has the strongest Fröhlich coupling of the family, the lowest optical phonon
 * energies (so the emission rate is high), and a thermal conductivity 5-25x
 * below GaN's (so the acoustic reservoir that absorbs the decayed LO energy is
 * itself slow to unload).
 *
 * This driver runs the same transport problem two ways:
 *   --use_hpb 0 : equilibrium phonon occupation (the published assumption)
 *   --use_hpb 1 : q-resolved non-equilibrium occupation N_q evolved
 *                 self-consistently with the carrier ensemble, optionally
 *                 coupled to an acoustic reservoir (three-temperature model)
 * and reports the velocity-field curve, the mean energy, and the phonon
 * populations for both.
 *
 * Usage:
 *   ./hotPhononGa2O3 --fields 10,50,100,200,300,400 --use_hpb 1 --tau_lo 5e-12
 *
 * Output (per run, into --outdir):
 *   ga2o3_vE_<tag>.txt        field [kV/cm], v_drift [cm/s], <E> [eV], N_LO, T_LO, T_ac
 *   ga2o3_trace_<tag>_F<..>.txt  time trace of the same quantities
 */

const SizeType Dim = 3;
using NumType = double;
using MaterialType = emcMaterial<NumType>;
using DeviceType = emcDevice<NumType, Dim>;
using ParticleHandler = basicBulkParticleHandler<NumType, DeviceType>;
using MapIdxTypeToPartType = ParticleHandler::MapIdxToParticleTypes;

// ---------------------------------------------------------------------------
// minimal CLI parsing: --key value
// ---------------------------------------------------------------------------
std::map<std::string, std::string> parseCLI(int argc, char **argv) {
  std::map<std::string, std::string> cli;
  for (int i = 1; i + 1 < argc; i += 2) {
    std::string key = argv[i];
    if (key.rfind("--", 0) == 0)
      cli[key.substr(2)] = argv[i + 1];
  }
  return cli;
}

NumType getArg(const std::map<std::string, std::string> &cli,
               const std::string &key, NumType fallback) {
  auto it = cli.find(key);
  return (it == cli.end()) ? fallback : std::stod(it->second);
}

std::string getArgStr(const std::map<std::string, std::string> &cli,
                      const std::string &key, const std::string &fallback) {
  auto it = cli.find(key);
  return (it == cli.end()) ? fallback : it->second;
}

std::vector<NumType> parseList(const std::string &s) {
  std::vector<NumType> out;
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, ','))
    if (!item.empty())
      out.push_back(std::stod(item));
  return out;
}

// ---------------------------------------------------------------------------
struct RunConfig {
  bool useHPB;
  bool useAcoustic;
  NumType tauLO;
  NumType tauAc;
  NumType temperature;
  NumType doping;   // [1 / m³]
  NumType boxSide;  // [m]
  NumType dt;       // [s]
  NumType totalTime;
  NumType steadyFrac; // fraction of the run averaged for the steady state
  NumType alpha;
  NumType sigmaAc;
  // Scatter tables are rebuilt every reinitEvery steps rather than every step.
  // N_q evolves on the scale of tau_LO (ps), which is 4-5 orders of magnitude
  // slower than dt, so a modest stride is physically harmless and removes the
  // table rebuild as the runtime bottleneck. Convergence must be checked
  // against reinitEvery = 1.
  SizeType reinitEvery;
  bool useMultimode; // 5-mode ab initio polar set vs 44 meV single-mode lump
  bool useScreening; // free-carrier (plasmon) screening of the Froehlich vertex
  NumType screenEps; // background permittivity entering the Debye wavevector
  bool useImpurity;  // ionized-impurity (Brooks-Herring) scattering
  // q-resolved phonon occupation in the rates: use the 1/q-weighted mean over
  // the allowed window [|k_i-k_f|, k_i+k_f] instead of the DOS-weighted mean
  // over all q. Polar transitions are distributed as dq/q, so this is what
  // carriers actually sample; the DOS mean understates it whenever the excess
  // piles up at small q (the generic polar case).
  bool useQResolved;
  bool useQResAngle;
  // Pauli exclusion / Burstein-Moss band filling. The occupancy grid holds
  // N_c = floor(2 V dk^3 / 8pi^3) particles per k-slot, so a usable N_c needs a
  // large simulation volume; the production box gives N_c = 0 and would block
  // nothing. Enabling this also forces the SEQUENTIAL particle loop.
  bool usePauli;
  NumType pauliDk;
  NumType pauliEmax; // [eV] top of the energy range the grid represents
  SizeType seed;
  std::string outdir;
  std::string tag;
};

struct RunResult {
  NumType vDrift;   // [m / s]
  NumType energy;   // [eV]
  NumType meanNq;   // LO occupation
  NumType tempLO;   // [K], from the Planck inversion of <N_q>
  NumType tempAc;   // [K]
  NumType N0;       // equilibrium LO occupation at the lattice temperature
};

//! Planck inversion: temperature at which occupation N would be the
//! equilibrium value for a mode of energy energyEV.
NumType planckTemp(NumType energyEV, NumType N) {
  if (N <= 0.)
    return 0.;
  return constants::q * energyEV /
         (constants::kB * std::log(1. + 1. / N));
}

// ---------------------------------------------------------------------------
RunResult runOneField(NumType fieldkVcm, const RunConfig &cfg) {
  const NumType field = fieldkVcm * 1e5; // kV/cm → V/m
  const std::array<NumType, 3> maxPos = {cfg.boxSide, cfg.boxSide, cfg.boxSide};
  const NumType h = cfg.boxSide / 2.;
  const std::array<NumType, 3> spacing = {h, h, h};
  const std::array<NumType, 3> fieldDirection = {-1, 0, 0};
  const NumType Vsim = cfg.boxSide * cfg.boxSide * cfg.boxSide;

  DeviceType device{Ga2O3::getGa2O3Material<NumType>(), maxPos, spacing,
                    cfg.temperature};
  device.addConstantDopingRegion({0, 0, 0}, maxPos, cfg.doping);

  MapIdxTypeToPartType particleTypes;
  // 2000 energy levels up to 5 eV: electrons reach ~2 eV at 500 kV/cm.
  particleTypes[0] =
      std::make_unique<emcElectron<NumType, DeviceType>>(2000, 5., false);

  Ga2O3::addGammaValley(particleTypes[0], cfg.alpha);

  // Non-polar channels are identical in both cases.
  Ga2O3::addAcousticScattering(0, particleTypes[0], device, {0}, cfg.sigmaAc);
  Ga2O3::addNonPolarOpticalScattering(0, particleTypes[0], device, {0});
  if (cfg.useImpurity)
    Ga2O3::addImpurityScattering(0, particleTypes[0], device, {0});

  // Polar channel. Two independent switches:
  //   useMultimode — 5-mode ab initio polar set vs. the 44 meV single-mode lump
  //   useScreening — free-carrier (plasmon) screening of the Froehlich vertex
  // Both the equilibrium reference and the hot-phonon case use the SAME polar
  // set and the SAME screening, so the comparison stays like-for-like.
  const std::vector<NumType> modeEnergy =
      cfg.useMultimode ? Ga2O3::modeEnergyMulti
                       : std::vector<NumType>{Ga2O3::phononEnergyPOP};
  const std::vector<NumType> modeEpsLo =
      cfg.useMultimode ? Ga2O3::modeEpsLoMulti()
                       : std::vector<NumType>{Ga2O3::epsLo};

  auto screening =
      Ga2O3::makeScreening(cfg.useScreening, cfg.screenEps);
  // Seed with the thermal value so the first steps are not unscreened.
  screening->update(cfg.doping, cfg.temperature);

  std::vector<std::shared_ptr<Ga2O3::PhononBath>> baths;
  if (cfg.useHPB) {
    baths = Ga2O3::makePhononBaths(modeEnergy, cfg.tauLO, cfg.temperature, Vsim,
                                   cfg.useAcoustic, cfg.tauAc);
    for (auto &bath : baths)
      bath->setScreeningQ2(screening->getQs2());
    Ga2O3::addScreenedHotPolarOpticalScattering(
        0, particleTypes[0], {0}, modeEnergy, modeEpsLo, baths, screening,
        cfg.useQResolved, cfg.useQResAngle);
  } else {
    Ga2O3::addScreenedEquilibriumPolarOpticalScattering(
        0, particleTypes[0], {0}, cfg.temperature, modeEnergy, modeEpsLo,
        screening);
  }

  // k-space occupancy grid for band filling. kMax covers the highest carrier
  // energy the rate table represents.
  std::unique_ptr<emcPauliExclusion<NumType>> pauli;
  if (cfg.usePauli) {
    // kMax must cover the populated part of the distribution, not the whole
    // rate-table range. The grid is rebuilt every timestep, so its cost scales
    // as (kMax/dk)^3 and an over-large kMax dominates the runtime.
    const NumType eMax = cfg.pauliEmax;
    const NumType gMax = eMax * (1. + cfg.alpha * eMax);
    const NumType kMax = std::sqrt(2. * Ga2O3::relEffMass * constants::me *
                                   gMax * constants::q) / constants::hbar;
    pauli = std::make_unique<emcPauliExclusion<NumType>>(cfg.pauliDk, kMax, Vsim);
  }

  ParticleHandler handler(device, particleTypes, fieldDirection);
  handler.setSeed(cfg.seed);
  handler.resetAppliedFieldStrength(field);
  handler.generateInitialParticles();

  const SizeType nrSteps = std::ceil(cfg.totalTime / cfg.dt);
  const SizeType firstSteadyStep =
      static_cast<SizeType>((1. - cfg.steadyFrac) * nrSteps);

  std::ostringstream traceName;
  traceName << cfg.outdir << "/ga2o3_trace_" << cfg.tag << "_F"
            << static_cast<int>(fieldkVcm) << ".txt";
  std::ofstream trace(traceName.str());
  trace << "# t[s]  v_drift[m/s]  energy[eV]  N_LO  T_LO[K]  T_ac[K]\n";

  NumType sumV = 0., sumE = 0., sumNq = 0., sumTac = 0.;
  double pauliRejected = 0., pauliAttempted = 0.;
  SizeType nSteady = 0;

  // Coupling weights used to form the reported mean phonon occupation.
  const std::vector<NumType> wMode =
      cfg.useMultimode ? Ga2O3::modeWeightMulti : std::vector<NumType>{1.};

  for (SizeType step = 1; step <= nrSteps; ++step) {
    if (cfg.usePauli) {
      handler.moveParticleTypeWithBandFilling(cfg.dt, 0, *pauli);
      // counters reset each step, so accumulate for the run-average fraction
      // of scatter events refused for want of a free final state
      pauliRejected += pauli->nRejected;
      pauliAttempted += pauli->nScattered;
    } else {
      handler.moveParticles(cfg.dt);
    }

    const NumType v = handler.getAvgDriftVelocity(0)[0];
    const NumType e = handler.getAvgEnergy(0)[0];

    // Screening follows the CARRIER temperature, which rises under field and
    // weakens screening: T_e = 2<E>/(3 k_B).
    bool tablesStale = false;
    if (cfg.useScreening) {
      const NumType Te = 2. * e * constants::q / (3. * constants::kB);
      screening->update(cfg.doping, Te);
      // Keep the baths' transition weight in step with q_s, so the q-resolved
      // rate and angular sampler both see the screened vertex. Serial section:
      // this mutates prefix sums the parallel particle loop reads.
      for (auto &bath : baths)
        bath->setScreeningQ2(screening->getQs2());
      tablesStale = true;
    }

    if (cfg.useHPB) {
      for (auto &bath : baths)
        bath->update(cfg.dt);
      tablesStale = true;
    }
    // Rate tables depend on N_q and on q_s; refresh as those evolve.
    if (tablesStale && step % cfg.reinitEvery == 0)
      particleTypes[0]->reinitScatterTables();

    NumType nq = 0.;
    if (cfg.useHPB) {
      for (SizeType m = 0; m < baths.size(); ++m)
        nq += wMode[m] * baths[m]->getMeanNq();
    } else {
      for (SizeType m = 0; m < modeEnergy.size(); ++m)
        nq += wMode[m] /
              (std::exp(constants::q * modeEnergy[m] /
                        (constants::kB * cfg.temperature)) - 1.);
    }
    const NumType tac =
        cfg.useHPB ? baths[0]->getAcousticTemp() : cfg.temperature;

    if (step % 200 == 0 || step == nrSteps)
      trace << step * cfg.dt << " " << v << " " << e << " " << nq << " "
            << planckTemp(modeEnergy[0], nq) << " " << tac << "\n";

    if (step > firstSteadyStep) {
      sumV += v;
      sumE += e;
      sumNq += nq;
      sumTac += tac;
      ++nSteady;
    }
  }
  trace.close();

  // Dump the resolved phonon spectrum N_q(q). This is what tests the
  // confinement counter-argument of Smith et al. (arXiv:2402.06349): they argue
  // that in a 2DEG, energy-momentum conservation confines emission to a small
  // volume of q-space so the occupancy never departs from equilibrium "across
  // the spectrum of wavevectors". The bulk spectrum here shows directly how
  // wide the departure is, and where it sits.
  //
  // It also audits an approximation in this model: the scatter rates use the
  // DOS-weighted MEAN N_q, not the local N_q at the transition wavevector. If
  // the excess is concentrated at small q -- where the DOS weight q^2 is
  // smallest but polar coupling (1/q^2) is strongest -- then the mean
  // UNDERSTATES what the carriers actually experience, and the reported
  // suppression is conservative.
  if (cfg.useHPB) {
    std::ostringstream specName;
    specName << cfg.outdir << "/ga2o3_nq_" << cfg.tag << "_F"
             << static_cast<int>(fieldkVcm) << ".txt";
    std::ofstream spec(specName.str());
    spec << "# q [1/m]  then one N_q column per polar mode\n";
    spec << "# mode energies [eV]:";
    for (auto e : modeEnergy)
      spec << " " << e;
    spec << "\n# equilibrium N_0 per mode:";
    for (auto &b : baths)
      spec << " " << b->getN0();
    spec << "\n";
    for (SizeType i = 0; i < baths[0]->nrBins; ++i) {
      spec << baths[0]->qCentre(i);
      for (auto &b : baths)
        spec << " " << b->Nq[i];
      spec << "\n";
    }
    spec.close();
  }

  RunResult res;
  res.vDrift = sumV / nSteady;
  res.energy = sumE / nSteady;
  res.meanNq = sumNq / nSteady;
  // Reported T_LO is that of the LOWEST mode (the one that dominates phonon
  // absorption and hence the room-temperature mobility), evaluated from the
  // coupling-weighted mean occupation.
  res.tempLO = planckTemp(modeEnergy[0], res.meanNq);
  res.tempAc = sumTac / nSteady;
  res.N0 = 0.;
  for (SizeType m = 0; m < modeEnergy.size(); ++m)
    res.N0 += wMode[m] / (std::exp(constants::q * modeEnergy[m] /
                                   (constants::kB * cfg.temperature)) - 1.);

  if (cfg.usePauli && pauliAttempted > 0)
    std::cout << "    Pauli-blocked " << 100. * pauliRejected / pauliAttempted
              << "% of scatter events\n";

  handler.deleteParticles();
  return res;
}

// ---------------------------------------------------------------------------
int main(int argc, char **argv) {
  auto cli = parseCLI(argc, argv);

  RunConfig cfg;
  cfg.useHPB = getArg(cli, "use_hpb", 1) > 0.5;
  cfg.useAcoustic = getArg(cli, "use_acoustic", 1) > 0.5;
  cfg.tauLO = getArg(cli, "tau_lo", Ga2O3::tauLODefault);
  cfg.tauAc = getArg(cli, "tau_ac", Ga2O3::tauAcDefault);
  cfg.temperature = getArg(cli, "temp", 300.);
  cfg.doping = getArg(cli, "doping", 1e23);      // 1e17 cm^-3
  cfg.boxSide = getArg(cli, "box", 3e-7);        // 300 nm
  cfg.dt = getArg(cli, "dt", 1e-16);
  cfg.totalTime = getArg(cli, "time", 5e-12);
  cfg.steadyFrac = getArg(cli, "steady_frac", 0.4);
  cfg.alpha = getArg(cli, "alpha", Ga2O3::alpha);
  cfg.sigmaAc = getArg(cli, "sigma_ac", Ga2O3::sigmaAc);
  cfg.reinitEvery =
      std::max<SizeType>(1, static_cast<SizeType>(getArg(cli, "reinit_every", 1)));
  cfg.useMultimode = getArg(cli, "use_multimode", 0) > 0.5;
  cfg.useScreening = getArg(cli, "use_screening", 0) > 0.5;
  // Background permittivity for the Debye wavevector; epsLo is the default,
  // epsHi gives the stronger (shorter-range) screening.
  cfg.screenEps = getArg(cli, "screen_eps", Ga2O3::epsLo);
  cfg.useImpurity = getArg(cli, "use_impurity", 0) > 0.5;
  cfg.useQResolved = getArg(cli, "use_qresolved", 0) > 0.5;
  cfg.useQResAngle = getArg(cli, "qres_angle", 1) > 0.5;
  cfg.usePauli = getArg(cli, "use_pauli", 0) > 0.5;
  cfg.pauliDk = getArg(cli, "pauli_dk", 1.5e8);
  cfg.pauliEmax = getArg(cli, "pauli_emax", 2.0);
  cfg.seed = static_cast<SizeType>(getArg(cli, "seed", 1));
  cfg.outdir = getArgStr(cli, "outdir", ".");
  cfg.tag = getArgStr(cli, "tag", cfg.useHPB ? "hpb" : "eq");

  const SizeType nrThreads = static_cast<SizeType>(getArg(cli, "threads", 4));
#ifdef _OPENMP
  omp_set_num_threads(nrThreads);
#endif

  auto fields = parseList(getArgStr(cli, "fields", "10,50,100,150,200,300,400"));

  std::cout << "β-Ga₂O₃ bulk transport — "
            << (cfg.useHPB ? "NON-EQUILIBRIUM (hot) phonons"
                           : "EQUILIBRIUM phonons (published assumption)")
            << "\n";
  std::cout << "  T = " << cfg.temperature << " K,  n = " << cfg.doping
            << " m^-3,  box = " << cfg.boxSide * 1e9 << " nm\n";
  std::cout << "  m* = " << Ga2O3::relEffMass << " m0,  alpha = " << cfg.alpha
            << " 1/eV\n";
  if (cfg.useMultimode) {
    auto eps = Ga2O3::modeEpsLoMulti();
    std::cout << "  polar set: 5-mode ab initio (Santia 2019 Table VI, "
              << "sum-rule normalized)\n";
    for (SizeType m = 0; m < Ga2O3::modeEnergyMulti.size(); ++m)
      std::cout << "    mode " << m << ": hw = " << Ga2O3::modeEnergyMulti[m] * 1e3
                << " meV,  w = " << Ga2O3::modeWeightMulti[m]
                << ",  eps_lo,eff = " << eps[m] << "\n";
  } else {
    std::cout << "  polar set: SINGLE-MODE lump, hw_POP = "
              << Ga2O3::phononEnergyPOP * 1e3 << " meV (Ma 2016)\n";
  }
  std::cout << "  plasmon screening: " << (cfg.useScreening ? "ON" : "OFF")
            << (cfg.useScreening ? " (static Debye at T_e - conservative)" : "")
            << "\n";
  std::cout << "  phonon occupation in rates: "
            << (cfg.useQResolved
                    ? (cfg.useQResAngle
                           ? "q-RESOLVED rate and angle"
                           : "q-RESOLVED rate, UNIFORM angle (inconsistent)")
                    : "mean-field (DOS-weighted over all q)")
            << "\n";
  std::cout << "  band filling: "
            << (cfg.usePauli ? "Pauli exclusion ON (sequential loop)" : "off")
            << "\n";
  std::cout << "  ionized impurities: "
            << (cfg.useImpurity ? "ON (Brooks-Herring, N_I = doping)" : "OFF")
            << "\n";
  if (cfg.useHPB) {
    std::cout << "  tau_LO = " << cfg.tauLO * 1e12 << " ps  [ESTIMATED]"
              << ",  acoustic bath " << (cfg.useAcoustic ? "ON" : "OFF");
    if (cfg.useAcoustic)
      std::cout << ",  tau_ac = " << cfg.tauAc * 1e12 << " ps  [ESTIMATED]";
    std::cout << "\n  q-resolved bath: " << Ga2O3::nrPhononBins << " bins x "
              << Ga2O3::dqBin << " 1/m\n";
  }
  std::cout << "  dt = " << cfg.dt << " s,  t_total = " << cfg.totalTime
            << " s,  seed = " << cfg.seed << "\n\n";

  std::ostringstream sumName;
  sumName << cfg.outdir << "/ga2o3_vE_" << cfg.tag << ".txt";
  std::ofstream summary(sumName.str());
  summary << "# beta-Ga2O3 velocity-field, "
          << (cfg.useHPB ? "non-equilibrium" : "equilibrium") << " phonons\n";
  summary << "# tau_LO=" << cfg.tauLO << " s  tau_ac=" << cfg.tauAc
          << " s  T=" << cfg.temperature << " K  n=" << cfg.doping << " m^-3\n";
  summary << "# F[kV/cm]  v[cm/s]  <E>[eV]  N_LO  N_LO/N_0  T_LO[K]  T_ac[K]\n";

  auto start = std::chrono::high_resolution_clock::now();
  for (auto F : fields) {
    auto res = runOneField(F, cfg);
    const NumType vcms = std::abs(res.vDrift) * 100.; // m/s → cm/s

    summary << std::setw(8) << F << " " << std::scientific << std::setprecision(4)
            << vcms << " " << std::fixed << std::setprecision(4) << res.energy
            << " " << res.meanNq << " " << res.meanNq / res.N0 << " "
            << std::setprecision(1) << res.tempLO << " " << res.tempAc << "\n";
    summary.flush();

    std::cout << "F = " << std::setw(6) << F << " kV/cm   v = "
              << std::scientific << std::setprecision(3) << vcms << " cm/s   "
              << std::fixed << std::setprecision(3) << "<E> = " << res.energy
              << " eV   N_LO = " << std::setprecision(3) << res.meanNq
              << " (x" << std::setprecision(2) << res.meanNq / res.N0 << " N_0)"
              << "   T_LO = " << std::setprecision(0) << res.tempLO << " K"
              << "   T_ac = " << res.tempAc << " K\n";
  }
  summary.close();

  auto end = std::chrono::high_resolution_clock::now();
  std::cout << "\nCPU time: "
            << std::chrono::duration_cast<std::chrono::seconds>(end - start)
                   .count()
            << " s\nWrote " << sumName.str() << "\n";
  return 0;
}
