/**
 * Hot-carrier dynamics and device metrics in metal-halide perovskites.
 *
 * Bipolar ensemble Monte Carlo of photo-excited carrier cooling in a bulk
 * metal-halide perovskite, built on a three-temperature carrier/LO/acoustic
 * phonon model. Extends the single-carrier hot-phonon-bottleneck (HPB)
 * treatment of
 *
 *   Faber, Filipovic, Koster, J. Phys. Chem. Lett. 2024, 15, 12601
 *
 * with holes and electron-hole scattering, radiative and Auger recombination,
 * Pauli band filling, the Ridley (LO -> LA + TO) decay channel, multi-mode
 * Froehlich coupling, and energy-selective extraction from which the
 * open-circuit voltage and power conversion efficiency are computed.
 *
 * Every mechanism and material parameter is set from the command line. The
 * defaults are MAPbI3; pass --material sn for CsSnI3 presets. Running with no
 * arguments gives the full default model. Examples:
 *
 *   ./hotCarrierMHP                                  # full default model
 *   ./hotCarrierMHP --use_hot_phonon 0               # equilibrium reference
 *   ./hotCarrierMHP --material sn --E_ex 0.2 --delta_E 0.05
 *
 * Output files carry a suffix encoding the enabled mechanisms (for example
 * carrierTempHPB_AC_ESC.txt): avgEnergyElectrons*, avgEnergyHoles*,
 * phononOccupation*, nrCarriers*, carrierTemp*. The console reports V_OC,
 * yield, and PCE. A nonzero --seed makes a run reproducible; the seed and the
 * OMP thread count together define the trajectory.
 */

#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <string>

// Reuse the basicBulkParticleHandler from the bulk simulation example.
// Copy or symlink that file, or adjust the path as needed.
#include "../bulkSimulation/basicBulkParticleHandler.hpp"

#include <emcPauliExclusion.hpp>
#include <emcHotCarrierOutput.hpp>

#include <ParticleType/emcElectron.hpp>
#include <ParticleType/emcHole.hpp>
#include <ScatterMechanisms/emcFroehlichInteraction.hpp>
#include <ScatterMechanisms/emcHotPhononFroehlichMechanism.hpp>
#include <ValleyTypes/emcNonParabolicIsotropValley.hpp>
#include <emcCarrierCarrierScatter.hpp>
#include <emcInterCarrierScatter.hpp>
#include <emcDevice.hpp>
#include <emcPhononBath.hpp>
#include <emcEnergySelectiveContact.hpp>
#include <emcRecombination.hpp>

// Physics flags are now runtime CLI arguments (--use_hot_phonon 0/1 etc.)
// Default values (matching the original constexprs) are set in main().

const SizeType Dim = 3;
using NumType = double;
using MaterialType = emcMaterial<NumType>;
using DeviceType = emcDevice<NumType, Dim>;
using ParticleHandler = basicBulkParticleHandler<NumType, DeviceType>;
using MapIdxTypeToPartType = ParticleHandler::MapIdxToParticleTypes;

// ---- Fixed simulation numerics (not swept) -----------------------------
static constexpr NumType dt = 5e-15;        // [s]  5 fs — resolves fastest scatter
static constexpr SizeType nrStepsBetweenOutput = 10; // every 0.05 ps
static constexpr SizeType nrPhononBins = 1;
static constexpr NumType dq = 2e9;          // [m^{-1}]  single-mode phonon bin
static constexpr NumType alpha_np = 0.;     // non-parabolicity (parabolic approx)
static constexpr NumType dk_pauli  = 4e7;   // [m^{-1}]  Pauli k-grid spacing
static constexpr NumType kMax_pauli = 4e9;  // [m^{-1}]  covers all photoexcited k

// Box size and Vsim are now set in main() after USE_BAND_FILLING is read from CLI.

// -----------------------------------------------------------------------
// CLI argument helpers: --key value pairs
// -----------------------------------------------------------------------
static std::map<std::string, std::string> parseArgs(int argc, char *argv[]) {
  std::map<std::string, std::string> m;
  for (int i = 1; i + 1 < argc; i += 2) {
    std::string k = argv[i];
    if (k.size() > 2 && k[0] == '-' && k[1] == '-')
      m[k.substr(2)] = argv[i + 1];
  }
  return m;
}

template <class T>
static T getArg(const std::map<std::string, std::string> &m,
                const std::string &key, T def) {
  auto it = m.find(key);
  if (it == m.end()) return def;
  if constexpr (std::is_same_v<T, std::string>) return it->second;
  else if constexpr (std::is_same_v<T, bool>)   return std::stoi(it->second) != 0;
  else return static_cast<T>(std::stod(it->second));
}

// -----------------------------------------------------------------------

int main(int argc, char *argv[]) {
  auto cli = parseArgs(argc, argv);

  // ---- Physics model flags (CLI-overridable) ---------------------------
  // Defaults reproduce the original behaviour (HPB + acoustic + ESC, no BF).
  const bool USE_HOT_PHONON    = getArg(cli, "use_hot_phonon", true);
  const bool USE_ESC           = getArg(cli, "use_esc",        true);
  // Radiative + Auger recombination. Disable (--use_recomb 0) for clean
  // phonon-cascade studies: Auger returns E_gap per event to the carrier
  // gas, which at slow acoustic sinks drives a cumulative bath heat-up.
  const bool USE_RECOMB        = getArg(cli, "use_recomb",     true);
  // Carrier-carrier (intra + inter species) scattering. Disable (--use_cc 0)
  // for the mechanism-ablation study: c-c scattering thermalizes the
  // distribution shape without removing energy to the lattice.
  const bool USE_CC            = getArg(cli, "use_cc",         true);
  const bool USE_ACOUSTIC_BATH = getArg(cli, "use_acoustic",   true);
  const bool USE_BAND_FILLING  = getArg(cli, "use_bf",         false);
  // Box side [m]; default 100 nm (~1000 particles/species at n = 1e24 m^-3).
  // Larger boxes reduce the finite-size shot noise that seeds the N_q
  // feedback loop into metastable branches at small phonon quanta.
  const NumType boxSideCLI     = getArg(cli, "box",            NumType(0));
  // Set false for single-carrier (electrons only) runs — disables hole
  // particle type, inter-species scatter, and hole ESC/recombination.
  const bool USE_HOLES         = getArg(cli, "use_holes",      true);

  // ---- Simulation box ---------------------------------------------------
  // Band-filling requires (464 nm)³ box so N_c ≥ 50 per k-slot.
  // Alternatively, --boxside overrides directly (in nm).
  const NumType boxSide =
      (boxSideCLI > NumType(0))
          ? boxSideCLI
          : (USE_BAND_FILLING ? NumType(464e-9) : NumType(100e-9));
  const NumType Vsim = boxSide * boxSide * boxSide;
  const std::array<NumType, 3> maxPos  = {boxSide, boxSide, boxSide};
  const std::array<NumType, 3> spacing = {boxSide/10, boxSide/10, boxSide/10};

  // ---- Material parameters (CLI-overridable) ---------------------------
  // Per-material spectroscopic presets (manuscript Table 1 lists the sources).
  // Fields: eps_hi, eps_lo, m_e, m_h, hw_LO [eV], tau_LO [s], E_gap [eV].
  //   pb      = MAPbI3   (Sendner2016 / Wright2016 / Umari2014)
  //   sn      = CsSnI3   (Huang2013)
  //   fapbi3  = FAPbI3   (Eperon2014 gap / NatCommun2018 LO / Galkowski2016 m*)
  //   cspbi3  = CsPbI3   (gamma-phase; dielectric and LO are approximate)
  //   cspbbr3 = CsPbBr3  (Iaru2021)
  // Override any individual value with --eps_hi, --mass_e, --tau_lo, ... .
  // The generic "perovskite-like" trilogy set (10 / 30 / 33 meV) remains
  // reachable via --eps_hi 10 --eps_lo 30 --omega_lo 0.033 for SI runs.
  const std::string matName = getArg<std::string>(cli, "material", "pb");
  struct MatPreset { NumType eps_hi, eps_lo, m_e, m_h, hw_lo, tau_lo, e_gap; };
  MatPreset mp;
  if      (matName == "sn")      mp = {NumType(7.0), NumType(26.0), NumType(0.10), NumType(0.10), NumType(0.0200), NumType(0.1e-12), NumType(1.30)};
  else if (matName == "fapbi3")  mp = {NumType(5.0), NumType(30.0), NumType(0.17), NumType(0.17), NumType(0.0107), NumType(2.0e-12), NumType(1.48)};
  else if (matName == "cspbi3")  mp = {NumType(5.0), NumType(18.0), NumType(0.20), NumType(0.25), NumType(0.0130), NumType(1.0e-12), NumType(1.70)};
  else if (matName == "cspbbr3") mp = {NumType(4.5), NumType(16.0), NumType(0.20), NumType(0.20), NumType(0.0200), NumType(1.0e-12), NumType(2.30)};
  else                            mp = {NumType(5.0), NumType(33.5), NumType(0.20), NumType(0.25), NumType(0.0115), NumType(2.0e-12), NumType(1.60)}; // pb (MAPbI3)

  const NumType eps_hi        = getArg(cli, "eps_hi",  mp.eps_hi);
  const NumType eps_lo        = getArg(cli, "eps_lo",  mp.eps_lo);
  const NumType relEffMassE   = getArg(cli, "mass_e",  mp.m_e);
  const NumType relEffMassH   = getArg(cli, "mass_h",  mp.m_h);
  const NumType phononEnergy  = getArg(cli, "omega_lo",mp.hw_lo);  // [eV]
  const NumType tauLO         = getArg(cli, "tau_lo",  mp.tau_lo); // [s]
  const NumType E_gap         = getArg(cli, "E_gap",   mp.e_gap);  // [eV]
  const NumType latticeTempK  = getArg(cli, "T_lat",   NumType(300));   // [K]
  const NumType carrierDensity= getArg(cli, "density", NumType(1e24));  // [m^-3]
  const NumType E_photon      = getArg(cli, "E_photon",NumType(3.1));   // [eV]
  const NumType totalTime     = getArg(cli, "total_time", NumType(10e-12)); // [s]
  const NumType tauAcoustic   = getArg(cli, "tau_ac",  NumType(30e-12)); // [s]
  // Klemens: ħω_ac ≈ ħω_LO/2, applied per LO mode below.

  // ---- Model extensions (SI robustness studies; defaults reproduce the
  //      single-mode, Klemens-only baseline exactly) ---------------------
  // Multi-mode Froehlich: two MAPbI3 LO branches from Sendner 2016.
  const bool USE_MULTIMODE   = getArg(cli, "use_multimode", false);
  // Ridley channel LO -> LA + TO, branching ratio w_R (0 = Klemens only).
  const NumType ridleyW      = getArg(cli, "ridley_w", NumType(0));
  // Daughter TO mode: MAPbI3 TO branches sit at 32 and 63 cm^-1 (Sendner
  // 2016); 63 cm^-1 = 7.81 meV is the stretch-derived branch.
  const NumType toPhononEnergy = getArg(cli, "omega_to", NumType(0.00781)); // [eV]
  const NumType tauTO          = getArg(cli, "tau_to",   tauAcoustic);      // [s]

  std::cout << "=== Hot Carrier Cooling in MHP"
            << (USE_HOT_PHONON ? " (HPB enabled)" : " (equilibrium phonons)")
            << " ===\n\n";

  // MHP-like material: use eps_lo as the dielectric constant of the device
  // (used only for Debye length; irrelevant for bulk / periodic BC).
  MaterialType material(eps_lo,    // relative dielectric constant
                        4000.,     // mass density [kg/m^3] (perovskite ~4 g/cm^3)
                        1e10,      // intrinsic carrier density [m^{-3}] (wide gap)
                        2000.,     // sound velocity [m/s]
                        E_gap);    // band gap [eV]

  // Photoexcitation energies: excess split by inverse mass ratio
  const NumType E_excess = E_photon - material.getBandGap();
  const NumType fh = relEffMassH / (relEffMassE + relEffMassH); // m_h / M
  const NumType fe = relEffMassE / (relEffMassE + relEffMassH); // m_e / M
  const NumType E_e = fh * E_excess; // electron kinetic energy above CBM [eV]
  const NumType E_h = fe * E_excess; // hole kinetic energy above VBM [eV]
  std::cout << "Photoexcitation: E_photon = " << E_photon << " eV"
            << "  E_excess = " << E_excess << " eV"
            << "  E_e = " << E_e << " eV"
            << "  E_h = " << E_h << " eV\n\n";

  // Device temperature only affects Debye length; use lattice temperature.
  DeviceType device{material, maxPos, spacing, latticeTempK};
  device.addConstantDopingRegion({0, 0, 0}, maxPos, carrierDensity);

  // ---- Particle types: conduction band electrons + valence band holes ----
  MapIdxTypeToPartType particleTypes;

  // idx 0: electrons — initialise at fixed photoexcitation energy E_e
  particleTypes[0] = std::make_unique<emcElectron<NumType, DeviceType>>(
      1000, 4., /*usePotentialForInit=*/false, E_e);
  particleTypes[0]->addValley(
      std::make_unique<emcNonParabolicIsotropValley<NumType>>(
          relEffMassE, constants::me, /*degFactor=*/1, alpha_np));

  // idx 1: holes (skipped when USE_HOLES = false → single-carrier mode)
  if (USE_HOLES) {
    particleTypes[1] = std::make_unique<emcHole<NumType, DeviceType>>(
        1000, 4., /*usePotentialForInit=*/false, E_h);
    particleTypes[1]->addValley(
        std::make_unique<emcNonParabolicIsotropValley<NumType>>(
            relEffMassH, constants::me, /*degFactor=*/1, alpha_np));
  }

  // ---- Phonon bath(s) --------------------------------------------------
  // Both electrons and holes emit into the same LO phonon bath(s).
  // Optionally coupled to an acoustic bath (Klemens cascade, USE_ACOUSTIC_BATH)
  // and a transverse-optical reservoir (Ridley channel, --ridley_w > 0).
  //
  // Multi-mode Froehlich (--use_multimode 1): MAPbI3 has two IR-active LO
  // branches (Sendner 2016, Table 1): 40 cm^-1 (4.96 meV, Pb-I-Pb rocking)
  // and 133 cm^-1 (16.49 meV, Pb-I stretch). The ionic screening
  // 1/eps* = 1/eps_inf - 1/eps_s is partitioned between them by the residues
  // of 1/eps(w) at the LO poles, giving weights 6.2% and 93.8%. Each mode
  // gets its own bath and its own (eps_hi, eps_lo_n) pair reproducing that
  // partial screening exactly; summed they recover Sendner's 1/eps* = 0.17.
  const bool enableAcoustic = USE_HOT_PHONON && USE_ACOUSTIC_BATH;

  std::vector<NumType> modeEnergy, modeEpsLo;
  if (USE_MULTIMODE) {
    if (matName != "pb") {
      std::cerr << "ERROR: --use_multimode uses the MAPbI3 two-mode "
                   "decomposition (Sendner 2016) and is not defined for "
                   "other materials. Aborting.\n";
      return 1;
    }
    modeEnergy = {NumType(0.004959), NumType(0.016490)};  // [eV]
    modeEpsLo  = {NumType(5.2799),   NumType(25.4312)};   // with eps_hi = 5.0
  } else {
    modeEnergy = {phononEnergy};
    modeEpsLo  = {eps_lo};
  }
  const SizeType nrModes = modeEnergy.size();

  std::vector<std::shared_ptr<emcPhononBath<NumType>>> phononBaths;
  for (SizeType m = 0; m < nrModes; ++m) {
    phononBaths.push_back(std::make_shared<emcPhononBath<NumType>>(
        nrPhononBins, dq, tauLO, modeEnergy[m], latticeTempK, Vsim,
        enableAcoustic, modeEnergy[m] / NumType(2), tauAcoustic,
        ridleyW, toPhononEnergy, tauTO));
  }
  auto phononBath = phononBaths[0]; // primary bath (single-mode: the only one)

  if (USE_MULTIMODE)
    std::cout << "Multi-mode Froehlich: " << nrModes << " LO branches "
              << "(Sendner 2016 MAPbI3)\n";
  for (SizeType m = 0; m < nrModes; ++m)
    std::cout << "Phonon bath " << m << ": ħω_LO = " << modeEnergy[m] * 1e3
              << " meV,  eps_lo = " << modeEpsLo[m]
              << ",  N_0 = " << phononBaths[m]->getN0() << "\n";
  if (enableAcoustic)
    std::cout << "  Acoustic bath ON: τ_ac = " << tauAcoustic * 1e12 << " ps\n";
  if (phononBath->ridleyEnabled())
    std::cout << "  Ridley channel ON: w_R = " << ridleyW
              << ",  ħω_TO = " << toPhononEnergy * 1e3 << " meV"
              << ",  τ_TO = " << tauTO * 1e12 << " ps\n";

  // ---- Scatter mechanisms: Froehlich absorption + emission for e and h --
  if (USE_HOT_PHONON) {
    for (SizeType m = 0; m < nrModes; ++m) {
      const NumType wLO = modeEnergy[m];
      const NumType eLo = modeEpsLo[m];
      particleTypes[0]->addScatterMechanism(
          {0}, std::make_unique<emcHotPhononFroehlichAbsorption3D<NumType>>(
                   0, wLO, relEffMassE, eps_hi, eLo, phononBaths[m],
                   "MHP-e" + std::to_string(m)));
      particleTypes[0]->addScatterMechanism(
          {0}, std::make_unique<emcHotPhononFroehlichEmission3D<NumType>>(
                   0, wLO, relEffMassE, eps_hi, eLo, phononBaths[m],
                   "MHP-e" + std::to_string(m)));
      if (USE_HOLES) {
        particleTypes[1]->addScatterMechanism(
            {0}, std::make_unique<emcHotPhononFroehlichAbsorption3D<NumType>>(
                     0, wLO, relEffMassH, eps_hi, eLo, phononBaths[m],
                     "MHP-h" + std::to_string(m)));
        particleTypes[1]->addScatterMechanism(
            {0}, std::make_unique<emcHotPhononFroehlichEmission3D<NumType>>(
                     0, wLO, relEffMassH, eps_hi, eLo, phononBaths[m],
                     "MHP-h" + std::to_string(m)));
      }
    }
  } else {
    particleTypes[0]->addScatterMechanism(
        {0}, std::make_unique<emcFroehlichAbsorption3D<NumType>>(
                 0, phononEnergy, relEffMassE, eps_hi, eps_lo, latticeTempK,
                 "MHP-e"));
    particleTypes[0]->addScatterMechanism(
        {0}, std::make_unique<emcFroehlichEmission3D<NumType>>(
                 0, phononEnergy, relEffMassE, eps_hi, eps_lo, latticeTempK,
                 "MHP-e"));
    if (USE_HOLES) {
      particleTypes[1]->addScatterMechanism(
          {0}, std::make_unique<emcFroehlichAbsorption3D<NumType>>(
                   0, phononEnergy, relEffMassH, eps_hi, eps_lo, latticeTempK,
                   "MHP-h"));
      particleTypes[1]->addScatterMechanism(
          {0}, std::make_unique<emcFroehlichEmission3D<NumType>>(
                   0, phononEnergy, relEffMassH, eps_hi, eps_lo, latticeTempK,
                   "MHP-h"));
    }
  }

  // ---- Carrier-carrier scattering (Debye-screened binary collisions) ----
  emcCarrierCarrierScatter<NumType> ccScatterE(eps_hi, relEffMassE, Vsim,
                                               latticeTempK);
  emcCarrierCarrierScatter<NumType> ccScatterH(eps_hi, relEffMassH, Vsim,
                                               latticeTempK);

  // Inter-species e-h scattering (bipolar only)
  emcInterCarrierScatter<NumType> ccScatterEH(eps_hi, relEffMassE, relEffMassH,
                                              Vsim, latticeTempK);

  // ---- Band-to-band recombination (radiative + Auger) -----------------
  // MAPbI3 literature values (SI units):
  //   B   ~ 1e-16  m^3/s   (radiative;  = 1e-10  cm^3/s, Stranks 2014)
  //   C   ~ 1e-39  m^6/s   (Auger;      = 1e-27  cm^6/s, Herz 2016)
  //   E_g ~ 1.6 eV
  // At n = p = 1e24 m^-3 both lifetimes are ~ns — beyond the 10 ps cooling
  // window, so recombination barely affects the carrier count here. The
  // implementation is provided so longer simulations give correct loss rates.
  const NumType B_rad    = 1e-16; // radiative coefficient [m^3/s]
  const NumType C_auger  = 1e-39; // Auger coefficient [m^6/s]
  emcRecombination<NumType> recombination(B_rad, C_auger, C_auger, E_gap,
                                          relEffMassE, relEffMassH, Vsim);

  // ---- Energy-selective contacts (ESC) --------------------------------
  // E_ex: extraction energy above the band edge [eV]. At T_c ~ 3000 K the
  //       mean carrier energy is ~0.4 eV, so E_ex = 0.4 eV targets the bulk
  //       of the hot distribution.
  // deltaE: full width of the transmission window [eV]. Narrower windows give
  //         a higher extracted chemical potential at the cost of lower current.
  // tauEx: extraction time constant [s]. A 100 nm device with carriers at
  //        v ~ 10^5 m/s has a transit time of ~1 ps; use that as the baseline.
  const NumType E_ex_e  = getArg(cli, "E_ex_e",  getArg(cli, "E_ex", NumType(0.4)));
  const NumType E_ex_h  = getArg(cli, "E_ex_h",  getArg(cli, "E_ex", NumType(0.4)));
  const NumType deltaE  = getArg(cli, "delta_E",  NumType(0.05));
  const NumType tauEx   = getArg(cli, "tau_ex",   NumType(1e-12));
  emcEnergySelectiveContact<NumType> escElectrons(E_ex_e, deltaE, tauEx);
  emcEnergySelectiveContact<NumType> escHoles(E_ex_h, deltaE, tauEx);

  // ---- Particle handler (bulk, periodic BC, no applied field) ---------
  std::array<NumType, 3> noField = {0., 0., 0.};
  // --seed 0 (default) draws from the wall clock; non-zero makes the run
  // bit-reproducible, which the regression tests and SI comparisons rely on.
  const auto rngSeed =
      static_cast<long unsigned int>(getArg(cli, "seed", NumType(0)));
  ParticleHandler handler(device, particleTypes, noField, 0., rngSeed);

  // Generate initial hot carrier distribution
  std::cout << "Creating hot carrier ensemble...\n";
  handler.generateInitialParticles();
  handler.printNrParticles();
  // True initial count — the yield denominator. Remaining + extracted at
  // end-of-run undercounts it because recombined carriers vanish from both.
  const SizeType nrInitElectrons = handler.getNrParticles(0);

  // ---- Pauli exclusion grids (band filling) ----------------------------
  // Instantiated unconditionally; only used when USE_BAND_FILLING = true.
  emcPauliExclusion<NumType> pauliE(dk_pauli, kMax_pauli, Vsim);
  emcPauliExclusion<NumType> pauliH(dk_pauli, kMax_pauli, Vsim);
  if (USE_BAND_FILLING)
    std::cout << "Band-filling (Pauli exclusion) enabled — sequential move loop.\n";

  // ---- Main simulation loop -------------------------------------------
  const SizeType nrSteps = static_cast<SizeType>(std::ceil(totalTime / dt));
  std::cout << "Simulation: " << nrSteps << " steps × " << dt * 1e15
            << " fs = " << totalTime * 1e12 << " ps\n\n";

  // Coupling-weighted mean LO occupation across modes. Weight = partial ionic
  // screening 1/eps_hi - 1/eps_lo_n, i.e. each mode's share of the Froehlich
  // interaction. Single-mode: reduces to that bath's getMeanNq().
  std::vector<NumType> modeWeight(nrModes);
  {
    NumType wSum = 0;
    for (SizeType m = 0; m < nrModes; ++m) {
      modeWeight[m] = NumType(1) / eps_hi - NumType(1) / modeEpsLo[m];
      wSum += modeWeight[m];
    }
    for (auto &w : modeWeight) w /= wSum;
  }
  auto meanNqWeighted = [&]() {
    NumType s = 0;
    for (SizeType m = 0; m < nrModes; ++m)
      s += modeWeight[m] * phononBaths[m]->getMeanNq();
    return s;
  };

  std::string suffix = USE_HOT_PHONON ? "HPB" : "EQ";
  if (enableAcoustic)    suffix += "_AC";
  if (USE_MULTIMODE)     suffix += "_MM";
  if (ridleyW > 0)       suffix += "_RID";
  if (USE_BAND_FILLING)  suffix += "_BF";
  if (!USE_HOLES)        suffix += "_1C";  // single-carrier
  if (USE_ESC)           suffix += "_ESC";
  std::ofstream energyFileE("avgEnergyElectrons" + suffix + ".txt");
  std::ofstream energyFileH("avgEnergyHoles" + suffix + ".txt");
  std::ofstream phononFile("phononOccupation" + suffix + ".txt");
  std::ofstream carrierFile("nrCarriers" + suffix + ".txt");
  std::ofstream tempFile("carrierTemp" + suffix + ".txt");
  // Per-step FD fit columns enable the extraction-weighted (working) V_OC to
  // be reconstructed in post-processing: the voltage is built from the carrier
  // temperature and chemical potential at the moment carriers are extracted,
  // weighted by the extraction rate, rather than from the end-of-run remnant.
  tempFile << "# time[s]  T_MB_e[K]  T_MB_h[K]  T_FD_e[K]  mu_e[eV]  "
              "T_FD_h[K]  mu_h[eV]\n";
  carrierFile << "# time[s]  N_electrons  N_holes  N_esc_e  N_esc_h\n";
  // Phonon file: always write N_LO (coupling-weighted across modes); add
  // N_ac/T_ac when the acoustic bath is on and N_TO/T_TO when Ridley is on.
  {
    phononFile << "# time[s]  N_LO";
    if (enableAcoustic)               phononFile << "  N_ac  T_ac[K]";
    if (phononBath->ridleyEnabled())  phononFile << "  N_TO  T_TO[K]";
    phononFile << "\n";
  }

  // Record initial state
  auto avgE0e = handler.getAvgEnergy(0);
  energyFileE << 0.0 << " " << avgE0e[0] << "\n";
  if (USE_HOLES) {
    auto avgE0h = handler.getAvgEnergy(1);
    energyFileH << 0.0 << " " << avgE0h[0] << "\n";
  }
  phononFile  << 0.0 << " " << meanNqWeighted();
  if (enableAcoustic)
    phononFile << " " << phononBath->getMeanNac()
               << " " << phononBath->getAcousticTemp();
  if (phononBath->ridleyEnabled())
    phononFile << " " << phononBath->getMeanNTO()
               << " " << phononBath->getTOTemp();
  phononFile << "\n";
  carrierFile << 0.0
              << " " << handler.getNrParticles(0)
              << " " << (USE_HOLES ? handler.getNrParticles(1) : SizeType(0))
              << " 0 0\n";

  auto tStart = std::chrono::high_resolution_clock::now();

  for (SizeType step = 1; step <= nrSteps; step++) {
    // 1. Drift + Froehlich scatter (+ Pauli exclusion when USE_BAND_FILLING)
    if (USE_BAND_FILLING) {
      handler.moveParticleTypeWithBandFilling(dt, 0, pauliE);
      if (USE_HOLES)
        handler.moveParticleTypeWithBandFilling(dt, 1, pauliH);
    } else {
      handler.moveParticles(dt);
    }

    // 2. Carrier-carrier thermalisation
    if (USE_CC) {
      handler.carrierCarrierScatter(ccScatterE, 0, dt);
      if (USE_HOLES) {
        handler.carrierCarrierScatter(ccScatterH, 1, dt);
        handler.interCarrierScatter(ccScatterEH, 0, 1, dt);
      }
    }

    // 3. Band-to-band recombination (bipolar only)
    if (USE_HOLES && USE_RECOMB)
      handler.recombine(recombination, 0, 1, dt);

    // 4. Energy-selective contact extraction
    if (USE_ESC) {
      handler.extractCarriers(escElectrons, 0, dt);
      if (USE_HOLES)
        handler.extractCarriers(escHoles, 1, dt);
    }

    // 5. Evolve phonon bath(s) and rebuild scatter tables (HPB only)
    if (USE_HOT_PHONON) {
      for (auto &bath : phononBaths)
        bath->update(dt);
      particleTypes[0]->reinitScatterTables();
      if (USE_HOLES)
        particleTypes[1]->reinitScatterTables();
    }

    if (step % nrStepsBetweenOutput == 0) {
      auto avgEe = handler.getAvgEnergy(0);
      NumType t = step * dt;
      energyFileE << t << " " << avgEe[0] << "\n";
      NumType T_MBe = emcHotCarrierOutput::getMBTemp(avgEe[0]);
      NumType T_MBh = T_MBe;
      NumType avgEh_val = avgEe[0];
      if (USE_HOLES) {
        auto avgEh = handler.getAvgEnergy(1);
        avgEh_val = avgEh[0];
        energyFileH << t << " " << avgEh_val << "\n";
        T_MBh = emcHotCarrierOutput::getMBTemp(avgEh_val);
      }

      // Per-step Fermi-Dirac fit on the live population (for extraction-
      // weighted V_OC in post-processing). Cheap: bisection, every 0.05 ps.
      NumType n_e_t = static_cast<NumType>(handler.getNrParticles(0)) / Vsim;
      auto fdEt = emcHotCarrierOutput::getFDFit(avgEe[0], n_e_t, relEffMassE,
                                                latticeTempK);
      NumType T_FDe_t = fdEt.T_FD, mu_e_t = fdEt.mu_eV;
      NumType T_FDh_t = T_FDe_t, mu_h_t = mu_e_t;
      if (USE_HOLES && handler.getNrParticles(1) > 0) {
        NumType n_h_t = static_cast<NumType>(handler.getNrParticles(1)) / Vsim;
        auto fdHt = emcHotCarrierOutput::getFDFit(avgEh_val, n_h_t, relEffMassH,
                                                  latticeTempK);
        T_FDh_t = fdHt.T_FD;
        mu_h_t  = fdHt.mu_eV;
      }
      phononFile  << t << " " << meanNqWeighted();
      if (enableAcoustic)
        phononFile << " " << phononBath->getMeanNac()
                   << " " << phononBath->getAcousticTemp();
      if (phononBath->ridleyEnabled())
        phononFile << " " << phononBath->getMeanNTO()
                   << " " << phononBath->getTOTemp();
      phononFile << "\n";
      carrierFile << t
                  << " " << handler.getNrParticles(0)
                  << " " << (USE_HOLES ? handler.getNrParticles(1) : SizeType(0))
                  << " " << escElectrons.getNrExtracted()
                  << " " << (USE_HOLES ? escHoles.getNrExtracted() : SizeType(0))
                  << "\n";
      tempFile << t << " " << T_MBe << " " << T_MBh
               << " " << T_FDe_t << " " << mu_e_t
               << " " << T_FDh_t << " " << mu_h_t << "\n";

      std::cout << "  t = " << t * 1e12 << " ps"
                << " | T_e = " << static_cast<int>(T_MBe) << " K";
      if (USE_HOLES)
        std::cout << " | T_h = " << static_cast<int>(T_MBh) << " K";
      std::cout << " | N_LO = "  << meanNqWeighted();
      if (enableAcoustic)
        std::cout << " | T_ac = " << static_cast<int>(phononBath->getAcousticTemp()) << " K";
      if (USE_BAND_FILLING)
        std::cout << " | Pauli_rej = "
                  << static_cast<int>(pauliE.getRejectionRate() * 100) << "%e";
      if (USE_HOLES && USE_ESC)
        std::cout << " | ESC: " << escElectrons.getNrExtracted() << "e/"
                  << escHoles.getNrExtracted() << "h";
      else if (USE_ESC)
        std::cout << " | ESC: " << escElectrons.getNrExtracted() << "e";
      std::cout << "\n";
    }
  }

  auto tEnd = std::chrono::high_resolution_clock::now();
  std::cout << "\nCPU time: "
            << std::chrono::duration_cast<std::chrono::seconds>(tEnd - tStart)
                   .count()
            << " s\n";

  energyFileE.close();
  energyFileH.close();
  phononFile.close();
  carrierFile.close();
  tempFile.close();

  // ---- Final carrier temperature and VOC estimate ---------------------
  {
    auto finalEe = handler.getAvgEnergy(0);
    NumType n_e = static_cast<NumType>(handler.getNrParticles(0)) / Vsim;
    NumType T_MBe = emcHotCarrierOutput::getMBTemp(finalEe[0]);
    std::cout << "\n── Final carrier temperatures ──────────────────────────\n"
              << "  T_MB (electrons) = " << static_cast<int>(T_MBe) << " K";
    NumType n_h = NumType(0);
    NumType T_MBh = T_MBe;
    NumType finalEh_avg = NumType(0);
    if (USE_HOLES) {
      auto finalEh = handler.getAvgEnergy(1);
      finalEh_avg = finalEh[0];
      n_h = static_cast<NumType>(handler.getNrParticles(1)) / Vsim;
      T_MBh = emcHotCarrierOutput::getMBTemp(finalEh_avg);
      std::cout << "  T_MB (holes) = " << static_cast<int>(T_MBh) << " K";
    }
    std::cout << "\n";

    if (USE_HOLES && n_e > 0 && n_h > 0) {
      auto fdE = emcHotCarrierOutput::getFDFit(finalEe[0], n_e, relEffMassE, latticeTempK);
      auto fdH = emcHotCarrierOutput::getFDFit(finalEh_avg, n_h, relEffMassH, latticeTempK);

      std::cout << "  FD fit (electrons): T = " << static_cast<int>(fdE.T_FD)
                << " K  μ = " << fdE.mu_eV << " eV  η = " << fdE.eta
                << (fdE.converged ? "" : " [NOT CONVERGED]") << "\n"
                << "  FD fit (holes):     T = " << static_cast<int>(fdH.T_FD)
                << " K  μ = " << fdH.mu_eV << " eV  η = " << fdH.eta
                << (fdH.converged ? "" : " [NOT CONVERGED]") << "\n";

      NumType T_avg = NumType(0.5) * (fdE.T_FD + fdH.T_FD);
      NumType voc = emcHotCarrierOutput::getVOC(
          fdE.mu_eV, fdH.mu_eV, fdE.T_FD, fdH.T_FD, latticeTempK,
          E_gap, finalEe[0], finalEh_avg);

      // Yield-adjusted J_SC: fraction of carriers actually extracted × 25 mA/cm²
      // When ESC is off, assume 100% yield (theoretical maximum).
      NumType n_esc_e = static_cast<NumType>(escElectrons.getNrExtracted());
      NumType n_init_e = static_cast<NumType>(nrInitElectrons);
      NumType yield_e = (USE_ESC && n_init_e > 0) ? n_esc_e / n_init_e : NumType(1);
      NumType J_eff = NumType(25) * yield_e;
      NumType pce = emcHotCarrierOutput::getPCE(J_eff, voc);

      std::cout << "\n── Hot-carrier VOC estimate (Paper 3 eq. 3) ────────────\n"
                << "  T_eh = " << static_cast<int>(T_avg) << " K"
                << "  (T_L = " << static_cast<int>(latticeTempK) << " K)\n"
                << "  V_OC = " << voc << " V\n"
                << "  yield_e = " << yield_e * 100 << " %\n"
                << "  PCE  = " << pce * 100 << " %"
                << "  (J_eff = " << J_eff << " mA/cm², FF = 0.85, AM1.5G)\n";
    }
  }

  if (USE_ESC) {
    NumType area = boxSide * boxSide;
    std::cout << "\n── ESC summary ─────────────────────────────────────────\n"
              << "  Electrons: " << escElectrons.getNrExtracted()
              << " extracted  |  <E_ex> = "
              << escElectrons.getMeanExtractedEnergy() << " eV"
              << "  |  J_e = "
              << escElectrons.getCurrentDensity(area, totalTime) * 1e-4
              << " A/cm^2\n"
              << "  Holes:     " << escHoles.getNrExtracted()
              << " extracted  |  <E_ex> = "
              << escHoles.getMeanExtractedEnergy() << " eV"
              << "  |  J_h = "
              << escHoles.getCurrentDensity(area, totalTime) * 1e-4
              << " A/cm^2\n"
              << "  (ESC window: [E_ex - " << deltaE / 2 << ", E_ex + "
              << deltaE / 2 << "] eV, tau_ex = " << tauEx * 1e12 << " ps)\n";
  }

  std::cout << "\nOutput: avgEnergyElectrons" << suffix << ".txt"
            << ", avgEnergyHoles" << suffix << ".txt"
            << ", phononOccupation" << suffix << ".txt"
            << ", nrCarriers" << suffix << ".txt"
            << ", carrierTemp" << suffix << ".txt\n";
  return 0;
}
