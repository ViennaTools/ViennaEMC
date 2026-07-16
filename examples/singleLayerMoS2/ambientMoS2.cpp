// ambientMoS2 — runs the single-layer MoS2 ambient/adsorbate model (layer C).
//
// Three demonstrations, all through the MoS2Ambient model:
//   (1) O2 partial-pressure sweep  -> sigma(P)/sigma0   (C1 depletion + C2 Coulomb)
//   (2) relative-humidity sweep    -> mu(RH)            (C3 dielectric screening)
//   (3) adsorption/desorption noise -> N_ads(t) series  (C4, 1/f + RTN)
//
// Mobility is obtained the same way as the other MoS2 examples: a short bulk EMC
// drift-velocity run at a low field. sigma = n e mu. Outputs three text files that
// the plot scripts in validation/ (plot_ambient.py, plot_noise.py) can read.
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

#include <emcDevice.hpp>
#include <emcUtil.hpp>

#include "../bulkSimulation/basicBulkParticleHandler.hpp"
#include "electron2D.hpp"
#include "parameterAmbient.hpp"

using NumType = double;
const SizeType Dim = 3;
using DeviceType = emcDevice<NumType, Dim>;
using ParticleHandler = basicBulkParticleHandler<NumType, DeviceType>;
using MapIdxTypeToPartType = ParticleHandler::MapIdxToParticleTypes;
using ValueVec = DeviceType::ValueVec;

// --- simulation knobs ---------------------------------------------------------
const NumType temperature = 300;   // K
const NumType dt = 1e-16;          // s   (EMC time step)
const NumType field = 3e5;         // V/m (low field, ohmic regime)
const ValueVec fieldDir = {1, 0, 0};
const ValueVec maxPos = {5e-7, 5e-7, 0.65e-9};
const ValueVec spacing = {1e-8, 1e-8, 0.65e-9};
const SizeType nrSteps = 25000, nrAvgFrom = 12500;

// baseline charged-impurity density for the humidity demo (fixed substrate
// charges being screened by the adsorbed-water dielectric), 5e12 cm^-2, each of
// unit charge (Z=1), unlike the partial-charge O2 acceptors.
const NumType baselineImpurityDensity = 5e16;
const NumType substrateImpurityCharge = 1.0;

const std::string prefix = "ambientMoS2";

// Low-field EMC mobility [cm^2/Vs] for a given ambient state, via the model.
NumType mobility(NumType n, NumType nImp, NumType epsAvg,
                 NumType impCharge = MoS2Ambient::scattererCharge) {
  emcMaterial<NumType> MoS2{1, 1, 1, 1, 1};
  DeviceType device{MoS2, maxPos, spacing};
  device.addConstantDopingRegion({0, 0, 0}, maxPos, 1);
  MapIdxTypeToPartType pts;
  pts[0] = std::make_unique<electron2D<NumType, DeviceType>>();
  MoS2Ambient::addAmbientScatterMechanisms(pts[0], temperature, n, nImp, epsAvg,
                                           impCharge);
  ParticleHandler handler(device, pts, fieldDir, field, 12345);
  handler.generateInitialParticles();
  NumType sum = 0;
  SizeType nAvg = 0;
  for (SizeType s = 1; s <= nrSteps; ++s) {
    handler.moveParticles(dt);
    if (s >= nrAvgFrom) {
      sum += handler.getAvgDriftVelocity(0)[0];
      ++nAvg;
    }
  }
  return std::fabs(sum / nAvg / field * 1e4); // m^2/Vs -> cm^2/Vs
}

// (1) O2 partial-pressure sweep: sigma(P) with C1 depletion + C2 scattering.
void runOxygenSweep() {
  std::cout << "\n=== (1) O2 sweep: C1 depletion + C2 Coulomb scattering ===\n";
  std::ofstream os(prefix + "_O2sweep.txt");
  os << "# P[Pa]   theta   n[cm^-2]   N_ads[cm^-2]   mu[cm2/Vs]   sigma/sigma0\n";
  std::cout << "# P[Pa]     theta   n[cm^-2]   mu[cm2/Vs]  sigma/sigma0\n";
  NumType epsAvg = MoS2Ambient::envScreeningPermittivity, sigma0 = 0;
  for (NumType P : {0.0, 1e3, 1e4, 3e4, 1e5}) {
    NumType n = MoS2Ambient::effectiveCarrierDensity(P);
    NumType nAds = MoS2Ambient::chargedScattererDensity(P);
    NumType mu = mobility(n, nAds, epsAvg);
    NumType sigma = n * constants::q * mu * 1e-4;
    if (P == 0.0)
      sigma0 = sigma;
    os << std::scientific << std::setprecision(2) << P << "  " << std::fixed
       << std::setprecision(3) << MoS2Ambient::coverage(P) << "  "
       << std::scientific << std::setprecision(2) << n / 1e4 << "  " << nAds / 1e4
       << "   " << std::fixed << std::setprecision(1) << mu << "        "
       << std::setprecision(3) << sigma / sigma0 << "\n";
    std::cout << "  " << std::scientific << std::setprecision(1) << P << "  "
              << std::fixed << std::setprecision(3) << MoS2Ambient::coverage(P)
              << "  " << std::scientific << std::setprecision(2) << n / 1e4
              << "  " << std::fixed << std::setprecision(1) << mu << "       "
              << std::setprecision(3) << sigma / sigma0 << "\n";
  }
  std::cout << "  wrote " << prefix << "_O2sweep.txt\n";
}

// (2) Relative-humidity sweep: mu(RH) via environmental dielectric screening.
void runHumiditySweep() {
  std::cout << "\n=== (2) Humidity sweep: C3 dielectric screening (EDS) ===\n";
  std::ofstream os(prefix + "_humidity.txt");
  os << "# RH   theta_w   eps_avg   mu[cm2/Vs]   mu/mu_dry\n";
  std::cout << "# RH    eps_avg   mu[cm2/Vs]   mu/mu_dry\n";
  NumType n = MoS2Ambient::baselineDensity, muDry = 0;
  for (NumType RH : {0.0, 0.25, 0.5, 0.75, 1.0}) {
    NumType epsAvg = MoS2Ambient::effectiveEnvPermittivity(RH);
    NumType mu = mobility(n, baselineImpurityDensity, epsAvg,
                          substrateImpurityCharge);
    if (RH == 0.0)
      muDry = mu;
    os << std::fixed << std::setprecision(2) << RH << "   "
       << MoS2Ambient::waterFilm().waterCoverage(RH) << "   " << std::setprecision(2)
       << epsAvg << "   " << std::setprecision(1) << mu << "        "
       << std::setprecision(3) << mu / muDry << "\n";
    std::cout << "  " << std::fixed << std::setprecision(2) << RH << "   "
              << epsAvg << "     " << std::setprecision(1) << mu << "       "
              << std::setprecision(3) << mu / muDry << "\n";
  }
  std::cout << "  wrote " << prefix << "_humidity.txt\n";
}

// (3) Adsorption/desorption conductivity noise (C4): N_ads(t) time series.
void runNoise() {
  std::cout << "\n=== (3) Adsorption/desorption noise: C4 (1/f + RTN) ===\n";
  const NumType dtN = 2e-6;           // s   (<< 1/gamma_max)
  const SizeType nStepsN = 300000;    // 0.6 s, fs = 5e5 Hz
  emcRNG rngGen(1234u);
  // 1/f: many sites, gamma log-uniform 30..6e4/s -> clean 1/f ~5 Hz..10 kHz.
  emcAdsorbateNoise<NumType> onef =
      MoS2Ambient::makeOnefNoiseModel(1500, 30.0, 60000.0, rngGen);
  // RTN: a single site, gamma = 200/s (corner ~32 Hz).
  emcAdsorbateNoise<NumType> rtn(std::vector<NumType>{100.0},
                                 std::vector<NumType>{100.0});
  emcRNG rngA(11u), rngB(22u);
  rtn.seedEquilibrium(rngA);
  onef.seedEquilibrium(rngB);
  std::ofstream os(prefix + "_noise.txt");
  os << dtN << " " << nStepsN << " " << onef.nSites() << "\n"; // header
  for (SizeType s = 0; s < nStepsN; ++s) {
    SizeType nR = rtn.step(dtN, rngA);
    SizeType nF = onef.step(dtN, rngB);
    os << nR << " " << nF << "\n";
  }
  std::cout << "  wrote " << prefix << "_noise.txt (" << nStepsN
            << " steps, RTN 1 site + 1/f " << onef.nSites() << " sites)\n";
}

int main() {
  std::cout << "single-layer MoS2 ambient model (layer C) — T = " << temperature
            << " K\n";
  auto t0 = std::chrono::high_resolution_clock::now();
  runOxygenSweep();
  runHumiditySweep();
  runNoise();
  auto t1 = std::chrono::high_resolution_clock::now();
  std::cout << "\nTotal time: "
            << std::chrono::duration_cast<std::chrono::seconds>(t1 - t0).count()
            << " s\n";
  return 0;
}
