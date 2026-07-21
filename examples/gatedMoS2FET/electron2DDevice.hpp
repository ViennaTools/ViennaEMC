#ifndef ELECTRON_2D_DEVICE_HPP
#define ELECTRON_2D_DEVICE_HPP

#include <ParticleType/emcParticleType.hpp>
#include <emcConstants.hpp>

/**
 * @brief Injected, device-compatible 2D-material electron.
 *
 * Combines the 2D transport of examples/singleLayerMoS2/electron2D.hpp (wave
 * vector in the in-plane x-y directions, k_z = 0 => carriers confined to the
 * monolayer) with the injection / density-from-potential machinery of
 * examples/mosfet2D/electronVWD.hpp, so it can be used in the self-consistent
 * emcSimulation (EMC + Poisson) with ohmic + gate contacts.
 */
template <class T, class DeviceType>
struct electron2DDevice : public emcParticleType<T, DeviceType> {
  typedef typename DeviceType::ValueVec ValueVec;
  typedef typename DeviceType::SizeVec SizeVec;
  static const SizeType Dim = DeviceType::Dimension;

  std::uniform_real_distribution<T> dist{0., 1.};
  std::uniform_real_distribution<T> distForLog{1e-6, 1.};

  //! 2D Fermi-Dirac degeneracy parameter (quantum capacitance). 0 = Boltzmann.
  //! MUST equal the value passed to emcSORSolver::setDegeneracyBeta so the
  //! particle density and the Poisson closure use the same statistics.
  T degeneracyBeta = 0;
  void setDegeneracyBeta(T beta) { degeneracyBeta = beta; }

  electron2DDevice(SizeType nrEnergyLevels = 1000, T maxEnergy = 4.)
      : emcParticleType<T, DeviceType>(nrEnergyLevels, maxEnergy) {}

  std::string getName() const { return "Electrons"; }
  T getMass() const { return constants::me; }
  T getCharge() const { return -constants::q; }
  bool isMoved() const { return true; }
  bool isInjected() const { return true; }

  //! normalized 2D-material electron volume density n(phi) [1/m^3] from the
  //! (normalized) potential: Boltzmann (beta<=0) or degenerate 2D Fermi-Dirac
  //! (beta>0, quantum capacitance). Shared by the equilibrium + Schottky paths.
  T eVolumeDensity(T phi, T Ni) const {
    if (degeneracyBeta <= 0)
      return std::exp(phi) * Ni;
    T e = std::exp(phi);
    return std::log1p(degeneracyBeta * e) / degeneracyBeta * Ni;
  }

  //! equilibrium electron number from the self-consistent potential
  T getInitialNrParticles(const SizeVec &coord, const DeviceType &device,
                          const emcGrid<T, Dim> &potential) {
    T eDensity = eVolumeDensity(potential[coord], device.getMaterial().getNi());
    for (SizeType d = 0; d < Dim; d++)
      if (coord[d] == 0 || coord[d] == potential.getSize(d) - 1)
        eDensity *= 0.5;
    return std::round(eDensity * device.getCellVolume());
  }

  //! expected reservoir number at a contact: from the doping at an OHMIC
  //! contact, or from the metal barrier at a SCHOTTKY contact (n_res =
  //! Ni*exp(phi_sb), phi_sb = -ln(beta) - phiB/Vth -- matching the solver BC).
  T getExpectedNrParticlesAtContact(const SizeVec &coord,
                                    const DeviceType &device) {
    const auto &surface = device.getSurface();
    T nr;
    if (surface.isSchottkyContact(coord)) {
      auto boundPos = surface.getBoundaryPos(coord);
      auto coordSurf = surface.getCoordBoundary(coord, boundPos);
      T barrier = surface.getContactFurtherParameter(coordSurf, boundPos, 0);
      T phiSb = -std::log(degeneracyBeta) - device.normalizeVoltage(barrier);
      nr = device.getCellVolume() *
           eVolumeDensity(phiSb, device.getMaterial().getNi());
    } else {
      nr = device.getCellVolume() * device.getDopingProfile().getDoping(coord);
    }
    const auto &extent = device.getGridExtent();
    for (SizeType d = 0; d < Dim; d++)
      if (coord[d] == 0 || coord[d] == extent[d] - 1)
        nr *= 0.5;
    return std::round(nr);
  }

  emcParticle<T> generateInitialParticle(const SizeVec &coord,
                                         const DeviceType &device, emcRNG &rng) {
    return make2DParticle(coord, device, rng);
  }
  emcParticle<T> generateInjectedParticle(const SizeVec &coord,
                                          const DeviceType &device,
                                          emcRNG &rng) {
    return make2DParticle(coord, device, rng);
  }

private:
  //! 2D in-plane wave vector (k_z = 0), thermal energy, region-based tau.
  emcParticle<T> make2DParticle(const SizeVec &coord, const DeviceType &device,
                                emcRNG &rng) {
    emcParticle<T> part;
    part.region = device.getDopingProfile().getDopingRegionIdx(coord);
    part.valley = std::floor(this->getNrValleys() * dist(rng));
    auto &valley = this->valleys[part.valley];
    part.subValley = std::floor(valley->getDegeneracyFactor() * dist(rng));
    part.energy = -device.getThermalVoltage() * std::log(distForLog(rng));
    auto normK = valley->getNormWaveVec(part.energy);
    T angle = 2 * constants::pi * dist(rng);
    part.k[0] = normK * std::cos(angle);
    part.k[1] = normK * std::sin(angle);
    part.k[2] = 0; // monolayer confinement
    part.tau = this->getNewTau(part.valley, part.region, rng);
    part.grainTau = this->getNewGrainTau(rng);
    return part;
  }
};

#endif // ELECTRON_2D_DEVICE_HPP
