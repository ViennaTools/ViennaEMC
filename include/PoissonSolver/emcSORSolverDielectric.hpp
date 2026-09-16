#ifndef EMC_SOR_SOLVER_DIELECTRIC_HPP
#define EMC_SOR_SOLVER_DIELECTRIC_HPP

#include <PoissonSolver/emcAbstractSolver.hpp>
#include <emcSurface.hpp>

#include <cmath>
#include <functional>
#include <numeric>

/*! \brief Successive over-relaxation Poisson solver for a semiconductor sheet
 * embedded in dielectrics (gate oxide, air).
 *
 * Additive alternative to emcSORSolver for devices whose grid contains regions
 * without mobile carriers around a thin semiconductor layer (a 2D-material
 * channel). Differences to emcSORSolver:
 *  - a per-node relative permittivity, given in units of the device material's
 *    eps_r (which sets the Debye-length normalisation); the face coefficients
 *    use the harmonic mean of the two node values;
 *  - the carrier closure (equilibrium) and the doping enter only at nodes flagged
 *    as semiconductor; dielectric nodes carry only the particle charge that the
 *    particle-mesh scheme assigns to them;
 *  - GATE contacts are Dirichlet planes (the oxide is meshed explicitly, so the
 *    oxide-capacitor boundary condition of emcSORSolver is not used) and the
 *    applied gate voltage is part of the equilibrium solve, so that the initial
 *    carrier distribution is the gated one;
 *  - OHMIC and SCHOTTKY contacts are handled as in emcSORSolver.
 * The FMM background-potential path is not supported. No existing class is
 * modified; models that do not instantiate this solver are unaffected.
 *
 * @param accuracy stopping criterion [V]: the solver stops when the largest
 * potential update of a sweep is below this value
 * @param omega relaxation factor in (0, 2)
 * @param epsRelOfNode relative permittivity of a node divided by the device
 * material's eps_r
 * @param isSemiconductorNode true for nodes that carry the carrier closure and
 * the doping (the sheet)
 */
template <class T, class DeviceType, class ParticleHandler>
class emcSORSolverDielectric
    : public emcAbstractSolver<T, DeviceType, ParticleHandler> {
  static const int Dim = DeviceType::Dimension;
  typedef typename DeviceType::SizeVec SizeVec;
  typedef typename DeviceType::SizeVecSurface SizeVecSurface;
  typedef typename DeviceType::ValueVec ValueVec;
  typedef emcGrid<T, Dim> GridType;

public:
  emcSORSolverDielectric(const DeviceType &inDevice, T inAccuracy, T inOmega,
                         std::function<T(const SizeVec &)> epsRelOfNode,
                         std::function<bool(const SizeVec &)> isSemiconductorNode)
      : accuracy(inDevice.normalizeVoltage(inAccuracy)), omega(inOmega),
        device(inDevice), surface(inDevice.getSurface()),
        doping(inDevice.getDopingProfile().getDoping(true)),
        epsRel(inDevice.getGridExtent(), 1), isSemi(inDevice.getGridExtent(), 0),
        isDir(inDevice.getGridExtent(), 0), h(inDevice.getSpacing(true)) {
    init(epsRelOfNode, isSemiconductorNode);
  }

  //! 2D Fermi-Dirac degeneracy parameter of the sheet nodes (0 = Boltzmann);
  //! beta = Ni * dz / N2D with dz the grid spacing normal to the sheet.
  void setDegeneracyBeta(T beta) { degeneracyBeta = beta; }
  //! include the applied gate voltage in the equilibrium solve (default true)
  void setGateVoltageInEquilibrium(bool on) { gateInEquilibrium = on; }
  void setMaxSweeps(SizeType n) { maxSweeps = n; }
  SizeType getLastNrSweeps() const { return lastSweeps; }

  void calcEquilibriumPotential(GridType &pot, const DeviceType & /*device*/,
                                bool resetBC = true) {
    if (resetBC)
      setDirichletValues(pot, false);
    solve(pot, nullptr);
  }

  void calcNonEquilibriumPotential(GridType &pot, const DeviceType & /*device*/,
                                   const GridType &eConc, bool resetBC = true) {
    if (resetBC)
      setDirichletValues(pot, true);
    solve(pot, &eConc);
  }

  void calcBackgroundPotential(GridType & /*pot*/, const DeviceType & /*device*/,
                               ParticleHandler & /*handler*/,
                               bool /*resetBC*/ = true) {
    emcMessage::getInstance()
        .addError("emcSORSolverDielectric does not support the FMM background "
                  "potential.")
        .print();
  }

private:
  T accuracy, omega;
  const DeviceType &device;
  const typename DeviceType::SurfaceType &surface;
  GridType doping, epsRel;
  emcGrid<int, Dim> isSemi, isDir;
  ValueVec h, hFactor;
  T hProduct = 1;
  SizeVec maxCoord;
  T degeneracyBeta = 0;
  T potStepClamp = 2.0; //!< max |update| per sweep [kT/q] for the degenerate closure
  bool gateInEquilibrium = true;
  SizeType maxSweeps = 500000, lastSweeps = 0;

  void init(const std::function<T(const SizeVec &)> &epsRelOfNode,
            const std::function<bool(const SizeVec &)> &isSemiconductorNode) {
    if (Dim == 2) {
      hFactor[0] = h[1] / h[0];
      hFactor[1] = h[0] / h[1];
    } else {
      hFactor[0] = h[1] * h[2] / h[0];
      hFactor[1] = h[0] * h[2] / h[1];
      hFactor[2] = h[0] * h[1] / h[2];
    }
    hProduct = std::accumulate(h.begin(), h.end(), T(1), std::multiplies<T>());
    maxCoord = device.getGridExtent();
    for (auto &val : maxCoord)
      val--;
    SizeVec coord;
    for (coord.fill(0); !epsRel.isEndCoord(coord); epsRel.advanceCoord(coord)) {
      epsRel[coord] = epsRelOfNode(coord);
      isSemi[coord] = isSemiconductorNode(coord) ? 1 : 0;
      isDir[coord] = (surface.isReservoirContact(coord) || isGateNode(coord)) ? 1 : 0;
    }
  }

  //! true if the node lies on a boundary face that carries a GATE contact
  bool isGateNode(const SizeVec &coord) const {
    for (SizeType d = 0; d < Dim; d++) {
      for (int side = 0; side < 2; side++) {
        if ((side == 0 && coord[d] != 0) || (side == 1 && coord[d] != maxCoord[d]))
          continue;
        auto boundPos = static_cast<emcBoundaryPos>(2 * d + side);
        auto coordSurf = surface.getCoordBoundary(coord, boundPos);
        if (surface.isContactType(coordSurf, boundPos, emcContactType::GATE))
          return true;
      }
    }
    return false;
  }

  //! Dirichlet values at reservoir (ohmic, Schottky) and gate contact nodes
  void setDirichletValues(GridType &pot, bool applied) {
    SizeVecSurface coordSurf;
    emcBoundaryPos boundPos;
    for (surface.initCoord(boundPos, coordSurf);
         !surface.isEndCoord(boundPos, coordSurf);
         surface.advanceCoord(boundPos, coordSurf)) {
      if (surface.isArtificialBoundary(coordSurf, boundPos))
        continue;
      auto coord = surface.getCoordDevice(coordSurf, boundPos);
      T vApp = surface.getContactVoltage(coordSurf, boundPos, true);
      if (surface.isContactType(coordSurf, boundPos, emcContactType::OHMIC)) {
        pot[coord] = builtInPotential(doping[coord]) + (applied ? vApp : 0);
      } else if (surface.isContactType(coordSurf, boundPos, emcContactType::SCHOTTKY)) {
        pot[coord] = schottkyBuiltInPotential(coordSurf, boundPos) + (applied ? vApp : 0);
      } else if (surface.isContactType(coordSurf, boundPos, emcContactType::GATE)) {
        T barrier = device.normalizeVoltage(
            surface.getContactFurtherParameter(coordSurf, boundPos, 2));
        pot[coord] = barrier + ((applied || gateInEquilibrium) ? vApp : 0);
      }
    }
  }

  static T harmonic(T a, T b) { return 2 * a * b / (a + b); }

  //! one SOR solve; eConc = nullptr selects the equilibrium closure
  void solve(GridType &pot, const GridType *eConc) {
    SizeVec coord;
    T error;
    lastSweeps = 0;
    do {
      error = 0;
      coord.fill(0);
      for (auto &currPot : pot) {
        if (isDir[coord]) {
          pot.advanceCoord(coord);
          continue;
        }
        const T epsC = epsRel[coord];
        T sumW = 0, sumWphi = 0;
        for (SizeType d = 0; d < Dim; d++) {
          // minus neighbour (mirror the plus neighbour at an artificial boundary)
          if (coord[d] != 0) {
            T w = hFactor[d] * harmonic(epsC, epsRel.getPrevValue(coord, d));
            sumW += w;
            sumWphi += w * pot.getPrevValue(coord, d);
          } else {
            T w = hFactor[d] * harmonic(epsC, epsRel.getNextValue(coord, d));
            sumW += w;
            sumWphi += w * pot.getNextValue(coord, d);
          }
          // plus neighbour
          if (coord[d] != maxCoord[d]) {
            T w = hFactor[d] * harmonic(epsC, epsRel.getNextValue(coord, d));
            sumW += w;
            sumWphi += w * pot.getNextValue(coord, d);
          } else {
            T w = hFactor[d] * harmonic(epsC, epsRel.getPrevValue(coord, d));
            sumW += w;
            sumWphi += w * pot.getPrevValue(coord, d);
          }
        }
        T nominator = sumWphi, denominator = sumW;
        if (isSemi[coord]) {
          T n, p;
          if (eConc) { // transport: electrons from the particles, holes in equilibrium
            n = (*eConc)[coord];
            p = std::exp(-currPot);
            nominator += hProduct * (p - n + doping[coord] + currPot * (p + n));
            denominator += hProduct * (n + p);
          } else if (degeneracyBeta <= 0) { // Boltzmann closure
            n = std::exp(currPot);
            p = 1. / n;
            nominator += hProduct * (p - n + doping[coord] + currPot * (p + n));
            denominator += hProduct * (n + p);
          } else { // degenerate 2D Fermi-Dirac closure (quantum capacitance)
            T e = std::exp(currPot);
            n = std::log1p(degeneracyBeta * e) / degeneracyBeta;
            T nPrime = e / (1 + degeneracyBeta * e);
            p = std::exp(-currPot);
            nominator += hProduct * (p - n + doping[coord] + currPot * (p + nPrime));
            denominator += hProduct * (p + nPrime);
          }
        } else if (eConc) { // dielectric node: fixed particle charge only
          nominator -= hProduct * (*eConc)[coord];
        }
        T deltaPot = omega * (nominator / denominator - currPot);
        if (degeneracyBeta > 0)
          deltaPot = std::max(-potStepClamp, std::min(potStepClamp, deltaPot));
        currPot += deltaPot;
        if (std::abs(deltaPot) > error)
          error = std::abs(deltaPot);
        pot.advanceCoord(coord);
      }
      lastSweeps++;
    } while (error > accuracy && lastSweeps < maxSweeps);
  }

  //! ohmic-contact built-in potential from charge neutrality (as emcSORSolver)
  T builtInPotential(T dopingNorm) const {
    if (degeneracyBeta <= 0 || dopingNorm <= 0)
      return std::asinh(0.5 * dopingNorm);
    return std::log(std::expm1(degeneracyBeta * dopingNorm) / degeneracyBeta);
  }

  T schottkyBuiltInPotential(const SizeVecSurface &coordSurf,
                             emcBoundaryPos boundPos) const {
    if (degeneracyBeta <= 0)
      emcMessage::getInstance()
          .addError("Schottky contact requires setDegeneracyBeta(beta>0).")
          .print();
    T barrier = surface.getContactFurtherParameter(coordSurf, boundPos, 0);
    return -std::log(degeneracyBeta) - device.normalizeVoltage(barrier);
  }
};

#endif // EMC_SOR_SOLVER_DIELECTRIC_HPP
