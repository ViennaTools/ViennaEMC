#ifndef EMC_SOR_SOLVER_HPP
#define EMC_SOR_SOLVER_HPP

#include <PoissonSolver/emcAbstractSolver.hpp>
#include <emcSurface.hpp>

/*! \brief Successive Over-Relaxation Poisson Solver (SOR-Solver)
 *
 * Assumes that device given in constructor is always used! Precalculates
 * and stores some variables based on that assumption.
 *
 * @param accuracy stopping criterion, solver stops when max. change in solution
 * is smaller than this value
 * @param omega relaxation factor omega (appropriate values in (0,2))
 */
template <class T, class DeviceType, class ParticleHandler>
class emcSORSolver : public emcAbstractSolver<T, DeviceType, ParticleHandler> {
  static const int Dim = DeviceType::Dimension;
  typedef typename DeviceType::SizeVec SizeVec;
  typedef typename DeviceType::SizeVecSurface SizeVecSurface;
  typedef typename DeviceType::ValueVec ValueVec;
  typedef emcGrid<T, Dim> GridType;

public:
  emcSORSolver(const DeviceType &inDevice, const T inAccuracy = 1e-5,
               const T inOmega = 1.5)
      : accuracy(inDevice.normalizeVoltage(inAccuracy)), omega(inOmega),
        device(inDevice), surface(inDevice.getSurface()),
        doping(inDevice.getDopingProfile().getDoping(true)),
        neumannBC(inDevice.getGridExtent(), 0),
        potPart(inDevice.getGridExtent(), 0), h(inDevice.getSpacing(true)),
        gateFactor(surface.getNrContacts()),
        gateVoltage(surface.getNrContacts()),
        gateBarrierHeight(surface.getNrContacts()) {
    initDeviceCharacteristics();
  }

  /// calculates equilibrium potential
  void calcEquilibriumPotential(GridType &pot, const DeviceType & /*device*/,
                                bool resetBC = true) {
    SizeVecSurface coordSurf;
    emcBoundaryPos boundPos;
    SizeVec coord;

    if (resetBC) {
      // set built-in potential at ohmic contacts
      for (surface.initCoord(boundPos, coordSurf);
           !surface.isEndCoord(boundPos, coordSurf);
           surface.advanceCoord(boundPos, coordSurf)) {
        if (surface.isContactType(coordSurf, boundPos, emcContactType::OHMIC)) {
          coord = surface.getCoordDevice(coordSurf, boundPos);
          pot[coord] = std::asinh(0.5 * doping[coord]);
        }
      }
    }

    // update potential in iterative way until error is small enough
    T error, deltaPot, p, n, nominator, denominator;
    do {
      error = 0;
      coord.fill(0);
      for (auto &currPot : pot) {
        // skip if point is at ohmic contacts
        if (surface.isOhmicContact(coord)) {
          pot.advanceCoord(coord);
          continue;
        }

        n = std::exp(currPot);
        p = 1. / n;
        nominator = hProduct * (p - n + doping[coord] + currPot * (p + n));
        denominator = 2 * hFactorSum + hProduct * (n + p);
        for (SizeType idxDim = 0; idxDim < Dim; idxDim++) {
          if (coord[idxDim] != 0)
            nominator += pot.getPrevValue(coord, idxDim) * hFactor[idxDim];
          else { // at Min-Boundary
            nominator += pot.getNextValue(coord, idxDim) * hFactor[idxDim];
            checkForGateContact(coord, idxDim, nominator, denominator, 0,
                                false);
          }
          if (coord[idxDim] != maxCoord[idxDim])
            nominator += pot.getNextValue(coord, idxDim) * hFactor[idxDim];
          else { // at Max-Boundary
            nominator += pot.getPrevValue(coord, idxDim) * hFactor[idxDim];
            checkForGateContact(coord, idxDim, nominator, denominator, 1,
                                false);
          }
        }

        // update potential + store max. change of potential
        deltaPot = omega * (nominator / denominator - currPot);
        currPot += deltaPot;
        if (std::abs(deltaPot) > error)
          error = std::abs(deltaPot);
        pot.advanceCoord(coord);
      }
    } while (error > accuracy);
  }

  /// calculates non-equilibrium potential
  void calcNonEquilibriumPotential(GridType &pot, const DeviceType & /*device*/,
                                   const GridType &eConc, bool resetBC = true) {
    SizeVecSurface coordSurf;
    emcBoundaryPos boundPos;
    SizeVec coord;

    // set applied + built-in potential at ohmic contacts
    if (resetBC) {
      for (surface.initCoord(boundPos, coordSurf);
           !surface.isEndCoord(boundPos, coordSurf);
           surface.advanceCoord(boundPos, coordSurf)) {
        if (surface.isContactType(coordSurf, boundPos, emcContactType::OHMIC)) {
          coord = surface.getCoordDevice(coordSurf, boundPos);
          pot[coord] = surface.getContactVoltage(coordSurf, boundPos, true) +
                       std::asinh(0.5 * doping[coord]);
        }
      }
    }

    // update potential in an iterative way until error is small enough
    T error, deltaPot, p, n, nominator, denominator;
    do {
      error = 0;
      coord.fill(0);
      for (auto &currPot : pot) {
        // skip if point is at ohmic contacts
        if (surface.isOhmicContact(coord)) {
          pot.advanceCoord(coord);
          continue;
        }

        p = std::exp(-currPot);
        n = eConc[coord];
        nominator = hProduct * (p - n + doping[coord] + currPot * (p + n));
        denominator = 2 * hFactorSum + hProduct * (n + p);
        for (SizeType idxDim = 0; idxDim < Dim; idxDim++) {
          if (coord[idxDim] != 0)
            nominator += pot.getPrevValue(coord, idxDim) * hFactor[idxDim];
          else {
            nominator += pot.getNextValue(coord, idxDim) * hFactor[idxDim];
            checkForGateContact(coord, idxDim, nominator, denominator, 0, true);
          }
          if (coord[idxDim] != maxCoord[idxDim])
            nominator += pot.getNextValue(coord, idxDim) * hFactor[idxDim];
          else {
            nominator += pot.getPrevValue(coord, idxDim) * hFactor[idxDim];
            checkForGateContact(coord, idxDim, nominator, denominator, 1, true);
          }
        }

        // update potential
        deltaPot = omega * (nominator / denominator - currPot);
        currPot += deltaPot;
        if (std::abs(deltaPot) > error)
          error = std::abs(deltaPot);
        pot.advanceCoord(coord);
      }
    } while (error > accuracy);
  }

  /// calculates background potential
  void calcBackgroundPotential(GridType &pot, const DeviceType & /*device*/,
                               ParticleHandler &handler, bool resetBC = true) {
    SizeVecSurface coordSurf;
    emcBoundaryPos boundPos;

    if (resetBC) {
      for (surface.initCoord(boundPos, coordSurf);
           !surface.isEndCoord(boundPos, coordSurf);
           surface.advanceCoord(boundPos, coordSurf)) {
        if (surface.isContactType(coordSurf, boundPos, emcContactType::OHMIC)) {
          // potential = applied potential - potential from particles
          auto coordDev = surface.getCoordDevice(coordSurf, boundPos);
          T partPot = device.normalizeVoltage(
              handler.getParticlePotential(device.coordToPos(coordDev)));
          pot[coordDev] =
              surface.getContactVoltage(coordSurf, boundPos, true) - partPot;
        } else {
          // neumann BC = approx. neg. derivate of particle potential
          SizeVec coordDev = surface.getCoordDevice(coordSurf, boundPos);
          ValueVec posContact = device.coordToPos(coordDev);
          ValueVec posNbrIn(posContact), posNbrOut(posContact);
          SizeType idxDim = static_cast<SizeType>(boundPos) / 2;
          SizeType dir = std::pow(-1, (static_cast<SizeType>(boundPos) % 2));
          posNbrIn[idxDim] += dir * device.getSpacing()[idxDim];
          posNbrOut[idxDim] -= dir * device.getSpacing()[idxDim];
          // approx. with central difference scheme
          neumannBC[coordDev] =
              0.5 *
              device.normalizeVoltage((handler.getParticlePotential(posNbrOut) -
                                       handler.getParticlePotential(posNbrIn)));
          if (surface.isContactType(coordSurf, boundPos,
                                    emcContactType::GATE)) {
            potPart[coordDev] = device.normalizeVoltage(
                handler.getParticlePotential(device.coordToPos(coordDev)));
          }
        }
      }
    }

    // update potential in iterative way
    SizeVec coord;
    T error;
    do {
      error = 0;
      coord.fill(0);
      for (auto &currPot : pot) {
        if (surface.isOhmicContact(coord)) {
          pot.advanceCoord(coord);
          continue;
        }

        T nominator = 0;
        T denominator = 2 * hSquaredFactorSum;
        for (SizeType idxDim = 0; idxDim < Dim; idxDim++) {
          if (coord[idxDim] != 0) {
            nominator +=
                pot.getPrevValue(coord, idxDim) * hSquaredFactor[idxDim];
          } else {
            nominator +=
                pot.getNextValue(coord, idxDim) * hSquaredFactor[idxDim];
            handleBoundaryBackground(coord, idxDim, nominator, denominator,
                                     false);
          }
          if (coord[idxDim] != maxCoord[idxDim]) {
            nominator +=
                pot.getNextValue(coord, idxDim) * hSquaredFactor[idxDim];
          } else {
            nominator +=
                pot.getPrevValue(coord, idxDim) * hSquaredFactor[idxDim];
            handleBoundaryBackground(coord, idxDim, nominator, denominator,
                                     true);
          }
        }

        T deltaPot = omega * (nominator / denominator - currPot);
        currPot += deltaPot;
        if (std::abs(deltaPot) > error)
          error = std::abs(deltaPot);
        pot.advanceCoord(coord);
      }
    } while (error > accuracy);
  }

private:
  T accuracy; //!< accuracy used for stopping criterion
  T omega;    //!< parameter of SOR-Solver

  const DeviceType &device; //!< used device
  const typename DeviceType::SurfaceType &surface;

  /// stored + precalculated characteristics of device
  ValueVec hFactor, hSquaredFactor, h;
  T hFactorSum, hProduct, hSquaredFactorSum;
  SizeVec maxCoord;
  std::vector<T> gateVoltage, gateFactor, gateBarrierHeight;
  GridType doping;
  GridType neumannBC, potPart; // TODO change container!

  ///\brief helper function that stores + calculates device characteristics
  void initDeviceCharacteristics() {
    // store device spacing + extent
    if (Dim == 2) {
      hFactor[0] = h[1] / h[0];
      hFactor[1] = h[0] / h[1];

      hSquaredFactor[0] = h[1] * h[1];
      hSquaredFactor[1] = h[0] * h[0];
    } else {
      hFactor[0] = h[1] * h[2] / h[0];
      hFactor[1] = h[0] * h[2] / h[1];
      hFactor[2] = h[0] * h[1] / h[2];

      hSquaredFactor[0] = h[1] * h[1] * h[2] * h[2];
      hSquaredFactor[1] = h[0] * h[0] * h[2] * h[2];
      hSquaredFactor[2] = h[0] * h[0] * h[1] * h[1];
    }
    hFactorSum = std::accumulate(hFactor.begin(), hFactor.end(), 0.);
    hSquaredFactorSum =
        std::accumulate(hSquaredFactor.begin(), hSquaredFactor.end(), 0.);
    hProduct = std::accumulate(h.begin(), h.end(), 1., std::multiplies<T>());
    maxCoord = device.getGridExtent();
    std::for_each(maxCoord.begin(), maxCoord.end(),
                  [](SizeType &val) { val--; });

    // store gate - contact characteristics
    T epsRSc = device.getMaterial().getEpsR();
    for (SizeType idxCont = 0; idxCont < gateFactor.size(); idxCont++) {
      gateVoltage[idxCont] = surface.getContactVoltage(idxCont, true);
      if (surface.getContactType(idxCont) == emcContactType::GATE) {
        T gamma = surface.getContactFurtherParameter(idxCont, 0) / epsRSc;
        T thicknessOx = device.normalizeLength(
            surface.getContactFurtherParameter(idxCont, 1));
        gateFactor[idxCont] = 2 * gamma / thicknessOx;
        gateBarrierHeight[idxCont] = device.normalizeVoltage(
            surface.getContactFurtherParameter(idxCont, 2));
      }
    }
  }

  /*! \brief Helper function that adapts parameter when current coordinate
   * is at boundary (for background potential calculation).
   * @param isMaxBound boolean that selects max. or min. boundary in that
   * direction
   */
  void handleBoundaryBackground(const SizeVec &coord, SizeType idxDim,
                                T &nominator, T &denominator, bool isMaxBound) {
    emcBoundaryPos boundPos =
        static_cast<emcBoundaryPos>(2 * idxDim + isMaxBound);
    SizeVecSurface coordSurf = surface.getCoordBoundary(coord, boundPos);
    if (surface.isArtificialBoundary(coordSurf, boundPos)) {
      nominator -= neumannBC[coord] * hSquaredFactor[idxDim];
    } else { // gate - contact
      auto idxCont = surface.getContactIdx(coordSurf, boundPos);
      nominator += hSquaredFactor[idxDim] *
                   (gateFactor[idxCont] * h[idxDim] *
                        (gateBarrierHeight[idxCont] + gateVoltage[idxCont] -
                         potPart[coord]) +
                    2 * neumannBC[coord]);
      denominator += gateFactor[idxCont] * h[idxDim] * hSquaredFactor[idxDim];
    }
  }

  /*! \brief Helper function for Gate Contact, adapts parameter if gate contact
   * is present.
   * @param isMaxBound boolean that selects max. or min. boundary in that
   * direction
   * @param useVApp boolean that tells if applied voltage should be added or not
   */
  void checkForGateContact(const SizeVec &coord, SizeType idxDim, T &nominator,
                           T &denominator, bool isMaxBound, bool useVApp) {
    emcBoundaryPos boundPos =
        static_cast<emcBoundaryPos>(2 * idxDim + isMaxBound);
    SizeVecSurface coordSurf = surface.getCoordBoundary(coord, boundPos);
    if (surface.isContactType(coordSurf, boundPos, emcContactType::GATE)) {
      auto idxCont = surface.getContactIdx(coordSurf, boundPos);
      T Vg = gateBarrierHeight[idxCont];
      if (useVApp)
        Vg += gateVoltage[idxCont];
      nominator += gateFactor[idxCont] * Vg * hFactor[idxDim] * h[idxDim];
      denominator += gateFactor[idxCont] * hFactor[idxDim] * h[idxDim];
    }
  }
};

#endif // EMC_SOR_SOLVER_HPP