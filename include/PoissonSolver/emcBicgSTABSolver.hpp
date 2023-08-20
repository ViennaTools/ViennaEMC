#ifndef EMC_BICGSTAB_SOLVER_HPP
#define EMC_BICGSTAB_SOLVER_HPP

#include <PoissonSolver/emcAbstractSolver.hpp>
#include <emcSurface.hpp>

#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "Eigen/IterativeLinearSolvers"
#include "Eigen/SparseCholesky"

/*! \brief Biconjugate gradient stabilized method (BISCGSTAB)
 *
 * Assumes that device given in constructor is always used! Precalculates
 * and stores some variables based on that assumption.
 *
 * @param accuracy stopping criterion, solver stops when the relative residual
 * is smaller than this value
 * @param inMaxIterations max. iterations for the solver
 */
template <class T, class DeviceType, class ParticleHandler>
class emcBicgSTABSolver : public emcAbstractSolver<T, DeviceType, ParticleHandler> {
  static const int Dim = DeviceType::Dimension;
  typedef typename DeviceType::SizeVec SizeVec;
  typedef typename DeviceType::SizeVecSurface SizeVecSurface;
  typedef typename DeviceType::ValueVec ValueVec;
  typedef emcGrid<T, Dim> GridType;

public:
  emcBicgSTABSolver(const DeviceType &inDevice, const T inAccuracy = 1e-5,
               const T inMaxIterations = 1e4)
      : accuracy(inDevice.normalizeVoltage(inAccuracy)), maxIterations(inMaxIterations),
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
  void calcEquilibriumPotential(GridType &pot, const DeviceType & /*device*/, bool resetBC = true) {
    SizeVecSurface coordSurf;
    emcBoundaryPos boundPos;
    SizeVec coord;

    // set built-in potential at ohmic contacts
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

    // initialize Matrix A + Solver
    if(initFlag) {
      GridType pot_temp = pot;
      initMatrix(pot_temp);
      initSolver();
      initFlag = false;
    }

    T p, n;
    coord.fill(0);
    size_t i = 0;

    // update the vector b
    for (auto &currPot : pot) {
      // skip if point is at ohmic contacts
      if (surface.isOhmicContact(coord)) {
        b[i] = currPot;
        pot.advanceCoord(coord);
        i++;
        continue;
      }

      n = std::exp(currPot);
      p = 1. / n;
      b[i] = hProduct * (p - n + doping[coord] + currPot * (p + n));

      for (SizeType idxDim = 0; idxDim < Dim; idxDim++) {
        if (coord[idxDim] == 0) {
          checkForGateContact(coord, idxDim, i, 0, false);
        }

        if (coord[idxDim] == maxCoord[idxDim]) {
          checkForGateContact(coord, idxDim, i, 1, false);
        }
      }
      pot.advanceCoord(coord);
      i++;
    }

    //printVector(b, "vector_b.txt");

    //BicgSTAB Solver
    x = Solver.solve(b);
    checkSolver(Solver, A, b, x);

    //Update the Grid
    for(size_t i = 0; i < pot.getSize(); i++) {
      pot[i] = x[i];
    }

    //printVector(x, "vector_x.txt");
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

    // initialize Matrix A + Solver
    if(initFlagNon) {
      GridType pot_temp = pot;
      initMatrixNon(pot_temp, eConc);
      initSolver();
      initFlagNon = false;
    }

    T p, n;
    coord.fill(0);
    size_t i = 0;

    // update the vector b
    for (auto &currPot : pot) {
      // skip if point is at ohmic contacts
      if (surface.isOhmicContact(coord)) {
        b[i] = currPot;
        pot.advanceCoord(coord);
        i++;
        continue;
      }

      p = std::exp(-currPot);
      n = eConc[coord];
      b[i] = hProduct * (p - n + doping[coord] + currPot * (p + n));

      for (SizeType idxDim = 0; idxDim < Dim; idxDim++) {
        if (coord[idxDim] == 0) {
          checkForGateContact(coord, idxDim, i, 0, true);
        }

        if (coord[idxDim] == maxCoord[idxDim]) {
          checkForGateContact(coord, idxDim, i, 1, true);
        }
      }
      pot.advanceCoord(coord);
      i++;
    }

    //BicgSTAB Solver
    x = Solver.solve(b);
    checkSolver(Solver, A, b, x);

    //Update the Grid
    for(size_t i = 0; i < pot.getSize(); i++) {
      pot[i] = x[i];
    }
  }

  /// calculates background potential
  void calcBackgroundPotential(GridType &pot, const DeviceType & /*device*/,
                               ParticleHandler &handler, bool resetBC = true) {
    SizeVecSurface coordSurf;
    emcBoundaryPos boundPos;
    SizeVec coord;

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
    
    // initialize Matrix A + Solver
    if(initFlagBackground) {
      GridType pot_temp = pot;
      initMatrixBackground(pot_temp);
      initSolver();
      initFlagBackground = false;
    }

    coord.fill(0);
    size_t i = 0;

    // update the vector b
    for (auto &currPot : pot) {
      // skip if point is at ohmic contacts
      if (surface.isOhmicContact(coord)) {
        b[i] = currPot;
        pot.advanceCoord(coord);
        i++;
        continue;
      }

      b[i] = 0;

      for (SizeType idxDim = 0; idxDim < Dim; idxDim++) {
        if (coord[idxDim] == 0) {
          handleBoundaryBackground(coord, idxDim, i, false);
        } 

        if (coord[idxDim] == maxCoord[idxDim]) {
          handleBoundaryBackground(coord, idxDim, i, true);
        }
      }
      i++;
      pot.advanceCoord(coord);
    }

    //BicgSTAB Solver
    x = Solver.solve(b);
    checkSolver(Solver, A, b, x);

    //Update the grid
    for(size_t i = 0; i < pot.getSize(); i++) {
      pot[i] = x[i];
    }
  }

private:
  T accuracy;       //!< accuracy used for stopping criterion
  T maxIterations;  //!< Max. Iterations 
  double counter1 = 0;
  double counter2 = 0;
  
  /// stored values for the system of equations
  bool initFlag = true;
  bool initFlagNon = true;
  bool initFlagBackground = true;
  Eigen::SparseMatrix<T> A;
  Eigen::VectorXd b;
  Eigen::VectorXd x;
  Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> Solver;

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
  
  ///\brief helper function that initialize the Matrix
  void initMatrixNon(GridType pot, const GridType eConc) {
    A.resize(pot.getSize(), pot.getSize());
    b.resize(pot.getSize());
    x.resize(pot.getSize());
    size_t i = 0;
    SizeVec coord;
    coord.fill(0);
    T n, p;
    for (auto currPot : pot) {
      // skip if point is at ohmic contacts
      if (surface.isOhmicContact(coord)) {
        A.coeffRef(i, i) = 1;
        pot.advanceCoord(coord);
        i++;
        continue;
      }

      p = std::exp(-currPot);
      n = eConc[coord];
      A.coeffRef(i, i) = 2 * hFactorSum + hProduct * (n + p);

      for (SizeType idxDim = 0; idxDim < Dim; idxDim++) {
        if (coord[idxDim] != 0) {
          A.coeffRef(i, pot.getPrevValue_Idx(coord, idxDim)) -= hFactor[idxDim];
        } // at Min-Boundary
        else {
          A.coeffRef(i, pot.getNextValue_Idx(coord, idxDim)) -= hFactor[idxDim];
          checkForGateContact_Matrix(coord, idxDim, i, 0);
        }

        if (coord[idxDim] != maxCoord[idxDim]) {
          A.coeffRef(i, pot.getNextValue_Idx(coord, idxDim)) -= hFactor[idxDim];
        } // at Max-Boundary
        else {
          A.coeffRef(i, pot.getPrevValue_Idx(coord, idxDim)) -= hFactor[idxDim];
          checkForGateContact_Matrix(coord, idxDim, i, 1);
        }
      }
      pot.advanceCoord(coord);
      i++;
    }
    //printMatrixData(A, pot, "matrix_data_mosfet2D.txt");
  }
  
  void initMatrix(GridType pot) {
    A.resize(pot.getSize(), pot.getSize());
    b.resize(pot.getSize());
    x.resize(pot.getSize());
    size_t i = 0;
    SizeVec coord;
    coord.fill(0);
    T n, p;
    for (auto currPot : pot) {
      // skip if point is at ohmic contacts
      if (surface.isOhmicContact(coord)) {
        A.coeffRef(i, i) = 1;
        pot.advanceCoord(coord);
        i++;
        continue;
      }

      n = std::exp(currPot);
      p = 1. / n;
      A.coeffRef(i, i) = 2 * hFactorSum + hProduct * (n + p);

      for (SizeType idxDim = 0; idxDim < Dim; idxDim++) {
        if (coord[idxDim] != 0) {
          A.coeffRef(i, pot.getPrevValue_Idx(coord, idxDim)) -= hFactor[idxDim];
        } // at Min-Boundary
        else {
          A.coeffRef(i, pot.getNextValue_Idx(coord, idxDim)) -= hFactor[idxDim];
          checkForGateContact_Matrix(coord, idxDim, i, 0);
        }

        if (coord[idxDim] != maxCoord[idxDim]) {
          A.coeffRef(i, pot.getNextValue_Idx(coord, idxDim)) -= hFactor[idxDim];
        } // at Max-Boundary
        else {
          A.coeffRef(i, pot.getPrevValue_Idx(coord, idxDim)) -= hFactor[idxDim];
          checkForGateContact_Matrix(coord, idxDim, i, 1);
        }
      }
      pot.advanceCoord(coord);
      i++;
    }
    //printMatrixData(A, pot, "matrix_data_mosfet2D.txt");
  }

  ///\brief helper function that initialize the Matrix for background potential
  void initMatrixBackground(GridType pot) {
    A.resize(pot.getSize(), pot.getSize());
    b.resize(pot.getSize());
    x.resize(pot.getSize());
    size_t i = 0;
    SizeVec coord;
    coord.fill(0);

    for (auto &currPot : pot) {
      // skip if point is at ohmic contacts
      if (surface.isOhmicContact(coord)) {
        A.coeffRef(i, i) = 1;
        pot.advanceCoord(coord);
        i++;
        continue;
      }

      A.coeffRef(i, i) = 2 * hSquaredFactorSum;

      for (SizeType idxDim = 0; idxDim < Dim; idxDim++) {
        if (coord[idxDim] != 0) {
          A.coeffRef(i, pot.getPrevValue_Idx(coord, idxDim)) -= hSquaredFactor[idxDim];
        }
        else {
          A.coeffRef(i, pot.getNextValue_Idx(coord, idxDim)) -= hSquaredFactor[idxDim];
          handleBoundaryBackground_Matrix(coord, idxDim, i, false);
        }

        if (coord[idxDim] != maxCoord[idxDim]) {
          A.coeffRef(i, pot.getNextValue_Idx(coord, idxDim)) -= hSquaredFactor[idxDim];
        }
        else {
          A.coeffRef(i, pot.getPrevValue_Idx(coord, idxDim)) -= hSquaredFactor[idxDim];
          handleBoundaryBackground_Matrix(coord, idxDim, i, true);
        }
      }
      i++;
      pot.advanceCoord(coord);
    }
  }

  ///\brief helper function that initialize the Solver
  void initSolver() {
    Solver.setMaxIterations(maxIterations);
    Solver.setTolerance(accuracy);
    Solver.compute(A);
    checkDecomposition(Solver);
  }

  //----------------------------------------------------------------------------------
  /*! \brief Helper function that adapts parameter when current coordinate
   * is at boundary (for background potential calculation).
   * @param isMaxBound boolean that selects max. or min. boundary in that
   * direction
   */
  void handleBoundaryBackground(const SizeVec &coord, SizeType idxDim,
                                size_t i, bool isMaxBound) {
    emcBoundaryPos boundPos =
        static_cast<emcBoundaryPos>(2 * idxDim + isMaxBound);
    SizeVecSurface coordSurf = surface.getCoordBoundary(coord, boundPos);
    if (surface.isArtificialBoundary(coordSurf, boundPos)) {
      b[i] -= neumannBC[coord] * hSquaredFactor[idxDim];
    } else { // gate - contact
      auto idxCont = surface.getContactIdx(coordSurf, boundPos);
      b[i] += hSquaredFactor[idxDim] *
                   (gateFactor[idxCont] * h[idxDim] *
                        (gateBarrierHeight[idxCont] + gateVoltage[idxCont] -
                         potPart[coord]) +
                    2 * neumannBC[coord]);
    }
  }
  
  void handleBoundaryBackground_Matrix(const SizeVec &coord, SizeType idxDim,
                                size_t i, bool isMaxBound) {
    emcBoundaryPos boundPos = static_cast<emcBoundaryPos>(2 * idxDim + isMaxBound);
    SizeVecSurface coordSurf = surface.getCoordBoundary(coord, boundPos);
    if (surface.isArtificialBoundary(coordSurf, boundPos) == false) {
      auto idxCont = surface.getContactIdx(coordSurf, boundPos);
      A.coeffRef(i, i) += gateFactor[idxCont] * h[idxDim] * hSquaredFactor[idxDim];
    }
  }
  
  //----------------------------------------------------------------------------------
  /*! \brief Helper function for Gate Contact, adapts parameter if gate contact
   * is present.
   * @param isMaxBound boolean that selects max. or min. boundary in that
   * direction
   * @param useVApp boolean that tells if applied voltage should be added or not
   */
  void checkForGateContact(const SizeVec &coord, SizeType idxDim,
                           size_t i, bool isMaxBound, bool useVApp) {
    emcBoundaryPos boundPos = static_cast<emcBoundaryPos>(2 * idxDim + isMaxBound);
    SizeVecSurface coordSurf = surface.getCoordBoundary(coord, boundPos);
    if (surface.isContactType(coordSurf, boundPos, emcContactType::GATE)) {
      auto idxCont = surface.getContactIdx(coordSurf, boundPos);
      T Vg = gateBarrierHeight[idxCont];
      if (useVApp) {
        Vg += gateVoltage[idxCont];
      }
      b[i] += gateFactor[idxCont] * Vg * hFactor[idxDim] * h[idxDim];
      //std::cout << i << std::endl;
    }
  }

  void checkForGateContact_Matrix(const SizeVec &coord, SizeType idxDim,
                           size_t i, bool isMaxBound) {
    emcBoundaryPos boundPos = static_cast<emcBoundaryPos>(2 * idxDim + isMaxBound);
    SizeVecSurface coordSurf = surface.getCoordBoundary(coord, boundPos);
    if (surface.isContactType(coordSurf, boundPos, emcContactType::GATE)) {
      auto idxCont = surface.getContactIdx(coordSurf, boundPos);
      A.coeffRef(i, i) += gateFactor[idxCont] * hFactor[idxDim] * h[idxDim];
    }
  }
  
  //----------------------------------------------------------------------------------
  ///\brief Helper functions to check for potential errors of the solver
  void checkSolver(const Eigen::BiCGSTAB<Eigen::SparseMatrix<T>> &Solver, 
                   const Eigen::SparseMatrix<T> &A, Eigen::VectorXd &b, Eigen::VectorXd &x) {

    //Calculate relative residual
    /*Eigen::VectorXd residual = b - A * x;
    double residual_norm = residual.norm()/b.norm();

    if(residual_norm > accuracy) {
      std::cerr << "Accuracy not reached!" << std::endl;
      std::cout << "Number of iterations: " << Solver.iterations() << std::endl;
      std::cout << "Accuracy: " << accuracy << std::endl;
      std::cout << "Residual: " << residual_norm << std::endl;
      abort();
    }*/

    if(Solver.info() == Eigen::NumericalIssue) {
      std::cerr << "Invalid Calculation!" << std::endl;
      abort();
    }

    if(Solver.info() == Eigen::InvalidInput) {
      std::cerr << "Invalid Input!" << std::endl;
      abort(); 
    }

    if(Solver.info() == Eigen::NoConvergence) {
      std::cerr << "Does not Converge!" << std::endl;
      abort(); 
    }

    

    if(Solver.info() != Eigen::Success) {
      std::cerr << "Solving failed!" << std::endl;
      abort(); 
    }
  }

  void checkDecomposition(const Eigen::BiCGSTAB<Eigen::SparseMatrix<T>> &Solver) {
    if(Solver.info() != Eigen::Success) {
      std::cerr << "Decomposition was not successful" << std::endl;
      abort();
    }
  }
  
  //----------------------------------------------------------------------------------
  ///\brief Helper functions to print all non zero elements of Matrix + Coord + Type (for debugging)
  void printMatrixData(const Eigen::SparseMatrix<double>& matrix, GridType pot, const std::string& filename) {
    SizeVec coord;
    coord.fill(0);

    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        std::cerr << "Failed to open the file." << std::endl;
        return;
    }

    outfile << "Points in x: " << pot.getSize(0) << std::endl;
    outfile << "Points in y: " << pot.getSize(1) << std::endl;
    outfile << "Number of all grid Points: " << pot.getSize() << std::endl;

    // For each row
    for (int i = 0; i < matrix.rows(); ++i) {
        // Print coordinates
        printCoord(coord, pot, outfile);

         auto [type, gateValue] = pointType(coord);
        outfile << ", Type: " << type;

        if (type == "Gate_Min" || type == "Gate_Max") {
          outfile << ", Gate Value: " << gateValue;
        }
        
        outfile << ", Row: " << i;

        // Check each column in the current row
        for (int j = 0; j < matrix.cols(); ++j) {
            double value = matrix.coeff(i, j);
            if (value != 0) {  // if it's not zero
                outfile << " [Column: " << j << ", Value: " << value << "]";
            }
        }
        
        outfile << std::endl; // Move to the next line after processing all entries of a row

        // Advance coordinates
        pot.advanceCoord(coord);
    }
    outfile.close();
  }
  
  void printCoord(const SizeVec& coord, GridType pot, std::ostream& os) {
    os << "(";
    for (SizeType j = 0; j < Dim; j++) {
      os << coord[j];
      if (j < Dim - 1) os << ", ";
    }
    os << ")";
  }

  std::pair<std::string, double> pointType(const SizeVec &coord) {
    if (surface.isOhmicContact(coord)) {
      return {"Ohmic", 0.0};
    }

    for(SizeType idxDim = 0; idxDim < Dim; idxDim++) {
      // Check for Min-Boundary
      if (coord[idxDim] == 0) {
        emcBoundaryPos boundPos = static_cast<emcBoundaryPos>(2 * idxDim);
        SizeVecSurface coordSurf = surface.getCoordBoundary(coord, boundPos);
        if (surface.isContactType(coordSurf, boundPos, emcContactType::GATE)) {
            auto idxCont = surface.getContactIdx(coordSurf, boundPos);
            double value = gateFactor[idxCont] * hFactor[idxDim] * h[idxDim];
            return {"Gate_Min", value};
        }
      }
      // Check for Max-Boundary
      if (coord[idxDim] == maxCoord[idxDim]) {
        emcBoundaryPos boundPos = static_cast<emcBoundaryPos>(2 * idxDim + 1);
        SizeVecSurface coordSurf = surface.getCoordBoundary(coord, boundPos);
        if (surface.isContactType(coordSurf, boundPos, emcContactType::GATE)) {
            auto idxCont = surface.getContactIdx(coordSurf, boundPos);
            double value = gateFactor[idxCont] * hFactor[idxDim] * h[idxDim];
            return {"Gate_Max", value};
        }
      }
    }
    return {"None", 0.0};
  }

  ///\brief Helper functions to print all elements of Vector
  void printVector(const Eigen::VectorXd& vector, const std::string& filename = "vector_data.txt") {
    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        std::cerr << "Failed to open the file." << std::endl;
        return;
    }
    
    // Iterate over each element of the vector and write its index and value to the file
    for (int i = 0; i < vector.size(); ++i) {
        outfile << "Index: " << i << ", Value: " << vector[i] << std::endl;
    }
    outfile.close();
  }

  ///\brief Helper functions to print doping profile
  void printDoping(const std::string& filename = "doping_data.txt") {
    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        std::cerr << "Failed to open the file." << std::endl;
        return;
    }
    
    SizeVec coord;
    coord.fill(0);

    // Assuming you iterate through the doping grid similar to how you do with the GridType pot
    for (auto currentDoping : doping) {
      outfile << "Coordinate: (";
      for (SizeType dim = 0; dim < Dim - 1; ++dim) {
          outfile << coord[dim] << ", ";
      }
      outfile << coord[Dim-1] << ") Value: " << doping[coord] << std::endl;
      doping.advanceCoord(coord);  // Assuming you have this function or similar to go to the next coordinate
    }

    outfile.close();
  }
};

#endif // EMC_BICGSTAB_SOLVER_HPP