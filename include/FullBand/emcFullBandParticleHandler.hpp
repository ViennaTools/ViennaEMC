#ifndef EMC_FULL_BAND_PARTICLE_HANDLER_HPP
#define EMC_FULL_BAND_PARTICLE_HANDLER_HPP

/*! \brief Phase 4: full-band carriers inside ViennaEMC's DEVICE machinery.
 *
 * WHAT IT KEEPS. Everything in the device layer that never touches a band
 * structure: emcDevice (geometry, doping regions, surface/contacts), the
 * Poisson solver, the particle-mesh scheme (force interpolation, charge
 * assignment), emcSimulation's step loop, contact-reservoir bookkeeping
 * (expected carriers per contact cell, injection/removal counters), results
 * and output. The analytic emcParticleType is kept ONLY for bookkeeping -
 * charge, name, expected carriers at a contact, initial carriers per cell -
 * none of which involves a valley; it must carry NO valleys and NO
 * mechanisms, and none of its valley/scatter methods is ever called.
 *
 * WHAT IT REPLACES. The carrier dynamics. In ViennaEMC the free-flight
 * position step is inline effective-mass math (emcParticleDrift.hpp:
 * dPos = hbar k dt / m_c), every final-state routine rebuilds k from an
 * analytic dispersion, and injection is an analytic Maxwellian. Here each
 * particle carries a Cartesian k [1/m], a BAND index and a tetrahedron hint,
 * and the per-particle loop is the one validated in the bulk driver
 * (examples/fullBandBulk): slab-wise self-scattering bound Gamma0(E), the
 * slab-crossing flight cap, k-resolved rates, |g|^2 final states with band
 * changes, thermal injection by DOS-weighted energy and thermal band split.
 * The only change from the bulk loop is that the field is the LOCAL force
 * from the particle-mesh scheme instead of a uniform F.
 *
 * DOPING. The impurity ladder level is selected from the device's doping
 * (setDoping on every band's scattering object). Milestone 1 assumes ONE
 * doping region; a multi-region device needs the level per region, which
 * the engine's isActive() gate supports but this handler does not yet wire.
 *
 * BOUNDARIES. An ohmic-contact boundary removes the particle (the reservoir
 * refills it); any other boundary reflects specularly: the k component
 * normal to the boundary flips and the position mirrors. That conserves
 * energy exactly only when the boundary normal is a crystal axis with E(k)
 * mirror-symmetric about it - true for cubic Si with device axes along
 * <100>, which is what the first milestone uses.
 *
 * HOLES. A hole package stores E = VBM - E_electron with negated velocities:
 * the carrier is a formal ELECTRON on the inverted band, and its k-space
 * trajectory - integrated with force -q E exactly as for electrons, which is
 * what the bulk driver validated - is the MIRROR of the physical hole's. So
 * the k-update, the energy, the scattering and the injection are carrier-
 * agnostic; only the real-space velocity flips: v_physical = -v_stored. The
 * framework's emcHole type supplies +q for Poisson and the current tally.
 * The carrier is read from the package's own `carrier` attribute.
 *
 * CONFIGURATION. emcSimulation constructs the handler itself with a fixed
 * signature, so the package path and run knobs are static members the
 * driver sets BEFORE constructing the simulation.
 */

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <memory>
#include <string>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

#include <FullBand/emcFullBandScattering.hpp>
#include <FullBand/emcNumericBandStructure.hpp>
#include <ParticleHandler/emcAbstractParticleHandler.hpp>
#include <emcConstants.hpp>
#include <emcParticleInitialization.hpp>

template <class T, class DeviceType, class PMScheme>
class emcFullBandParticleHandler
    : public emcAbstractParticleHandler<T, DeviceType, PMScheme,
                                        DeviceType::Dimension> {
  static constexpr SizeType Dim = DeviceType::Dimension;
  typedef emcAbstractParticleHandler<T, DeviceType, PMScheme, Dim> Base;
  typedef typename Base::ParticleType ParticleType;
  typedef typename Base::SizeVec SizeVec;
  typedef typename Base::ValueVec ValueVec;
  typedef typename Base::MapIdxToParticleTypes MapIdxToParticleTypes;
  typedef typename Base::NettoParticleCounter NettoParticleCounter;
  typedef emcNumericBandStructure<T> NBS;
  typedef emcFullBandScattering<T> FBS;

public:
  /// --- configuration, set by the driver BEFORE the simulation is built ---
  static inline std::string packagePath;
  static inline T temperature = T(300);
  static inline T dopingCm3 = T(0);      ///< 0: take |doping| from the device
  static inline T binWidth = T(2e-3);    ///< energy-bin width [eV]
  static inline int mode = 0;            ///< 0 = best available (phaseB + g2)
  static inline bool useSlabs = true;    ///< slab-wise Gamma0(E) (as bulk)
  static inline T slabWidth = T(0.010);  ///< [eV]
  static inline T dtCapScale = T(0.5);   ///< slab-crossing cap, as bulk
  static inline SizeType avgFromStep = 0;///< accumulate averages from here
  static inline int nrBands = 0;         ///< 0: every band in the package

  struct FBParticle {
    std::array<T, 3> k{0, 0, 0}; ///< Cartesian wave vector [1/m]
    int band = 0;
    std::int64_t hint = -1;      ///< tetrahedron hint for locateTet
    T energy = 0;                ///< absolute band energy [eV]
    int level = -1;              ///< impurity ladder level of the current doping region
  };

private:
  std::vector<std::vector<FBParticle>> particles;
  std::vector<std::vector<std::array<T, Dim>>> positions;
  ValueVec deviceSpacing;
  SizeVec deviceExtent;
  ValueVec deviceMaxPos;

  std::unique_ptr<NBS> bs;
  std::vector<std::unique_ptr<FBS>> scatBand;
  T cbm = 0;
  std::vector<T> gamma0B, slabWB, vMaxB, zB;
  std::vector<int> levelOfRegion;   ///< doping region index -> ladder level
  T carrierSign = 1;                ///< +1 electron package, -1 hole package
  T zTot = 0;
  int nb = 1;

  // time-weighted averages over the run (from avgFromStep on)
  SizeType stepCounter = 0;
  // exact double accumulation: per-thread partials reduced once per step
  double wV[3]{0, 0, 0}, wE = 0, wSum = 0, wBand[8]{};
  // x-resolved profiles (one bin per grid column): time-weighted vx, E, count
  std::vector<double> pxV, pxE, pxN;
  std::atomic<long long> nReal{0}, nSelf{0}, nRemoved{0}, nReflected{0};

public:
  emcFullBandParticleHandler() = delete;

  emcFullBandParticleHandler(const DeviceType &inDevice, PMScheme &inPMScheme,
                             MapIdxToParticleTypes &inTypes,
                             SizeType inNrCarriersPerParticle, SizeType inSeed)
      : Base(inDevice, inPMScheme, inNrCarriersPerParticle, inTypes, inSeed),
        deviceSpacing(inDevice.getSpacing()),
        deviceExtent(inDevice.getGridExtent()),
        deviceMaxPos(inDevice.getMaxPos()) {
    particles.resize(Base::getNrParticleTypes());
    positions.resize(Base::getNrParticleTypes());
    loadPackage();
  }

  bool calcsPartPartInteraction() const { return false; }

  SizeType getNrParticles(SizeType idxType) const {
    return positions[idxType].size();
  }

  void printNrParticles() const {
    for (const auto &type : Base::idxTypeToPartType)
      std::cout << "\t" << positions[type.first].size() << " "
                << type.second->getName() << " (full-band)\n";
  }

  /// ---------------------------------------------------------------- step
  NettoParticleCounter
  driftScatterParticles(T tStep, std::vector<emcGrid<T, Dim>> &eField) {
    auto &rngs = Base::rngs;
    auto &surface = Base::device.getSurface();
    auto nettoNrPart = Base::initNettoParticleCounter();
    const bool accumulate = stepCounter >= avgFromStep;
    stepCounter++;
    for (const auto &partType : Base::idxTypeToPartType) {
      const SizeType idxType = partType.first;
      auto &type = partType.second;
      if (!type->isMoved())
        continue;
      // k-space dynamics ALWAYS with -q: the package's carrier is a formal
      // electron (see HOLES above); the type's charge sign is for Poisson
      const T charge = -constants::q;
      std::vector<SizeType> idxToRemove;
#pragma omp parallel
      {
        const SizeType idxThread = this->getCurrentThreadIdx();
        auto &rng = rngs[idxThread];
        std::uniform_real_distribution<T> U01(0, 1);
        double tV[3]{0, 0, 0}, tE = 0, tS = 0, tB[8]{};
        std::vector<double> pV(deviceExtent[0], 0.0), pE(deviceExtent[0], 0.0), pN(deviceExtent[0], 0.0);
#pragma omp for
        for (SizeType ip = 0; ip < getNrParticles(idxType); ip++) {
          auto &p = particles[idxType][ip];
          auto &pos = positions[idxType][ip];
          bool removed = false;
          T t = 0;
          while (t < tStep && !removed) {
            const auto force = Base::pmScheme.interpolateForce(
                eField, pos, deviceSpacing, charge);          // [N]
            const T fMag = std::sqrt(force[0] * force[0] + force[1] * force[1] +
                                     force[2] * force[2]) / constants::q; // [V/m]
            FBS &sc = *scatBand[p.band < nb ? p.band : 0];
            const T eCur = p.energy;
            const T g0 = useSlabs ? sc.gamma0At(eCur) : gamma0B[p.band];
            const T tau = -std::log(1 - U01(rng)) / g0;
            const T dtCap =
                (useSlabs && slabWB[p.band] > 0 && fMag * vMaxB[p.band] > 0)
                    ? dtCapScale * slabWB[p.band] / (fMag * vMaxB[p.band])
                    : std::numeric_limits<T>::infinity();
            const T dt = std::min(std::min(tau, dtCap), tStep - t);
            const bool capped = (dtCap < tau) && (dtCap < tStep - t);
            const auto vCur = bs->getVelocity(p.k, p.band, p.hint);
            for (int i = 0; i < 3; i++)
              p.k[i] += force[i] * dt / constants::hbar;
            const auto vEnd = bs->getVelocity(p.k, p.band, p.hint);
            T vm[3];
            for (int c = 0; c < 3; c++) vm[c] = carrierSign * T(0.5) * (vCur[c] + vEnd[c]);
            for (SizeType c = 0; c < Dim; c++) pos[c] += vm[c] * dt;
            if (accumulate) {
              for (int c = 0; c < 3; c++) tV[c] += vm[c] * dt;
              tE += (eCur - cbm) * dt;
              tS += dt;
              if (p.band < 8) tB[p.band] += dt;
              const SizeType ix = std::min<SizeType>(deviceExtent[0] - 1,
                  static_cast<SizeType>(std::max(T(0), pos[0]) / deviceSpacing[0] + T(0.5)));
              pV[ix] += vm[0] * dt; pE[ix] += (eCur - cbm) * dt; pN[ix] += dt;
            }
            t += dt;
            // ---- boundaries: ohmic contact removes, anything else reflects
            removed = handleBoundary(p, pos, surface);
            if (removed) break;
            p.energy = bs->getEnergy(p.k, p.band, p.hint);
            if (capped || t >= tStep)
              continue; // flight interrupted: redraw (memoryless) next time
            // ---- end of a free flight: real or self-scattering
            // doping level of the region the particle is in NOW
            {
              const auto coord = Base::device.posToCoord(pos);
              const int reg = Base::device.getDopingProfile().getDopingRegionIdx(coord);
              p.level = (reg >= 0 && reg < (int)levelOfRegion.size()) ? levelOfRegion[reg] : -1;
            }
            const T rateNow = sc.getTotalRate(p.k, p.energy, p.hint, p.level);
            if (rateNow > g0) sc.noteG0Violation();
            if (U01(rng) * g0 < rateNow) {
              nReal.fetch_add(1, std::memory_order_relaxed);
              const std::size_t m = sc.selectMechanism(p.k, p.energy, p.hint, rng, p.level);
              const T Ef = p.energy + sc.getDeltaE(m);
              typename NBS::Vec3 kNew;
              std::int64_t hNew = p.hint;
              int bNew = p.band;
              if (sc.sampleFinalStateG2(m, p.k, Ef, rng, kNew, hNew, &bNew)) {
                p.k = kNew; p.hint = hNew; p.band = bNew;
                p.energy = bs->getEnergy(p.k, p.band, p.hint);
              }
            } else {
              nSelf.fetch_add(1, std::memory_order_relaxed);
            }
          }
          if (removed) {
            const SizeType idxCont =
                surface.getOhmicContactIdx(Base::device.posToCoord(pos));
#pragma omp critical
            {
              nettoNrPart[idxType][idxCont]++;
              idxToRemove.push_back(ip);
            }
          }
        }
#pragma omp critical
        {
          for (int c = 0; c < 3; c++) wV[c] += tV[c];
          wE += tE; wSum += tS;
          for (int b = 0; b < 8; b++) wBand[b] += tB[b];
          if (pxV.empty()) { pxV.assign(deviceExtent[0], 0.0); pxE = pxV; pxN = pxV; }
          for (SizeType i = 0; i < deviceExtent[0]; i++) { pxV[i] += pV[i]; pxE[i] += pE[i]; pxN[i] += pN[i]; }
        }
      }
      removeParticles(idxType, idxToRemove);
    }
    return nettoNrPart;
  }

  void assignParticlesToMesh(SizeType idxType, emcGrid<T, Dim> &gridNrParticles) {
    Base::pmScheme.assignToMesh(positions[idxType], Base::nrCarriersPerPart,
                                deviceSpacing, gridNrParticles);
  }

  /// same reservoir logic as emcBasicParticleHandler, on our containers
  NettoParticleCounter handleOhmicContacts() {
    SizeVec coord;
    auto nrRemPart = Base::initNettoParticleCounter();
    auto nrInjPart = Base::initNettoParticleCounter();
    auto &surface = Base::device.getSurface();
    emcGrid<T, Dim> nrPart(deviceExtent, 0);
    for (const auto &[idxType, partType] : Base::idxTypeToPartType) {
      if (!partType->isInjected()) continue;
      nrPart.fill(0);
      for (SizeType ip = 0; ip < getNrParticles(idxType); ip++) {
        coord = Base::device.posToCoord(positions[idxType][ip]);
        if (surface.isReservoirContact(coord)) {
          if (nrPart[coord] < Base::expNrPart[idxType][coord]) {
            nrPart[coord] += Base::nrCarriersPerPart;
          } else {
            removeParticle(idxType, ip);
            ip--;
            nrRemPart[idxType][surface.getOhmicContactIdx(coord)]--;
          }
        }
      }
      nrInjPart[idxType] = Base::generateInjectedParticles(idxType, nrPart);
      std::transform(nrRemPart[idxType].begin(), nrRemPart[idxType].end(),
                     nrInjPart[idxType].begin(), nrRemPart[idxType].begin(),
                     std::plus<int>());
    }
    return nrRemPart;
  }

  void print(std::string namePrefix, std::string nameSuffix) {
    for (const auto &type : Base::idxTypeToPartType) {
      std::ofstream os(namePrefix + type.second->getName() + nameSuffix + ".txt");
      os << Base::device.getMaxPos() << "\n";
      const auto n = positions[type.first].size();
      for (SizeType ip = 0; ip < n; ip++) {
        os << ip << " " << positions[type.first][ip];
        if (type.second->isMoved()) {
          const auto &p = particles[type.first][ip];
          os << " " << p.k[0] << " " << p.k[1] << " " << p.k[2] << " "
             << p.energy - cbm << " " << p.band;
        }
        if (ip + 1 < n) os << "\n";
      }
    }
  }

  /// per-cell observables through emcSimulationResults: instantaneous means
  /// over the particles in each cell (physical velocity [m/s], energy above
  /// the band edge [eV], share of carriers in bands >= 1), zero in empty
  /// cells; the results object averages them over the final-average steps
  /// and writes <prefix>FB{Vx,Vy,Energy,BandShare}Avg.txt.
  std::vector<std::string> fieldNames() const override {
    return {"FBVx", "FBVy", "FBEnergy", "FBBandShare"};
  }
  void fillFields(std::map<std::string, emcGrid<T, Dim>> &fields) override {
    auto &fvx = fields.at("FBVx"); auto &fvy = fields.at("FBVy");
    auto &fe = fields.at("FBEnergy"); auto &fb = fields.at("FBBandShare");
    emcGrid<T, Dim> cnt(deviceExtent, 0);
    fvx.fill(0); fvy.fill(0); fe.fill(0); fb.fill(0);
    for (const auto &type : Base::idxTypeToPartType) {
      if (!type.second->isMoved()) continue;
      const SizeType it = type.first;
      for (SizeType i = 0; i < getNrParticles(it); i++) {
        auto &p = particles[it][i];
        const auto c = Base::device.posToCoord(positions[it][i]);
        const auto v = bs->getVelocity(p.k, p.band, p.hint);
        fvx[c] += carrierSign * v[0]; fvy[c] += carrierSign * v[1];
        fe[c] += p.energy - cbm; fb[c] += p.band > 0 ? T(1) : T(0); cnt[c] += 1;
      }
    }
    SizeVec c;
    for (c.fill(0); !cnt.isEndCoord(c); cnt.advanceCoord(c))
      if (cnt[c] > 0) { fvx[c] /= cnt[c]; fvy[c] /= cnt[c]; fe[c] /= cnt[c]; fb[c] /= cnt[c]; }
  }

  T getChannelDriftCurrent(SizeType idxType, T x0, T x1, T channelLength) const {
    T sumVx = 0;
    for (SizeType i = 0; i < getNrParticles(idxType); ++i) {
      const auto &pos = positions[idxType][i];
      if (pos[0] < x0 || pos[0] > x1) continue;
      auto &p = const_cast<FBParticle &>(particles[idxType][i]);
      sumVx += carrierSign * bs->getVelocity(p.k, p.band, p.hint)[0];
    }
    return constants::q * Base::nrCarriersPerPart * sumVx / channelLength;
  }

  /// ------------------------------------------------ run-level observables
  /// time-weighted mean velocity [m/s] (component c), energy above the edge
  /// [eV], band occupancies, and event counts - the bulk driver's numbers.
  T meanVelocity(int c) const { return wSum > 0 ? T(wV[c] / wSum) : T(0); }
  T meanEnergy() const { return wSum > 0 ? T(wE / wSum) : T(0); }
  /// x-profile at grid column i: time-weighted <vx> [m/s], <E>-CBM [eV], and
  /// the time-weighted particle count (density = count * carriersPerPart /
  /// cellVolume once divided by the averaging time)
  SizeType profileSize() const { return pxN.size(); }
  T profileVx(SizeType i) const { return pxN[i] > 0 ? T(pxV[i] / pxN[i]) : T(0); }
  T profileE(SizeType i) const { return pxN[i] > 0 ? T(pxE[i] / pxN[i]) : T(0); }
  T profileWeight(SizeType i) const { return T(pxN[i]); }
  T bandOccupancy(int b) const { return wSum > 0 && b < 8 ? T(wBand[b] / wSum) : T(0); }
  long long realEvents() const { return nReal; }
  long long selfEvents() const { return nSelf; }
  long long removedAtContacts() const { return nRemoved; }
  long long reflections() const { return nReflected; }
  const FBS &scattering(int b) const { return *scatBand[b]; }
  T carrier() const { return carrierSign; }   ///< +1 electrons, -1 holes
  T bandEdge() const { return cbm; }
  int bands() const { return nb; }

protected:
  /// thermal injection on the NUMERIC band: band by thermal weight, energy
  /// DOS-weighted, k on the iso-surface - the bulk driver's initialisation
  void addParticle(SizeType idxType, const SizeVec &coord, emcRNG &rng,
                   bool /*isInitial*/) {
    positions[idxType].push_back(
        initParticlePos(coord, deviceExtent, deviceSpacing, rng));
    if (!Base::idxTypeToPartType.at(idxType)->isMoved()) return;
    FBParticle p;
    std::uniform_real_distribution<T> U01(0, 1);
    if (nb > 1 && zTot > 0) {
      T u = U01(rng) * zTot, acc = 0;
      for (int b = 0; b < nb; b++) { acc += zB[b]; if (u <= acc) { p.band = b; break; } }
    }
    FBS &sc = *scatBand[p.band];
    const T kT = constants::kB * temperature / constants::q;
    T E;
    do E = sc.sampleThermalEnergy(kT, rng);
    while (!sc.sampleFinalState(E, rng, p.k, p.hint));
    p.energy = bs->getEnergy(p.k, p.band, p.hint);
    particles[idxType].push_back(p);
  }

private:
  bool handleBoundary(FBParticle &p, std::array<T, Dim> &pos,
                      const typename DeviceType::SurfaceType &surface) {
    bool out = false;
    for (SizeType d = 0; d < Dim; d++)
      if (pos[d] < 0 || pos[d] > deviceMaxPos[d]) out = true;
    if (!out) return false;
    std::array<T, Dim> clamped = pos;
    for (SizeType d = 0; d < Dim; d++)
      clamped[d] = std::max(T(0), std::min(pos[d], deviceMaxPos[d]));
    const auto coordB = Base::device.posToCoord(clamped);
    if (surface.isOhmicContact(coordB)) {
      pos = clamped;
      nRemoved.fetch_add(1, std::memory_order_relaxed);
      return true;
    }
    // specular reflection on every crossed face
    for (SizeType d = 0; d < Dim; d++) {
      if (pos[d] < 0) { pos[d] = -pos[d]; p.k[d] = -p.k[d]; }
      else if (pos[d] > deviceMaxPos[d]) { pos[d] = 2 * deviceMaxPos[d] - pos[d]; p.k[d] = -p.k[d]; }
      pos[d] = std::max(T(0), std::min(pos[d], deviceMaxPos[d]));
    }
    nReflected.fetch_add(1, std::memory_order_relaxed);
    return false;
  }

  void removeParticle(SizeType idxType, SizeType ip) {
    positions[idxType].erase(positions[idxType].begin() + ip);
    if (Base::idxTypeToPartType.at(idxType)->isMoved())
      particles[idxType].erase(particles[idxType].begin() + ip);
  }
  void removeParticles(SizeType idxType, std::vector<SizeType> idx) {
    std::sort(idx.begin(), idx.end());
    const bool moved = Base::idxTypeToPartType.at(idxType)->isMoved();
    for (SizeType i = 0; i < idx.size(); i++) {
      const SizeType cur = idx[i] - i;
      positions[idxType].erase(positions[idxType].begin() + cur);
      if (moved) particles[idxType].erase(particles[idxType].begin() + cur);
    }
  }

  void loadPackage() {
    if (packagePath.empty())
      throw std::runtime_error("emcFullBandParticleHandler: set packagePath first");
    bs.reset(new NBS(packagePath));
    cbm = bs->getBandMinimum(0);
    {   // carrier attribute on /bands/electron ("hole" -> flipped axis)
      hid_t f = H5Fopen(packagePath.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
      if (f >= 0) {
        hid_t g = H5Gopen2(f, "/bands/electron", H5P_DEFAULT);
        if (g >= 0) {
          if (H5Aexists(g, "carrier") > 0) {
            hid_t a = H5Aopen(g, "carrier", H5P_DEFAULT); hid_t t = H5Aget_type(a);
            char buf[64] = {0};
            if (H5Tis_variable_str(t)) { char *q = nullptr; if (H5Aread(a, t, &q) >= 0 && q) { std::snprintf(buf, sizeof buf, "%s", q); free(q); } }
            else H5Aread(a, t, buf);
            if (std::string(buf).find("hole") != std::string::npos) carrierSign = T(-1);
            H5Tclose(t); H5Aclose(a);
          }
          H5Gclose(g);
        }
        H5Fclose(f);
      }
    }
    nb = nrBands > 0 ? nrBands : static_cast<int>(bs->getNrBands());
    for (int b = 0; b < nb; b++)
      scatBand.emplace_back(new FBS(packagePath, *bs, b, temperature,
                                    cbm + T(1.0), binWidth, mode));
    // doping -> impurity ladder level (one region for milestone 1)
    T N = dopingCm3;
    if (N <= 0) {
      const auto &dp = Base::device.getDopingProfile();
      SizeVec c; c.fill(0);
      T mx = 0;
      for (; !Base::device.isEndCoord(c); Base::device.advanceCoord(c))
        mx = std::max(mx, std::fabs(dp.getDoping(c)));
      N = mx * T(1e-6);                              // m^-3 -> cm^-3
    }
    for (auto &sb : scatBand) {
      if (sb->hasDopingLadder()) sb->setDoping(N);   // global default level
      if (useSlabs) sb->buildGamma0Slabs(slabWidth); // bounds the WORST level
      gamma0B.push_back(sb->getGamma0());
      slabWB.push_back(sb->getSlabWidth());
    }
    // ladder level per DOPING REGION: scan the grid once, record each region's
    // |doping|, map it to the nearest level. Particles read their region's
    // level every scattering decision (the region index is refreshed from the
    // position), so an n+/n/n+ profile scatters at three different strengths
    // while one Gamma0 - built over the worst level - bounds all of them.
    if (scatBand[0]->hasDopingLadder()) {
      const auto &dp = Base::device.getDopingProfile();
      std::vector<T> regDoping;
      SizeVec c; c.fill(0);
      for (; !Base::device.isEndCoord(c); Base::device.advanceCoord(c)) {
        const int r = dp.getDopingRegionIdx(c);
        if (r < 0) continue;
        if ((int)regDoping.size() <= r) regDoping.resize(r + 1, T(0));
        regDoping[r] = std::max(regDoping[r], std::fabs(dp.getDoping(c)));
      }
      levelOfRegion.resize(regDoping.size());
      for (std::size_t r = 0; r < regDoping.size(); r++)
        levelOfRegion[r] = scatBand[0]->levelFor(regDoping[r] * T(1e-6));
      std::printf("# doping regions -> ladder levels:");
      for (std::size_t r = 0; r < regDoping.size(); r++)
        std::printf("  region %zu: %.3g cm^-3 -> level %d (%.3g)", r, regDoping[r] * 1e-6,
                    levelOfRegion[r], scatBand[0]->getDopingLadder()[levelOfRegion[r]]);
      std::printf("\n");
    }
    const T kT = constants::kB * temperature / constants::q;
    for (int b = 0; b < nb; b++) {
      vMaxB.push_back(bs->getMaxSpeed(b));
      zB.push_back(scatBand[b]->thermalWeight(kT)); zTot += zB.back();
    }
    std::printf("# full-band handler: %s  %s  %d band(s)  T=%.0f K  doping level %.3g cm^-3%s\n",
                packagePath.c_str(), carrierSign < 0 ? "HOLES (flipped axis, v_physical = -v_stored)" : "electrons", nb, temperature,
                scatBand[0]->hasDopingLadder() ? scatBand[0]->getActiveDoping() : N,
                scatBand[0]->hasDopingLadder() ? " (ladder)" : " (package as built)");
  }
};

#endif
