#ifndef EMC_BULK_AVERAGES_HPP
#define EMC_BULK_AVERAGES_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include <emcUtil.hpp>

/*! \brief Time-weighted bulk-transport moments and the observables derived
 * from them.
 *
 * WHY THIS EXISTS. ViennaEMC's `emcSimulationResults` is device-side: contact
 * currents, carrier concentrations, potential. There was no bulk counterpart,
 * so every bulk observable lived inline in `examples/fullBandBulk` - each new
 * one meant editing the example, and diffusion / MSD / noise could not be
 * computed at all.
 *
 * PORTED FROM ViennaMC `src/averages/` IN STRUCTURE, NOT IN RUNTIME.
 * The donor keeps a string-keyed map of `Average` objects and, per sample,
 * makes a virtual call that returns `std::vector<Quantity>` BY VALUE (a heap
 * allocation), then looks the average up by name; its `AverageData`
 * constructor also appends to a debug file on every construction. ViennaMC is
 * serial and never noticed. Here the samples arrive from
 *   for (particle) -> while (free flight) -> per substep
 * millions of times per replica, inside OpenMP particle loops, so that design
 * is not portable as written. What is kept is the part that carries the
 * value: named observables, moments accumulated ONCE and shared, and derived
 * outputs computed centrally from them rather than re-derived per caller.
 * What is dropped: virtual dispatch, allocation and unit objects on the hot
 * path. `addSample` is a handful of fused multiply-adds and allocates nothing.
 *
 * The donor's selective ACTIVATION is dropped with it, deliberately: it existed
 * because each average cost a virtual call and an allocation. Here the moments
 * are six additions, so computing them unconditionally is cheaper than testing
 * whether to. Only the energy histogram, which costs MEMORY rather than time,
 * stays opt-in.
 *
 * UNITS ARE THE CALLER'S AND ARE NOT CONVERTED. Two different Boltzmann
 * constants are live in this codebase - `emcConstants::kB` is 1.38066e-23 J/K
 * while the full-band example uses 8.617333262e-5 eV/K - so this class takes
 * `kT` as an argument in the SAME units as the energies it was fed, and never
 * reaches for a constant of its own. The full-band convention is energies in
 * eV, velocities in m/s, rates in 1/s, times in s.
 *
 * @param nPoints    grid points to resolve (1 for a homogeneous bulk run)
 * @param nEnergyBins  0 disables the energy histogram
 * @param energyBinWidth  bin width, in the caller's energy units
 */
template <class T> class emcBulkAverages {
private:
  SizeType nPoints;
  SizeType nBins;
  T binWidth;

  /// zeroth moment: total weight (sum of dt). The normaliser for every other
  /// moment, and the donor's "M0".
  std::vector<T> m0;
  /// first moment of velocity along the field direction ("M1v")
  std::vector<T> m1v;
  /// first moment of energy ("M2E" in the donor's naming, which counts the
  /// power of the VELOCITY the quantity came from, not the moment order)
  std::vector<T> m1E;
  /// second moment of velocity: <|v|^2>, the final-state measure check
  std::vector<T> m2v;
  /// <v^2 / Gamma>, the SERTA integrand evaluated on the MC's own ensemble
  std::vector<T> m2vOverRate;
  /// first moment of the total scattering rate
  std::vector<T> m1Rate;
  /// flat [point * nBins + bin] energy histogram, empty when disabled
  std::vector<T> hist;

public:
  emcBulkAverages(SizeType inNPoints = 1, SizeType inNEnergyBins = 0,
                  T inEnergyBinWidth = 0)
      : nPoints(std::max<SizeType>(1, inNPoints)), nBins(inNEnergyBins),
        binWidth(inEnergyBinWidth), m0(nPoints, 0), m1v(nPoints, 0),
        m1E(nPoints, 0), m2v(nPoints, 0), m2vOverRate(nPoints, 0),
        m1Rate(nPoints, 0), hist(nBins ? nPoints * nBins : 0, 0) {}

  /*! \brief Accumulate one time-weighted sample. HOT PATH.
   *
   * Everything is weighted by `weight` (the free-flight time slice), which is
   * what makes these time averages rather than per-event averages: a state the
   * carrier occupies for longer counts proportionally more.
   *
   * @param vProjected velocity projected onto the field direction
   * @param energy     energy measured from the band edge
   * @param vSquared   |v|^2 at the same point `rate` was evaluated
   * @param rate       total scattering rate; <=0 is skipped for the SERTA
   *                   moment only, since 1/Gamma is undefined there
   * @param weight     time slice this sample represents
   * @param point      grid point index (0 for bulk)
   */
  inline void addSample(T vProjected, T energy, T vSquared, T rate, T weight,
                        SizeType point = 0) noexcept {
    m0[point] += weight;
    m1v[point] += vProjected * weight;
    m1E[point] += energy * weight;
    m2v[point] += vSquared * weight;
    m1Rate[point] += rate * weight;
    if (rate > 0)
      m2vOverRate[point] += vSquared / rate * weight;
    if (nBins) {
      const std::ptrdiff_t b = static_cast<std::ptrdiff_t>(energy / binWidth);
      if (b >= 0 && b < static_cast<std::ptrdiff_t>(nBins))
        hist[point * nBins + static_cast<SizeType>(b)] += weight;
    }
  }

  /*! \brief Fold another accumulator in, for per-thread or per-replica
   * reduction. Every moment is a plain sum, so merging N thread-local
   * accumulators reproduces one serial accumulation to floating-point
   * roundoff - agreement is to ~1e-12 relative, NOT bit-exact, because the
   * summation order differs. Do not build a bit-reproducibility claim on it;
   * reproducibility of this engine is asserted over identical binary + inputs
   * + thread count, which fixes the order.
   */
  void merge(const emcBulkAverages<T> &other) {
    if (other.nPoints != nPoints)
      return;
    for (SizeType i = 0; i < nPoints; i++) {
      m0[i] += other.m0[i];
      m1v[i] += other.m1v[i];
      m1E[i] += other.m1E[i];
      m2v[i] += other.m2v[i];
      m2vOverRate[i] += other.m2vOverRate[i];
      m1Rate[i] += other.m1Rate[i];
    }
    if (nBins && other.nBins == nBins)
      for (SizeType i = 0; i < hist.size(); i++)
        hist[i] += other.hist[i];
  }

  void reset() {
    std::fill(m0.begin(), m0.end(), T{0});
    std::fill(m1v.begin(), m1v.end(), T{0});
    std::fill(m1E.begin(), m1E.end(), T{0});
    std::fill(m2v.begin(), m2v.end(), T{0});
    std::fill(m2vOverRate.begin(), m2vOverRate.end(), T{0});
    std::fill(m1Rate.begin(), m1Rate.end(), T{0});
    std::fill(hist.begin(), hist.end(), T{0});
  }

  SizeType getNrPoints() const { return nPoints; }
  SizeType getNrEnergyBins() const { return nBins; }
  T getTotalWeight(SizeType pt = 0) const { return m0[pt]; }

  // --- derived observables -------------------------------------------------
  // Each is a pure function of the moments above. They are gathered here so a
  // second caller cannot re-derive one slightly differently, which is exactly
  // how the example accumulated its own private copies.

  /// <v.F^> [caller's velocity units]
  T getDriftVelocity(SizeType pt = 0) const {
    return m0[pt] > 0 ? m1v[pt] / m0[pt] : T{0};
  }
  /// <E - E_edge>
  T getMeanEnergy(SizeType pt = 0) const {
    return m0[pt] > 0 ? m1E[pt] / m0[pt] : T{0};
  }
  /// <|v|^2>
  T getMeanSquaredVelocity(SizeType pt = 0) const {
    return m0[pt] > 0 ? m2v[pt] / m0[pt] : T{0};
  }
  /// <Gamma>
  T getMeanScatterRate(SizeType pt = 0) const {
    return m0[pt] > 0 ? m1Rate[pt] / m0[pt] : T{0};
  }
  /*! \brief Carrier temperature from the mean energy, 2<E>/(3 kB).
   * `kB` must be in the caller's energy units per kelvin. */
  T getCarrierTemperature(T kB, SizeType pt = 0) const {
    return kB > 0 ? 2 * getMeanEnergy(pt) / (3 * kB) : T{0};
  }

  /*! \brief Drift mobility, from the drift estimator.
   *
   * SIGN IS BY CARRIER, and this is not cosmetic. Electrons drift AGAINST the
   * field, so -vd/F is positive for them. A hole package stores a flipped axis
   * (E = VBM - E_e) and its carriers drift WITH the field, so the same formula
   * would report a negative mobility. Pass +1 for electrons, -1 for holes -
   * the value the package's carrier attribute implies.
   *
   * @param field       field strength [V/cm], matching `scale`
   * @param carrierSign +1 electrons, -1 holes
   * @param scale       unit conversion; 1e4 takes m^2/Vs to cm^2/Vs
   */
  T getMobility(T field, T carrierSign = 1, T scale = 1e4,
                SizeType pt = 0) const {
    return field != 0 ? -carrierSign * getDriftVelocity(pt) / field * scale
                      : T{0};
  }

  /*! \brief Mobility from the SERTA integrand, mu = <v^2/Gamma>/(3 kT).
   *
   * An INDEPENDENT estimator of the same quantity: same states, same rates, no
   * drift estimator involved. Comparing it against getMobility() is the
   * check that has repeatedly caught real defects - a corrupt valence package
   * split the two by 19x, where healthy silicon packages agree to 3-4%.
   *
   * @param kT in the SAME energy units as the energies fed to addSample
   */
  T getMobilitySERTA(T kT, T scale = 1e4, SizeType pt = 0) const {
    if (m0[pt] <= 0 || kT <= 0)
      return T{0};
    return m2vOverRate[pt] / m0[pt] / (3 * kT) * scale;
  }

  /*! \brief SERTA mobility on the ensemble's OWN temperature.
   *
   * getMobilitySERTA divides by 3kT of the LATTICE, which is only the right
   * normalisation if the ensemble is thermal at that temperature. It is not
   * when the band-edge DOS is under-resolved: the missing low-energy states
   * push the equilibrium distribution up, the ensemble sits hot at ANY field
   * (measured <E> = 0.0446 eV at 25 V/cm against 1.5kT = 0.0388 on a REFINE=6
   * silicon package), and the lattice-kT SERTA inflates by T_eff / T_lattice.
   * That alone accounted for most of a 16-20% MC-vs-SERTA gap (2026-09-02).
   *
   * This variant replaces 3kT with 2<E>, i.e. the ensemble's own effective
   * temperature, so the two estimators are compared on the SAME distribution.
   * The residual between it and getMobility() is then transport physics
   * (in-scattering, inelasticity, momentum retention), not a normalisation.
   */
  T getMobilitySERTASelfConsistent(T scale = 1e4, SizeType pt = 0) const {
    const T e = getMeanEnergy(pt);
    if (m0[pt] <= 0 || e <= 0)
      return T{0};
    return m2vOverRate[pt] / m0[pt] / (2 * e) * scale;
  }

  /// time-weighted energy histogram for one point; empty when disabled
  std::vector<T> getEnergyHistogram(SizeType pt = 0) const {
    if (!nBins)
      return {};
    return std::vector<T>(hist.begin() + pt * nBins,
                          hist.begin() + (pt + 1) * nBins);
  }
};

/*! \brief Mean-squared displacement at fixed checkpoint times, and the
 * diffusion coefficient and Einstein mobility derived from it.
 *
 * WHY. This is a THIRD, independent route to the mobility: no drift under a
 * field, and no SERTA normalisation by a temperature. Both of the other two
 * were shown on 2026-09-02 to inherit the mesh's band-edge error through
 * exactly those channels. D from the slope of the position variance depends
 * on neither, and has an exact answer on the parabolic test package through
 * the Einstein relation, D = mu kT / e.
 *
 * PORTED FROM ViennaMC `AverageMSD` / `AverageDiffusion` (dev_wendelin_final)
 * in intent: MSD against each trajectory's own reference position, D from
 * the time derivative of the variance about the mean (the donor's second
 * method, "better for shorter simulation times"). The donor sampled every
 * particle's position list per step; here the caller records one displacement
 * per particle per CHECKPOINT time, which is what makes it usable inside the
 * event-driven flight loop, and D is a least-squares slope over the
 * checkpoints rather than a single-time ratio, which removes the ballistic
 * transient.
 *
 * Units are the caller's: displacements in m and times in s give D in m^2/s;
 * getEinsteinMobility takes kT in the same energy unit the caller uses (eV
 * here) and returns m^2/Vs times `scale`.
 *
 * @param nTimes       number of checkpoints, index 0 being the reference
 * @param tStart       time of checkpoint 0
 * @param dtSample     checkpoint spacing
 */
template <class T> class emcDisplacementStatistics {
private:
  SizeType nTimes;
  T tStart, dtSample;
  std::vector<std::array<T, 3>> s1, s2;   // sum dr, sum dr^2 per checkpoint
  std::vector<T> n;                        // particles counted per checkpoint

public:
  emcDisplacementStatistics(SizeType inNTimes = 16, T inTStart = 0,
                            T inDtSample = 1)
      : nTimes(std::max<SizeType>(2, inNTimes)), tStart(inTStart),
        dtSample(inDtSample), s1(nTimes, {T{0}, T{0}, T{0}}),
        s2(nTimes, {T{0}, T{0}, T{0}}), n(nTimes, T{0}) {}

  SizeType getNrTimes() const { return nTimes; }
  T getTime(SizeType i) const { return tStart + static_cast<T>(i) * dtSample; }

  /// displacement of one particle from its own reference, at checkpoint i
  inline void addSample(SizeType i, const T dr[3]) noexcept {
    for (int c = 0; c < 3; c++) {
      s1[i][c] += dr[c];
      s2[i][c] += dr[c] * dr[c];
    }
    n[i] += 1;
  }

  void merge(const emcDisplacementStatistics<T> &o) {
    if (o.nTimes != nTimes)
      return;
    for (SizeType i = 0; i < nTimes; i++) {
      for (int c = 0; c < 3; c++) {
        s1[i][c] += o.s1[i][c];
        s2[i][c] += o.s2[i][c];
      }
      n[i] += o.n[i];
    }
  }

  T getMeanDisplacement(SizeType i, int c) const {
    return n[i] > 0 ? s1[i][c] / n[i] : T{0};
  }
  /// variance about the mean - the drift (vd t) is removed, so this is the
  /// diffusive spread even under a field
  T getVariance(SizeType i, int c) const {
    if (n[i] <= 0)
      return T{0};
    const T m = s1[i][c] / n[i];
    return s2[i][c] / n[i] - m * m;
  }

  /*! \brief Diffusion coefficient along component c: half the least-squares
   * slope of the variance against time over checkpoints [iFirst, nTimes).
   * Skipping the first checkpoints removes the ballistic regime (variance
   * grows as t^2 for t below the momentum relaxation time). */
  T getDiffusion(int c, SizeType iFirst = 1) const {
    T st = 0, sv = 0, stt = 0, stv = 0, cnt = 0;
    for (SizeType i = std::max<SizeType>(1, iFirst); i < nTimes; i++) {
      if (n[i] <= 0)
        continue;
      const T t = getTime(i) - getTime(0), v = getVariance(i, c);
      st += t; sv += v; stt += t * t; stv += t * v; cnt += 1;
    }
    if (cnt < 2)
      return T{0};
    const T den = cnt * stt - st * st;
    return den != 0 ? (cnt * stv - st * sv) / den / 2 : T{0};
  }

  /*! \brief Einstein mobility mu = e D / kT. With kT in eV the charge cancels:
   * D [m^2/s] / kT [eV] is already m^2/Vs; `scale` = 1e4 gives cm^2/Vs. */
  static T getEinsteinMobility(T D, T kT, T scale = 1e4) {
    return kT > 0 ? D / kT * scale : T{0};
  }
};

/*! \brief Mean and standard error of a named quantity across replicas.
 *
 * The donor had no equivalent: ViennaMC reports one run. Independent replicas
 * are what turn a number into a number with an error bar, and several
 * conclusions here have rested on that - the interband result was called at
 * 12 sigma, and the anisotropy gate passes at 1.02 sigma, neither of which is
 * meaningful without this.
 *
 * Uses the standard error of the mean over replicas, which makes no assumption
 * about the within-replica correlation structure - the samples inside one
 * replica are heavily autocorrelated along a trajectory, so a naive
 * sample-count SE would be badly optimistic.
 */
template <class T> class emcReplicaStatistics {
private:
  std::vector<std::string> names;
  std::vector<std::vector<T>> values;

public:
  /// record one replica's value of a named quantity
  void add(const std::string &name, T value) {
    for (SizeType i = 0; i < names.size(); i++) {
      if (names[i] == name) {
        values[i].push_back(value);
        return;
      }
    }
    names.push_back(name);
    values.push_back({value});
  }

  SizeType getNrReplicas(const std::string &name) const {
    for (SizeType i = 0; i < names.size(); i++)
      if (names[i] == name)
        return values[i].size();
    return 0;
  }

  T getMean(const std::string &name) const {
    for (SizeType i = 0; i < names.size(); i++) {
      if (names[i] == name) {
        if (values[i].empty())
          return T{0};
        T s{0};
        for (T v : values[i])
          s += v;
        return s / static_cast<T>(values[i].size());
      }
    }
    return T{0};
  }

  /// standard error of the mean; 0 with fewer than two replicas, which is
  /// honest rather than a placeholder - one replica carries no error estimate
  T getStandardError(const std::string &name) const {
    for (SizeType i = 0; i < names.size(); i++) {
      if (names[i] == name) {
        const SizeType n = values[i].size();
        if (n < 2)
          return T{0};
        const T mean = getMean(name);
        T var{0};
        for (T v : values[i])
          var += (v - mean) * (v - mean);
        return std::sqrt(var / static_cast<T>(n - 1) / static_cast<T>(n));
      }
    }
    return T{0};
  }

  const std::vector<std::string> &getNames() const { return names; }
};

#endif // EMC_BULK_AVERAGES_HPP
