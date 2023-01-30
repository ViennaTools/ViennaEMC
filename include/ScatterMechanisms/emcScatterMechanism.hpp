#ifndef EMC_SCATTER_MECHANISM_HPP
#define EMC_SCATTER_MECHANISM_HPP

#include <ValleyTypes/emcAbstractValley.hpp>
#include <emcParticle.hpp>
#include <emcUtil.hpp>

#include <memory>

/*! \brief Base Scatter Mechanism that has to be used for the implementation
 * of a scatter mechanism.
 *
 * @param idxValley index of the band-diagram valley used for scattering
 * @param ptrValley vector that points to all the valleys of the current
 * particleType
 */
template <class T> struct emcScatterMechanism {
  typedef emcAbstractValley<T> AbstractValley;

protected:
  SizeType idxValley;
  std::vector<AbstractValley *> ptrValley;

public:
  explicit emcScatterMechanism(SizeType inIdxValley) : idxValley(inIdxValley) {}

  virtual ~emcScatterMechanism() = default;

  /// @brief returns the scatter rate of the corresponding scatter mechanism
  /// @param energy energy of the particle [eV]
  /// @param idxRegion index of the corresponding region
  /// @return rate of the scatter mechanism [1 / s]
  virtual T getScatterRate(T energy, SizeType idxRegion) const = 0;

  /// @brief function that is called when this scatter event happens.
  /// @param particle particle that is scattered
  /// @param rng random number generator
  virtual void scatterParticle(emcParticle<T> &particle, emcRNG &rng) const = 0;

  /// \brief name of the scatter mechanism (used in outputFiles)
  virtual std::string getName() const = 0;

  void setPtrValley(std::vector<std::unique_ptr<AbstractValley>> &inPtrValley) {
    for (auto &ptr : inPtrValley)
      ptrValley.emplace_back(ptr.get());
  }

  SizeType getIdxValley() const { return idxValley; }

  /// @brief function that can be used to check if the implemented mechanism
  /// is right.
  virtual void check() {}
};

#endif // EMC_SCATTER_MECHANISM_HPP