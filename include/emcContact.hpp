#ifndef EMC_CONTACT_HPP
#define EMC_CONTACT_HPP

#include <emcMessage.hpp>
#include <emcUtil.hpp>

/// different types of contacts
enum struct emcContactType : SizeType { OHMIC, SCHOTTKY, GATE };

/*! \brief Base Class for Contacts.
 *
 * @tparam T Numeric Type
 * @param appliedVoltage applied voltage [V]
 */
template <class T> class emcContact {
protected:
  T appliedVoltage; // voltage applied at contact [V]

public:
  emcContact(T inAppliedVoltage) : appliedVoltage(inAppliedVoltage) {}

  virtual ~emcContact() = default;

  T getAppliedVoltage() const { return appliedVoltage; }

  /// returns further parameter for that contact
  virtual T getFurtherParameter(SizeType idxParameter) const = 0;

  /// returns the contactType
  virtual emcContactType getType() const = 0;
};

/*! \brief Class for Ideal Ohmic Contact.
 *
 * It is assumed that there is no voltage drop at the ideal ohmic contact, for
 * that assumption to be reasonable we need to make sure that
 * the contact is in thermal equilibrium (free-carrier concentration
 * near contact should be held constant).
 *
 */
template <class T> class emcOhmicContact : public emcContact<T> {
public:
  emcOhmicContact(T inAppliedVoltage) : emcContact<T>(inAppliedVoltage) {}

  T getFurtherParameter(SizeType idxParameter) const {
    emcMessage::getInstance()
        .addError("Ohmic Contact has no further parameter.")
        .print();
    return 0;
  }

  emcContactType getType() const { return emcContactType::OHMIC; }
};

/*! \brief Class for Gate Contact.
 *
 * Gate Contact represents a dielectric material at the  boundary
 * with a given thickness and an ohmic contact on top of that
 * material. The contact is modelled as a Schottky Contact
 * with a known barrier height.
 *
 * @param epsROxide relative dielectric constant of oxide
 * @param thickness thickness of the oxide [m]
 * @param barrierHeight barrierHeight that is assumed at contact [V]
 *
 */
template <class T> class emcGateContact : public emcContact<T> {
  T epsROxide;
  T thickness;
  T barrierHeight;

public:
  emcGateContact(T inAppliedVoltage, T inEpsR, T inThickness, T inBarrierHeight)
      : emcContact<T>(inAppliedVoltage), epsROxide(inEpsR),
        thickness(inThickness), barrierHeight(inBarrierHeight) {}

  /// returns oxide epsROxide / thickness dependent
  /// on given idxParameter
  T getFurtherParameter(SizeType idxParameter) const {
    if (idxParameter == 0)
      return epsROxide;
    else if (idxParameter == 1)
      return thickness;
    else if (idxParameter == 2)
      return barrierHeight;
    else {
      emcMessage::getInstance()
          .addError("Gate Contact only has 3 further parameter.")
          .print();
    }
    return 0;
  }

  emcContactType getType() const { return emcContactType::GATE; }
};

/// Class for Schottky Contact (not implemented yet)
template <class T> class emcSchottkyContact : public emcContact<T> {

public:
  emcSchottkyContact(T inAppliedVoltage) : emcContact<T>(inAppliedVoltage) {}

  T getFurtherParameter(SizeType idxParameter) const {
    emcMessage::getInstance()
        .addError("Schottky Contact has no further parameter.")
        .print();
    return 0;
  }

  emcContactType getType() const { return emcContactType::SCHOTTKY; }
};

#endif // EMC_CONTACT_HPP