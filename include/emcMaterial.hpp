#ifndef EMC_MATERIAL_HPP
#define EMC_MATERIAL_HPP

#include <emcUtil.hpp>

/// class that defines material constants of a semiconductor material
/// @tparam T Numeric Type
/// @param epsR relative dielectric constant
/// @param rho density [kg / m^3]
/// @param Ni intrinsic carrier density [ 1 / m^3]
/// @param velSound velocity of sound in material [m / s]
/// @param bandGap bandGap-energy [eV]
template <class T> class emcMaterial {
  T epsR;
  T rho;
  T Ni;
  T velSound;
  T bandGap;

public:
  emcMaterial() = delete;

  emcMaterial(T inEpsR, T inRho, T inNi, T inVelSound, T inBandGap)
      : epsR(inEpsR), rho(inRho), Ni(inNi), velSound(inVelSound),
        bandGap(inBandGap) {}

  T getEpsR() const { return epsR; }

  T getDielectricConstant() const { return epsR * constants::eps0; }

  T getRho() const { return rho; }

  T getNi() const { return Ni; }

  T getVelSound() const { return velSound; }

  T getBandGap() const { return bandGap; }
};

#endif // EMC_MATERIAL_HPP