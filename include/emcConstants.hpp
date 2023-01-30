#ifndef EMC_CONSTANTS_HPP
#define EMC_CONSTANTS_HPP

#include <math.h>

/// global physical constants using Si-Units
namespace constants {
/// Pi
constexpr double pi = M_PI;

/// elementary charge [C]
constexpr double q = 1.60219e-19;

/// boltzmann constant [m² kg s⁻² K⁻¹]
constexpr double kB = 1.38066e-23;

/// red. planck constant [J s]
constexpr double hbar = 1.05459e-34;

/// el. field constant [A s V⁻¹ m⁻¹]
constexpr double eps0 = 8.85419e-12;

/// mass of electron [kg]
constexpr double me = 9.11e-31;

/// coulomb constant
constexpr double ke = 1. / (4 * pi * eps0);

} // namespace constants

#endif // EMC_CONSTANTS_HPP