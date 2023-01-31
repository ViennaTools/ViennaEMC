#ifndef EMC_UTIL_HPP
#define EMC_UTIL_HPP

#include <algorithm>
#include <array>
#include <iostream>
#include <random>
#include <type_traits>
#include <vector>

typedef size_t SizeType;

typedef std::mt19937_64 emcRNG; //!< random number generator type

/// cast enum to underlying data structure
template <typename Enum>
constexpr typename std::underlying_type<Enum>::type toUnderlying(Enum e) {
  return static_cast<typename std::underlying_type<Enum>::type>(e);
}

// ----- array operations ----------------------------------------------
/// performs the inner product of an array with itself
template <class T, SizeType Dim> T square(const std::array<T, Dim> &vec) {
  T res = 0;
  std::for_each(vec.begin(), vec.end(),
                [&res](const T &entry) { res += entry * entry; });
  return res;
}

/// calculates the euclidian/L2 norm of an array
template <class T, SizeType Dim> T norm(const std::array<T, Dim> &vec) {
  return std::sqrt(square(vec));
}

/// returns a normalized array, that points in the same direction than
/// input
template <class T, SizeType Dim> void normalize(std::array<T, Dim> &vec) {
  T normL2 = norm(vec);
  std::for_each(vec.begin(), vec.end(),
                [&normL2](T &entry) { entry /= normL2; });
}

/// calculates the inner product between two 3-dimensional arrays
template <class T>
T innerProduct(const std::array<T, 3> &lvec, const std::array<T, 3> &rvec) {
  return lvec[0] * rvec[0] + lvec[1] * rvec[1] + lvec[2] * rvec[2];
}

/// multiplies each component of the given array with factor
template <class T, SizeType Dim>
std::array<T, Dim> scale(const std::array<T, Dim> &vec, T factor) {
  std::array<T, Dim> res;
  std::transform(vec.begin(), vec.end(), res.begin(),
                 [&factor](const T &val) -> T { return val * factor; });
  return res;
}

/// calculates the component-wise sum of the two arrays
template <class T, SizeType Dim>
std::array<T, Dim> add(const std::array<T, Dim> &lvec,
                       const std::array<T, Dim> &rvec) {
  std::array<T, Dim> res;
  std::transform(lvec.begin(), lvec.end(), rvec.begin(), res.begin(),
                 [](const T &lval, const T &rval) -> T { return lval + rval; });
  return res;
}

/// calculates the component-wise subtraction (lvec - rvec) of the two arrays
template <class T, SizeType Dim>
std::array<T, Dim> subtract(const std::array<T, Dim> &lvec,
                            const std::array<T, Dim> &rvec) {
  std::array<T, Dim> res;
  std::transform(lvec.begin(), lvec.end(), rvec.begin(), res.begin(),
                 [](const T &lval, const T &rval) -> T { return lval - rval; });
  return res;
}

/// writes input-array into ofstream (whitespace between components)
template <class T, SizeType Dim>
std::ostream &operator<<(std::ostream &os, const std::array<T, Dim> &vec) {
  bool printWS = false;
  std::for_each(vec.begin(), vec.end(), [&printWS, &os](const T &val) {
    if (printWS)
      os << " ";
    os << val;
    printWS = true;
  });
  return os;
}

/// writes input-vector into ofstream (whitespace between components)
template <class T>
std::ostream &operator<<(std::ostream &os, const std::vector<T> &vec) {
  bool printWS = 0;
  std::for_each(vec.begin(), vec.end(), [&printWS, &os](const T &val) {
    if (printWS)
      os << " ";
    os << val;
    printWS = true;
  });
  return os;
}

// ----- additional utility functions -----------------------------
/// calculates the nr. of grid points of a grid with the given spacing,
/// that starts at origin and goes on until maxPos
template <class T, SizeType Dim>
std::array<SizeType, Dim> maxPosToExtent(const std::array<T, Dim> &maxPos,
                                         const std::array<T, Dim> &spacing) {
  std::array<SizeType, Dim> gridExtent;
  std::transform(maxPos.begin(), maxPos.end(), spacing.begin(),
                 gridExtent.begin(),
                 [&spacing](const T &pos, const T &delta) -> SizeType {
                   return std::round(pos / delta) + 1;
                 });
  return gridExtent;
}

/// initialize vector in random direction with given norm
/// using the 2 given random nrs between 0 and 1
template <class T>
std::array<T, 3> initRandomDirection(T norm, T rand1, T rand2) {
  T phi = 2 * M_PI * rand1;
  T cosTheta = 1 - 2 * rand2;
  std::array<T, 3> res;
  res[0] = norm * std::sqrt(1 - cosTheta * cosTheta) * std::cos(phi);
  res[1] = norm * std::sqrt(1 - cosTheta * cosTheta) * std::sin(phi);
  res[2] = norm * cosTheta;
  return res;
}

template <class T>
std::array<T, 3>
initRandomDirectionWithRespectToCurrentK(const std::array<T, 3> &k, T cosTheta,
                                         T rand) {
  std::array<T, 3> res;
  /*=== Calculate the rotation angles ===*/
  T kxy = std::sqrt(k[0] * k[0] + k[1] * k[1]);
  T norm = std::sqrt(kxy * kxy + k[2] * k[2]);
  T ct0 = k[2] / norm;
  T st0 = kxy / norm;
  T cfi0 = k[0] / kxy;
  T sfi0 = k[1] / kxy;

  /*=== Randomize momentum in the rotated coordinate system ===*/
  T st = std::sqrt(1.0 - cosTheta * cosTheta);
  T phi = 2.0 * M_PI * rand;
  T kxp = norm * st * cos(phi);
  T kyp = norm * st * sin(phi);
  T kzp = norm * cosTheta;

  /*=== Return to the original coordinate system ===*/
  res[0] = kxp * cfi0 * ct0 - kyp * sfi0 + kzp * cfi0 * st0;
  res[1] = kxp * sfi0 * ct0 + kyp * cfi0 + kzp * sfi0 * st0;
  res[2] = -kxp * st0 + kzp * ct0;
  return res;
}

#endif // EMC_UTIL_HPP