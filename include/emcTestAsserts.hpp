#ifndef EMC_TEST_HPP
#define EMC_TEST_HPP

#include <math.h>
#include <string>

#ifdef _MSC_VER
#define __PRETTY_FUNCTION__ __FUNCSIG__
#endif

#define EMCTEST_ASSERT(condition)                                              \
  {                                                                            \
    if (!(condition)) {                                                        \
      throw std::runtime_error(std::string(__FILE__) + std::string(":") +      \
                               std::to_string(__LINE__) +                      \
                               std::string(" in ") +                           \
                               std::string(__PRETTY_FUNCTION__) +              \
                               std::string(" Condition not fulfilled"));       \
    }                                                                          \
  }

#define EMCTEST_ASSERT_ISCLOSE(first, second, eps)                             \
  {                                                                            \
    if ((std::fabs(first - second) > eps)) {                                   \
      throw std::runtime_error(                                                \
          std::string(__FILE__) + std::string(":") +                           \
          std::to_string(__LINE__) + std::string(" in ") +                     \
          std::string(__PRETTY_FUNCTION__) +                                   \
          std::string(" Numbers not close ") + std::to_string(first) +         \
          std::string(" ") + std::to_string(second));                          \
    }                                                                          \
  }

#endif // EMC_TEST_HPP