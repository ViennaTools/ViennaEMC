// from ViennaLS/lsMessage.hpp

#ifndef EMC_MESSAGE_HPP
#define EMC_MESSAGE_HPP

#include <iostream>

/// Singleton class for thread-safe logging.
class emcMessage {
  std::string message;

  bool error = false;
  const unsigned tabWidth = 4;

  emcMessage() {}

public:
  // delete constructors to result in better error messages by compilers
  emcMessage(const emcMessage &) = delete;
  void operator=(const emcMessage &) = delete;

  static emcMessage &getInstance() {
    static emcMessage instance;
    return instance;
  }

  emcMessage &add(std::string s) {
#pragma omp critical
    {
      message += /*"\n" +*/ std::string(tabWidth, ' ') + "WARNING: " + s + "\n";
    }
    return *this;
  }

  emcMessage &addWarning(std::string s) {
#pragma omp critical
    {
      message += /*"\n" +*/ std::string(tabWidth, ' ') + "WARNING: " + s + "\n";
    }
    return *this;
  }

  emcMessage &addError(std::string s, bool shouldAbort = true) {
#pragma omp critical
    { message += "\n" + std::string(tabWidth, ' ') + "ERROR: " + s + "\n"; }
    // always abort once error message should be printed
    error = true;
    // abort now if asked
    if (shouldAbort)
      print();
    return *this;
  }

  emcMessage &addDebug(std::string s) {
#pragma omp critical
    { message += std::string(tabWidth, ' ') + "DEBUG: " + s + "\n"; }
    return *this;
  }

  void print(std::ostream &out = std::cout) {
    out << message;
    message.clear();
    if (error)
      abort();
  }
};

#endif // EMC_MESSAGE_HPP