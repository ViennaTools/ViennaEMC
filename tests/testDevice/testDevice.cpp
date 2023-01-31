#include <emcDevice.hpp>
#include <emcTestAsserts.hpp>

#include "testBasicDevice2D.hpp"
#include "testBasicDevice3D.hpp"

int main() {
  testBasicDevice2D(/*true*/);
  testBasicDevice3D(/*true*/);
  return 0;
}