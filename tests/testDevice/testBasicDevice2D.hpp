#include <emcDevice.hpp>
#include <emcTestAsserts.hpp>

#include <iostream>

using NumType2D = double;
using MaterialType2D = emcMaterial<NumType2D>;
using DeviceType2D = emcDevice<NumType2D, 2>;
using ValueVec2D = DeviceType2D::ValueVec;
using SizeVec2D = DeviceType2D::SizeVec;
using SizeVecSurface2D = DeviceType2D::SizeVecSurface;

/*
Schematic of device that is tested:

  < OHMIC >< GATE >< OHMIC >
    0   1   2   3   4   5
  ------------------------
0 |          |           |
1 | Region 0 | Region 2  |
2 |          |           |
  ------------------------
3 |          |           |
4 | Region 1 | Region 3  |
5 |          |           |
  ------------------------
    0   1   2   3   4   5
  <        OHMIC         >

TESTS:
  - add regions
  - test grid of added regions
  - add contacts
  - test grid of added contacts
*/
void testBasicDevice2D(bool printGrids = false) {
  NumType2D epsilon = 1e-9;
  NumType2D electronMass = constants::me;

  // create silicon material
  MaterialType2D silicon{11.8, 2329, 1.45e16, 9040, 1.14};

  // create device structure (material + extent + spacing)
  ValueVec2D deviceMaxPos = {5, 10};
  DeviceType2D device{silicon, deviceMaxPos, {1, 2}};

  // TEST DOPINGREGIONS -------------------------------------------
  // add different doping regions
  device.addConstantDopingRegion({0, 0}, {3, 6}, 1e20);
  device.addConstantDopingRegion({0, 6}, {3, 10}, 2e20);
  device.addConstantDopingRegion({3, 0}, {5, 6}, 3e20);
  device.addConstantDopingRegion({3, 6}, {5, 10}, 4e20);

  // test: Doping Regions
  SizeVec2D coord;
  SizeType expectedIdx;
  NumType2D expectedDoping;
  for (coord.fill(0); !device.isEndCoord(coord); device.advanceCoord(coord)) {
    if (coord[0] < 3 && coord[1] < 3) {
      expectedIdx = 0;
      expectedDoping = 1e20;
    } else if (coord[0] < 3 && coord[1] >= 3) {
      expectedIdx = 1;
      expectedDoping = 2e20;
    } else if (coord[0] >= 3 && coord[1] < 3) {
      expectedIdx = 2;
      expectedDoping = 3e20;
    } else {
      expectedIdx = 3;
      expectedDoping = 4e20;
    }
    EMCTEST_ASSERT(device.getDopingProfile().getDopingRegionIdx(coord) ==
                   expectedIdx);
    EMCTEST_ASSERT_ISCLOSE(device.getDopingProfile().getDoping(coord),
                           expectedDoping, epsilon);
  }

  // TEST CONTACTS -----------------------------------------------
  // add contacts
  auto &surface = device.getSurface();
  EMCTEST_ASSERT(surface.getNrContacts() == 0);
  device.addOhmicContact(emcBoundaryPos::YMAX, -2, {0}, {1});
  device.addGateContact(emcBoundaryPos::YMAX, 4, {2}, {3}, 1, 1, 1);
  device.addOhmicContact(emcBoundaryPos::YMAX, 2, {4}, {5});
  device.addOhmicContact(emcBoundaryPos::YMIN, 0, {0}, {5});
  // device.addSchottkyContact(emcBoundaryPos::YMIN, 0, {0}, {5});

  // test: nr of contacts
  EMCTEST_ASSERT(surface.getNrContacts() == 4);

  // test: type of contact
  EMCTEST_ASSERT(surface.getContactType(0) == emcContactType::OHMIC);
  EMCTEST_ASSERT(surface.getContactType(1) == emcContactType::GATE);
  EMCTEST_ASSERT(surface.getContactType(2) == emcContactType::OHMIC);
  EMCTEST_ASSERT(surface.getContactType(3) == emcContactType::OHMIC);

  EMCTEST_ASSERT(
      surface.isContactType({3}, emcBoundaryPos::YMIN, emcContactType::OHMIC));
  EMCTEST_ASSERT(surface.isArtificialBoundary({3}, emcBoundaryPos::XMIN));
  EMCTEST_ASSERT(surface.isArtificialBoundary({3}, emcBoundaryPos::XMAX));
  EMCTEST_ASSERT(
      surface.isContactType({5}, emcBoundaryPos::YMAX, emcContactType::OHMIC));
  EMCTEST_ASSERT(
      surface.isContactType({2}, emcBoundaryPos::YMAX, emcContactType::GATE));
  EMCTEST_ASSERT(
      surface.isContactType({4}, emcBoundaryPos::YMAX, emcContactType::OHMIC));

  emcBoundaryPos boundPos;
  SizeVecSurface2D coordSurf;
  for (surface.initCoord(boundPos, coordSurf);
       !surface.isEndCoord(boundPos, coordSurf);
       surface.advanceCoord(boundPos, coordSurf)) {
    auto idxContact = surface.getContactIdx(coordSurf, boundPos);
    if (boundPos == emcBoundaryPos::YMIN) {
      EMCTEST_ASSERT(idxContact == 3);
      EMCTEST_ASSERT(surface.getContactVoltage(idxContact) == 0);
    } else if (boundPos == emcBoundaryPos::YMAX) {
      if (coordSurf[0] < 2) {
        EMCTEST_ASSERT(idxContact == 0);
        EMCTEST_ASSERT_ISCLOSE(surface.getContactVoltage(coordSurf, boundPos),
                               -2, epsilon);
        EMCTEST_ASSERT_ISCLOSE(surface.getContactVoltage(idxContact), -2,
                               epsilon);
      } else if (coordSurf[0] < 4) {
        EMCTEST_ASSERT(idxContact == 1);
        EMCTEST_ASSERT_ISCLOSE(surface.getContactVoltage(idxContact), 4,
                               epsilon);
        EMCTEST_ASSERT_ISCLOSE(surface.getContactVoltage(coordSurf, boundPos),
                               4, epsilon);
      } else {
        EMCTEST_ASSERT(idxContact == 2);
        EMCTEST_ASSERT_ISCLOSE(surface.getContactVoltage(idxContact), 2,
                               epsilon);
        EMCTEST_ASSERT_ISCLOSE(surface.getContactVoltage(coordSurf, boundPos),
                               2, epsilon);
      }
    } else if (boundPos == emcBoundaryPos::XMIN) {
      if (coordSurf[0] == 0) {
        EMCTEST_ASSERT(idxContact == 3);
      } else if (coordSurf[0] == 5) {
        EMCTEST_ASSERT(idxContact == 0);
      } else {
        EMCTEST_ASSERT(idxContact == -1);
      }
    } else {
      if (coordSurf[0] == 0) {
        EMCTEST_ASSERT(idxContact == 3);
      } else if (coordSurf[0] == 5) {
        EMCTEST_ASSERT(idxContact == 2);
      } else {
        EMCTEST_ASSERT(idxContact == -1);
      }
    }
  }

  if (printGrids) {
    std::cout << "Index of Doping Regions: \n";
    auto dopingIdx = device.getDopingProfile().getDopingRegionIdx();
    dopingIdx.print();
    std::cout << "Doping of Doping Regions: \n";
    auto doping = device.getDopingProfile().getDoping();
    doping.print();
    std::cout << "Index of Contact at each boundary: \n";
    surface.print();
  }
}