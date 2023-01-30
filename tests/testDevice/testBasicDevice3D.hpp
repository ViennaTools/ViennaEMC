#include <emcDevice.hpp>
#include <emcTestAsserts.hpp>

#include <iostream>

using NumType3D = double;
using MaterialType3D = emcMaterial<NumType3D>;
using DeviceType3D = emcDevice<NumType3D, 3>;
using ValueVec3D = DeviceType3D::ValueVec;
using SizeVec3D = DeviceType3D::SizeVec;
using SizeVecSurface3D = DeviceType3D::SizeVecSurface;
/*
Schematic of device that is tested:
z = 0
            <        OHMIC        >
              0   1   2   3   4   5
            ------------------------
          0 |          |           |
          1 | Region 0 | Region 1  |
          2 |          |           |
            ------------------------
          3 |          |           |
          4 | Region 2 | Region 3  |
          5 |          |           |
            ------------------------
              0   1   2   3   4   5
            <        GATE          >

z = 1
            <        OHMIC         >
              0   1   2   3   4   5
            ------------------------
          0 |                      |
          1 |                      |
          2 |      Region 4        |
          3 |                      |
          4 |                      |
          5 |                      |
            ------------------------
              0   1   2   3   4   5
            <        GATE          >

TESTS:
  - add regions
  - test grid of added regions
  - add contacts
  - test grid of added contacts
*/

void testBasicDevice3D(bool printGrids = false) {
  NumType3D epsilon = 1e-9;
  NumType3D electronMass = constants::me;

  // create silicon material (material characteristics + valleys)
  MaterialType3D silicon{11.8, 2329, 1.45e16, 9040, 1.14};

  // create device structure (material + extent + spacing)
  ValueVec3D maxPos = {5e-9, 10e-9, 5e-9};
  DeviceType3D device{silicon, maxPos, {1e-9, 2e-9, 5e-9}};

  // test: extent calculation
  EMCTEST_ASSERT(device.getGridExtent()[0] == 6);
  EMCTEST_ASSERT(device.getGridExtent()[1] == 6);
  EMCTEST_ASSERT(device.getGridExtent()[2] == 2);

  // TEST DOPINGREGIONS -------------------------------------------
  // add different doping regions
  device.addConstantDopingRegion({0, 0, 0}, {maxPos[0] / 2, maxPos[1] / 2, 0},
                                 1e6);
  device.addConstantDopingRegion({0, maxPos[1] / 2, 0},
                                 {maxPos[0] / 2, maxPos[1], 0}, 2e6);
  device.addConstantDopingRegion({maxPos[0] / 2, 0, 0},
                                 {maxPos[0], maxPos[1] / 2, 0}, 3e6);
  device.addConstantDopingRegion({maxPos[0] / 2, maxPos[1] / 2, 0},
                                 {maxPos[0], maxPos[1], 0}, 4e6);
  device.addConstantDopingRegion({0, 0, maxPos[2]},
                                 {maxPos[0], maxPos[1], maxPos[2]}, 5e6);

  // test: Doping Regions
  SizeVec3D coord;
  SizeType expectedIdx;
  NumType3D expectedDoping;
  for (coord.fill(0); !device.isEndCoord(coord); device.advanceCoord(coord)) {
    if (coord[2] == 0) {
      if (coord[0] < 3 && coord[1] < 3) {
        expectedIdx = 0;
        expectedDoping = 1e06;
      } else if (coord[0] < 3 && coord[1] >= 3) {
        expectedIdx = 1;
        expectedDoping = 2e06;
      } else if (coord[0] >= 3 && coord[1] < 3) {
        expectedIdx = 2;
        expectedDoping = 3e06;
      } else {
        expectedIdx = 3;
        expectedDoping = 4e6;
      }
    } else {
      expectedIdx = 4;
      expectedDoping = 5e6;
    }
    EMCTEST_ASSERT(device.getDopingProfile().getDopingRegionIdx(coord) ==
                   expectedIdx);
    EMCTEST_ASSERT_ISCLOSE(device.getDopingProfile().getDoping(coord),
                           expectedDoping, epsilon);
  }

  // TEST CONTACTS -----------------------------------------------
  // add contacts
  EMCTEST_ASSERT(device.getSurface().getNrContacts() == 0);
  device.addGateContact(emcBoundaryPos::YMIN, 0, {0, 0}, {maxPos[0], maxPos[2]},
                        1, 1, 1);
  device.addOhmicContact(emcBoundaryPos::YMAX, 4, {0, 0},
                         {maxPos[0], maxPos[2]});

  // test: nr of contacts
  EMCTEST_ASSERT(device.getSurface().getNrContacts() == 2);

  // test: type of contact
  auto &surface = device.getSurface();
  EMCTEST_ASSERT(surface.getContactType(0) == emcContactType::GATE);
  EMCTEST_ASSERT(surface.getContactType(1) == emcContactType::OHMIC);
  EMCTEST_ASSERT(surface.isArtificialBoundary({3, 1}, emcBoundaryPos::XMIN));
  EMCTEST_ASSERT(surface.isContactType({3, 0}, emcBoundaryPos::YMIN,
                                       emcContactType::GATE));
  EMCTEST_ASSERT(surface.isArtificialBoundary({3, 0}, emcBoundaryPos::ZMIN));
  EMCTEST_ASSERT(surface.isContactType({3, 0}, emcBoundaryPos::YMAX,
                                       emcContactType::OHMIC));

  emcBoundaryPos boundPos;
  SizeVecSurface3D coordSurf;
  for (surface.initCoord(boundPos, coordSurf);
       !surface.isEndCoord(boundPos, coordSurf);
       surface.advanceCoord(boundPos, coordSurf)) {
    auto idxContact = surface.getContactIdx(coordSurf, boundPos);
    if (boundPos == emcBoundaryPos::YMIN) {
      EMCTEST_ASSERT(idxContact == 0);
      EMCTEST_ASSERT_ISCLOSE(surface.getContactVoltage(idxContact), 0, epsilon);
      EMCTEST_ASSERT_ISCLOSE(surface.getContactVoltage(coordSurf, boundPos), 0,
                             epsilon);
    } else if (boundPos == emcBoundaryPos::YMAX ||
               (boundPos == emcBoundaryPos::XMIN && coordSurf[0] == 5) ||
               (boundPos == emcBoundaryPos::XMAX && coordSurf[0] == 5) ||
               (boundPos == emcBoundaryPos::ZMIN && coordSurf[1] == 5) ||
               (boundPos == emcBoundaryPos::ZMAX && coordSurf[1] == 5)) {
      EMCTEST_ASSERT(idxContact == 1);
      EMCTEST_ASSERT_ISCLOSE(surface.getContactVoltage(idxContact), 4, epsilon);
      EMCTEST_ASSERT_ISCLOSE(surface.getContactVoltage(coordSurf, boundPos), 4,
                             epsilon);
    } else {
      EMCTEST_ASSERT(idxContact == -1);
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