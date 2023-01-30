#include <ValleyTypes/emcNonParabolicAnistropValley.hpp>
#include <emcConstants.hpp>
#include <emcTestAsserts.hpp>

using NumType = double;
using ValleyType = emcNonParabolicAnisotropValley<NumType>;

static const NumType eps = 1e-9;

void testSubValley(ValleyType &valley, SizeType idxSubValley,
                   std::array<NumType, 3> dir1, std::array<NumType, 3> dir2,
                   std::array<NumType, 3> dir3) {
  std::array<NumType, 3> testVec, testVecEllipse, testVecDevice;
  testVec = {{1, 2, 3}};
  testVecEllipse = valley.transformToEllipseCoord(idxSubValley, testVec);
  testVecDevice = valley.transformToDeviceCoord(idxSubValley, testVecEllipse);

  // ellipse vector has to be rightly transformed to the ECS
  normalize(dir1);
  normalize(dir2);
  normalize(dir3);
  EMCTEST_ASSERT_ISCLOSE(testVecEllipse[0], innerProduct(testVec, dir1), eps);
  EMCTEST_ASSERT_ISCLOSE(testVecEllipse[1], innerProduct(testVec, dir2), eps);
  EMCTEST_ASSERT_ISCLOSE(testVecEllipse[2], innerProduct(testVec, dir3), eps);

  // original + resulting device vector have to be the same
  EMCTEST_ASSERT_ISCLOSE(testVec[0], testVecDevice[0], eps);
  EMCTEST_ASSERT_ISCLOSE(testVec[1], testVecDevice[1], eps);
  EMCTEST_ASSERT_ISCLOSE(testVec[2], testVecDevice[2], eps);
}

int main() {
  NumType electronMass = constants::me;
  NumType eps = 1e-9;

  // test: random X-valley
  // 6 subvalleys with longitudinal mass in following directions:
  // (1,0,0), (-1,0,0), (0,1,0), (0,-1,0), (0,0,1), (0,0,-1)
  NumType mlX = 0.9;
  NumType mtX = 0.1;
  ValleyType xValley({mlX, mtX, mtX}, electronMass, 6, 0.5);

  // dir1 + dir2 + dir3 have to be an orthogonal basis!
  xValley.setSubValleyEllipseCoordSystem(0, {1, 0, 0}, {0, 1, 0}, {0, 0, 1});
  xValley.setSubValleyEllipseCoordSystem(1, {-1, 0, 0}, {0, 1, 0}, {0, 0, 1});
  xValley.setSubValleyEllipseCoordSystem(2, {0, 1, 0}, {1, 0, 0}, {0, 0, 1});
  xValley.setSubValleyEllipseCoordSystem(3, {0, -1, 0}, {1, 0, 0}, {0, 0, 1});
  xValley.setSubValleyEllipseCoordSystem(4, {0, 0, 1}, {0, 1, 0}, {1, 0, 0});
  xValley.setSubValleyEllipseCoordSystem(5, {0, 0, -1}, {0, 1, 0}, {1, 0, 0});

  testSubValley(xValley, 0, {1, 0, 0}, {0, 1, 0}, {0, 0, 1});
  testSubValley(xValley, 1, {-1, 0, 0}, {0, 1, 0}, {0, 0, 1});
  testSubValley(xValley, 2, {0, 1, 0}, {1, 0, 0}, {0, 0, 1});
  testSubValley(xValley, 3, {0, -1, 0}, {1, 0, 0}, {0, 0, 1});
  testSubValley(xValley, 4, {0, 0, 1}, {0, 1, 0}, {1, 0, 0});
  testSubValley(xValley, 5, {0, 0, -1}, {0, 1, 0}, {1, 0, 0});

  // test: random L-valley
  // 8 subvalleys with longitudinal mass in following directions:
  // longitudinal mass is in (1,1,1), (-1,-1,-1), (-1,1,1), (1,-1,1), (1,1,-1),
  // (-1,-1,1), (-1,1,-1), (1,-1,-1),
  NumType mlL = 0.8;
  NumType mtL = 0.2;
  ValleyType lValley({mlL, mtL, mtL}, electronMass, 8, 0.3);

  // dir1 + dir2 + dir3 have to be an orthogonal basis!
  lValley.setSubValleyEllipseCoordSystem(0, {1, 1, 1}, {-1, 1, 0}, {-1, -1, 2});
  lValley.setSubValleyEllipseCoordSystem(1, {-1, -1, -1}, {-1, 1, 0},
                                         {-1, -1, 2});
  lValley.setSubValleyEllipseCoordSystem(2, {-1, 1, 1}, {1, 1, 0}, {1, -1, 2});
  lValley.setSubValleyEllipseCoordSystem(3, {1, -1, 1}, {1, 1, 0}, {-1, 1, 2});
  lValley.setSubValleyEllipseCoordSystem(4, {1, 1, -1}, {1, 0, 1}, {-1, 2, 1});
  lValley.setSubValleyEllipseCoordSystem(5, {-1, -1, 1}, {1, 0, 1}, {-1, 2, 1});
  lValley.setSubValleyEllipseCoordSystem(6, {-1, 1, -1}, {1, 1, 0}, {-1, 1, 2});
  lValley.setSubValleyEllipseCoordSystem(7, {-1, -1, 1}, {1, 0, 1}, {-1, 2, 1});

  testSubValley(lValley, 0, {1, 1, 1}, {-1, 1, 0}, {-1, -1, 2});
  testSubValley(lValley, 1, {-1, -1, -1}, {-1, 1, 0}, {-1, -1, 2});
  testSubValley(lValley, 2, {-1, 1, 1}, {1, 1, 0}, {1, -1, 2});
  testSubValley(lValley, 3, {1, -1, 1}, {1, 1, 0}, {-1, 1, 2});
  testSubValley(lValley, 4, {1, 1, -1}, {1, 0, 1}, {-1, 2, 1});
  testSubValley(lValley, 5, {-1, -1, 1}, {1, 0, 1}, {-1, 2, 1});
  testSubValley(lValley, 6, {-1, 1, -1}, {1, 1, 0}, {-1, 1, 2});
  testSubValley(lValley, 7, {-1, -1, 1}, {1, 0, 1}, {-1, 2, 1});

  return 0;
}