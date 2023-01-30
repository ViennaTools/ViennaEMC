#include <array>

#include <emcGrid.hpp>
#include <emcTestAsserts.hpp>

using NumType = double;

int main() {
  // TEST DIMENSION 3
  emcGrid<NumType, 3> grid3D{{4, 4, 3}, 1};
  // alternative:
  // emcGrid<NumType, 3> grid3D{4,4,3};
  // grid3D.fill(1);

  // // test: size of grid
  EMCTEST_ASSERT(grid3D.getSize() == 48);
  EMCTEST_ASSERT(grid3D.getSize(0) == 4);
  EMCTEST_ASSERT(grid3D.getSize(1) == 4);
  EMCTEST_ASSERT(grid3D.getSize(2) == 3);

  // test: fill-function
  SizeType counter = 0;
  for (auto itGrid = grid3D.begin(); itGrid < grid3D.end(); itGrid++) {
    if (*itGrid == 1)
      counter++;
  }
  EMCTEST_ASSERT(counter == 48);

  // test: element access
  std::array<SizeType, 3> currCoord;
  grid3D.iota(0);
  counter = 0;
  for (SizeType i = 0; i < grid3D.getSize(2); i++) {
    for (SizeType j = 0; j < grid3D.getSize(1); j++) {
      for (SizeType k = 0; k < grid3D.getSize(0); k++) {
        currCoord = {{k, j, i}};
        EMCTEST_ASSERT(counter == grid3D[currCoord]);
        counter++;
      }
    }
  }
  EMCTEST_ASSERT(counter == 48);

  // test: access neighbour values
  std::array<SizeType, 3> tmpCoord;
  for (SizeType i = 0; i < grid3D.getSize(2); i++) {
    for (SizeType j = 0; j < grid3D.getSize(1); j++) {
      for (SizeType k = 0; k < grid3D.getSize(0); k++) {
        currCoord = {{k, j, i}};
        if (!grid3D.onBoundary(currCoord)) {
          for (SizeType idxDim = 0; idxDim < 3; idxDim++) {
            tmpCoord = currCoord;
            tmpCoord[idxDim]++;
            EMCTEST_ASSERT(grid3D.getNextValue(currCoord, idxDim) ==
                           grid3D[tmpCoord]);
            tmpCoord[idxDim] -= 2;
            EMCTEST_ASSERT(grid3D.getPrevValue(currCoord, idxDim) ==
                           grid3D[tmpCoord]);
          }
        }
      }
    }
  }

  // test iteration over data-structure with coordinates
  counter = 0;
  for (currCoord.fill(0); !grid3D.isEndCoord(currCoord);
       grid3D.advanceCoord(currCoord)) {
    EMCTEST_ASSERT(counter == grid3D[currCoord]);
    counter++;
  }
  EMCTEST_ASSERT(counter == 48);

  counter = 0;
  currCoord.fill(0);
  for (auto component : grid3D) {
    EMCTEST_ASSERT(counter == grid3D[currCoord]);
    grid3D.advanceCoord(currCoord);
    counter++;
  }
  EMCTEST_ASSERT(counter == 48);

  // test: partial fill (fill cube with 8 entries)
  grid3D.fill(-1, {0, 0, 0}, {1, 1, 1});
  counter = 0;
  // grid3D.print();
  for (currCoord.fill(0); !grid3D.isEndCoord(currCoord);
       grid3D.advanceCoord(currCoord)) {
    if (grid3D[currCoord] == -1)
      counter++;
  }
  EMCTEST_ASSERT(counter == 8);

  // ---------------------------------------------------------------------------
  // TEST DIMENSION 2
  emcGrid<NumType, 2> grid2D{{2, 4}};

  // test: size of grid
  EMCTEST_ASSERT(grid2D.getSize() == 8);
  EMCTEST_ASSERT(grid2D.getSize(0) == 2);
  EMCTEST_ASSERT(grid2D.getSize(1) == 4);
  grid2D.iota(1);
  // grid2D.print();

  // test: access neighbour values
  std::array<SizeType, 2> currCoord2D, tmpCoord2D;
  for (SizeType j = 0; j < grid2D.getSize(1); j++) {
    for (SizeType k = 0; k < grid2D.getSize(0); k++) {
      currCoord2D = {{k, j}};
      if (!grid2D.onBoundary(currCoord2D)) {
        for (SizeType idxDim = 0; idxDim < 2; idxDim++) {
          tmpCoord2D = currCoord2D;
          tmpCoord2D[idxDim]++;
          EMCTEST_ASSERT(grid2D.getNextValue(currCoord2D, idxDim) ==
                         grid2D[tmpCoord2D]);
          tmpCoord[idxDim] -= 2;
          EMCTEST_ASSERT(grid2D.getPrevValue(currCoord2D, idxDim) ==
                         grid2D[tmpCoord2D]);
        }
      }
    }
  }

  // test: copy construction
  emcGrid<NumType, 2> grid2Dcopy = grid2D;
  std::array<SizeType, 2> coord2D;
  counter = 0;
  EMCTEST_ASSERT(grid2D.getSize() == grid2Dcopy.getSize());
  for (coord2D.fill(0); !grid2D.isEndCoord(coord2D);
       grid2D.advanceCoord(coord2D)) {
    EMCTEST_ASSERT(grid2D[coord2D] == grid2Dcopy[coord2D]);
    counter++;
  }
  EMCTEST_ASSERT(counter == 8);

  // test: deep copy (copy constr.) + operator
  grid2D += grid2D;
  for (coord2D.fill(0); !grid2D.isEndCoord(coord2D);
       grid2D.advanceCoord(coord2D)) {
    EMCTEST_ASSERT(grid2D[coord2D] == 2 * grid2Dcopy[coord2D]);
  }

  // test: copy assignment
  grid2D = grid2Dcopy;
  for (coord2D.fill(0); !grid2D.isEndCoord(coord2D);
       grid2D.advanceCoord(coord2D)) {
    EMCTEST_ASSERT(grid2D[coord2D] == grid2Dcopy[coord2D]);
  }

  // test: deep copy (copy assignment) + operator
  grid2D += grid2D;
  for (coord2D.fill(0); !grid2D.isEndCoord(coord2D);
       grid2D.advanceCoord(coord2D)) {
    EMCTEST_ASSERT(grid2D[coord2D] == 2 * grid2Dcopy[coord2D]);
  }

  // ---------------------------------------------------------------------------
  // TEST DIMENSION 1
  emcGrid<NumType, 1> grid1D{{4}};
  grid1D.fill(2);

  // test: extent
  EMCTEST_ASSERT(grid1D.getSize() == 4);
  EMCTEST_ASSERT(grid1D.getSize(0) == 4);

  // test: operators
  grid1D += grid1D;
  grid1D = grid1D - grid1D;
  for (const auto &component : grid1D)
    EMCTEST_ASSERT(component == 0);

  emcGrid<NumType, 1> grid1Dtest{{4}};
  grid1Dtest.fill(5);
  grid1D = grid1D + grid1Dtest;
  grid1Dtest = grid1Dtest + grid1Dtest;
  for (const auto &component : grid1D)
    EMCTEST_ASSERT(component == 5);
  for (const auto &component : grid1Dtest)
    EMCTEST_ASSERT(component == 10);

  return 0;
}
