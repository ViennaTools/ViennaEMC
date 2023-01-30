cmake_minimum_required(VERSION 3.4)

project(VIENNAEMCTests)

include("${CMAKE_SOURCE_DIR}/cmake/prepare.cmake")

find_package(ViennaEMC CONFIG PATHS ${VIENNAEMC_BINARY_DIR})
if(NOT ViennaEMC_FOUND)
  message(
    FATAL_ERROR
    "One of the dependencies hasn't been built yet and thus can't be used by the buildTests target.\n\
    Try the following steps to fix the problem:\n\
    1. Turn off tests by disabling the VIENNAEMC_BUILD_TESTS option when configuring the project\n\
    (e.g. with the `-D VIENNAEMC_BUILD_TESTS=OFF` flag)\n\
    2. Build the dependencies target: e.g. `make buildDependencies` when using Makefiles\n\
    3. Turn on tests again by re-enabling the VIENNAEMC_BUILD_TESTS option\n\
    (e.g. with the `-D VIENNAEMC_BUILD_TESTS=ON` flag)\n\
    4. Now the configure step should pick up the built dependencies and proceed as intended.")
endif()

list(PREPEND VIENNAEMC_INCLUDE_DIRS ${VIENNAEMC_BUILD_INCLUDE_DIRS})

add_custom_target(buildTests)

SUBDIRLIST(TEST_DIRS ${CMAKE_CURRENT_LIST_DIR})

# add all examples
foreach(subdir ${TEST_DIRS})
  add_subdirectory(${subdir})
endforeach()