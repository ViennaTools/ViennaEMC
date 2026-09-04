cmake_minimum_required(VERSION 3.4)

project("fullBandResistor")

# Full-band examples read .matpkg (HDF5) packages: they need the HDF5 C
# library in addition to the ViennaEMC dependencies. The examples build used
# to fail for the whole tree because this directory had no CMakeLists.txt
# while examples/CMakeLists.txt adds every subdirectory.
find_package(HDF5 REQUIRED COMPONENTS C)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")

add_executable(${PROJECT_NAME} ${PROJECT_NAME}.cpp)
target_include_directories(${PROJECT_NAME} PUBLIC ${VIENNAEMC_INCLUDE_DIRS} ${HDF5_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} PRIVATE ${VIENNAEMC_LIBRARIES} ${HDF5_C_LIBRARIES})

# Add to example target
add_dependencies(buildExamples ${PROJECT_NAME})
