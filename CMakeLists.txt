cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

project("ViennaEMC"   VERSION 1.0.0
                      DESCRIPTION "Multi-Valley Ensemble Monte Carlo Code for Device Simulation.")

# provides plattform independent installation directories
include(GNUInstallDirs)

# set cmake standard to c++17 (to enable structured bindings)
set(CMAKE_CXX_STANDARD "17")

# set default build type
SET(DEFAULT_BUILD_TYPE "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${DEFAULT_BUILD_TYPE}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${DEFAULT_BUILD_TYPE}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel"
                                                        "RelWithDebInfo")
endif()

set(LOCAL_CONFIG_DIR "lib/cmake/${PROJECT_NAME}")
set(VIENNAEMC_BUILD_INCLUDE_DIRS "${${PROJECT_NAME}_SOURCE_DIR}/include")
set(VIENNAEMC_INCLUDE_DIRS "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}")
set(VIENNAEMC_BINARY_DIR ${PROJECT_BINARY_DIR} CACHE PATH "Path to local ViennaEMC installation" FORCE)
set(VIENNAEMC_CONFIG_DIR "${CMAKE_INSTALL_PREFIX}/${LOCAL_CONFIG_DIR}")

# ##################################################################################################
# DEPENDENCIES
# ##################################################################################################
set(DEPENDENCIES_DIR ${CMAKE_SOURCE_DIR}/dependencies)

# add stamp-file (which enforces rerun of cmake every time buildDependencies is executed)
# found at: https://github.com/yozoon/ViennaPS-LibTorch-Example/blob/main/CMakeLists.txt
set(STAMPFILE ${CMAKE_BINARY_DIR}/target.stamp)
# Create the initial stamp file
file(TOUCH ${STAMPFILE})
# Include the stamp file, so that cmake is forced to re-run once the file has been touched
include(${STAMPFILE})

include(ExternalProject)
add_custom_target(buildDependencies COMMAND ${CMAKE_COMMAND} -E touch ${STAMPFILE})
add_subdirectory(external/upstream)

find_package(OpenMP)
if(OpenMP_FOUND)
  message(STATUS "OpenMP dependency found.")
  list(APPEND VIENNAEMC_LIBRARIES OpenMP::OpenMP_CXX)
else()
  message(STATUS "OpenMP dependency not found.")
endif()

# ##################################################################################################
# CONFIGURATION LIBRARY
# ##################################################################################################
add_library(${PROJECT_NAME} INTERFACE)

target_include_directories(
    ${PROJECT_NAME}
    INTERFACE   $<BUILD_INTERFACE:${VIENNAEMC_BUILD_INCLUDE_DIRS}>
                $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

# add expected compiler feature (C++17-standard)
target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_17)

include(CMakePackageConfigHelpers)
# create package version file
write_basic_package_version_file("${PROJECT_NAME}ConfigVersion.cmake"
                                    VERSION ${PROJECT_VERSION}
                                    COMPATIBILITY AnyNewerVersion)
# create config file for project
configure_package_config_file(
    "${PROJECT_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake.in"
    "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
    INSTALL_DESTINATION ${LOCAL_CONFIG_DIR})

# ##################################################################################################
# BUILD EXAMPLES
# ##################################################################################################
option(VIENNAEMC_BUILD_EXAMPLES "Build examples." OFF)
if(VIENNAEMC_BUILD_EXAMPLES)
  add_subdirectory(examples)
endif(VIENNAEMC_BUILD_EXAMPLES)

# ##################################################################################################
# BUILD TESTS
# ##################################################################################################
option(VIENNAEMC_BUILD_TESTS "Build tests." OFF)
if(VIENNAEMC_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif(VIENNAEMC_BUILD_TESTS)

# ##################################################################################################
# INSTALL
# ##################################################################################################
install(TARGETS ${PROJECT_NAME}
        EXPORT ${PROJECT_NAME}_Targets
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})


# file containing target is created + copied to install folder
install(EXPORT ${PROJECT_NAME}_Targets
    FILE ${PROJECT_NAME}Targets.cmake
    NAMESPACE ${PROJECT_NAME}::
    DESTINATION ${LOCAL_CONFIG_DIR})

# cmake config- und versionfile copied to install folder
install(FILES   "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
                "${PROJECT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
        DESTINATION ${LOCAL_CONFIG_DIR})

# copy header files in installation folder
install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/ DESTINATION include/${PROJECT_NAME})

# ##################################################################################################
# UNINSTALL TARGET
# ##################################################################################################
if(NOT TARGET uninstall)
  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in"
                 "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake" IMMEDIATE @ONLY)

  add_custom_target(uninstall COMMAND ${CMAKE_COMMAND} -P
                                      ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
endif()
