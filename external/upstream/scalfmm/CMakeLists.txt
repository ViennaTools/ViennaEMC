
# ##################################################################################################
# Check scalfmm dependency
# ##################################################################################################

# Check if the variable is a cache variable
get_property(
  is_cached
  CACHE scalfmm_DIR
  PROPERTY TYPE)

# If it is not from cache and defined at this point, the variable was passed as a -D command line
# argument
if("${is_cached}" STREQUAL "" AND DEFINED scalfmm_DIR)
    find_package(scalfmm REQUIRED PATHS ${scalfmm_DIR} NO_DEFAULT_PATH)
    message(STATUS "Found scalfmm at ${scalfmm_DIR}.")
    add_library(scalfmm_external INTERFACE) # dummy target
else()
    # First try to locate any local build
    find_package(scalfmm QUIET PATHS "${DEPENDENCIES_DIR}/Install/scalfmm_external" NO_DEFAULT_PATH)

    # If that fails look for a system installation
    if(NOT scalfmm_FOUND)
        find_package(scalfmm QUIET)
    endif()

    if(scalfmm_FOUND)
        message(STATUS "Found scalfmm at ${scalfmm_DIR}.")
        add_library(scalfmm_external INTERFACE) # dummy target
    else(scalfmm_FOUND)
        message(
            STATUS
            "Couldn't find scalfmm. Dependency will be built with the buildDependencies target (e.g. `make buildDependencies`).")

        # install scalFMM in dependencies folder
        ExternalProject_Add(
            scalfmm_external
            GIT_REPOSITORY https://gitlab.inria.fr/solverstack/ScalFMM.git
            GIT_TAG master
            GIT_SHALLOW true
            INSTALL_DIR "${DEPENDENCIES_DIR}/Install/scalfmm_external"
            STAMP_DIR "${DEPENDENCIES_DIR}/Stamp/scalfmm_external"
            GIT_SUBMODULES  "CMakeModules/morse_cmake"
                            "inastemp"
            CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
                    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                    -DCMAKE_CXX_STANDARD=${CMAKE_CXX_STANDARD}
                    -DCMAKE_CXX_EXTENSIONS=${CMAKE_CXX_EXTENSIONS}
                    -DCMAKE_CXX_STANDARD_REQUIRED=${CMAKE_CXX_STANDARD_REQUIRED}
                    -DCMAKE_BUILD_TYPE=Release
                    -DSCALFMM_BUILD_EXAMPLES=OFF
            BUILD_COMMAND ${CMAKE_COMMAND} --build .
            USES_TERMINAL_DOWNLOAD 1
            USES_TERMINAL_UPDATE 1
            GIT_PROGRESS 1
            BUILD_ALWAYS 0)

        ExternalProject_Get_Property(scalfmm_external INSTALL_DIR)
        set(scalfmm_DIR "${INSTALL_DIR}")

        add_dependencies(buildDependencies scalfmm_external)
    endif()
    set(scalfmm_DIR "${scalfmm_DIR}" CACHE PATH "Path to scalfmm installation" FORCE)
endif()
