cmake_minimum_required(VERSION 3.5)
                                                                                                                  

INCLUDE_DIRECTORIES(include/)

# ---------------------------------------------------------
# Compile the simulator sources into a library
#
SET(SIMULATOR_SOURCES
   src/acoustic.c
   src/apply_voltage.c
   src/check_source_drain.c
   src/core.c
   src/coulomb_scattering.c
   src/relaxation_scattering.c
   src/device_structure.c
   src/drift.c
   src/electric_field_update.c
   src/electrons_initialization.c
   src/free_flight_scatter.c
   src/init_kspace.c
   src/integrate.c
   src/intervalley.c
   src/mat_par_initialization.c
   src/NEC_scheme.c
   src/output.c
   src/poisson_SOR.c
   src/scatter_carrier.c
   src/scattering_table.c
   src/velocity_energy_cumulative.c
)

# ---------------------------------------------------------
#
SET(TARGET cemcsim)
ADD_EXECUTABLE(${TARGET} src/main.c ${SIMULATOR_SOURCES})
TARGET_LINK_LIBRARIES(${TARGET} m -fopenmp)



# ---------------------------------------------------------
# Build the tests
#
#ADD_SUBDIRECTORY(tests)
