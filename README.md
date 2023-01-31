# ViennaEMC

ViennaEMC is a Semiconductor Device Simulation library based on the Multi-Valley Ensemble Monte Carlo Method. 

The aim of this project is to be able to simulate different kind of devices with different materials and also different particle types (e.g. electrons). Also it should be possible to include the calculation of the particle-particle interactions in the simulation, this is done with the Fast Multipole Method (FMM) and the library [ScalFMM](https://gitlab.inria.fr/solverstack/ScalFMM).

Finally, also parts of the library can be reused for the simulation of bulk-material with an applied electric background field, as shown in the examples.

## Building

### Supported Operating Systems

- Linux (g++ / clang)

### System Requirements

- C++17 compiler with OpenMP support
- CMake (Version 3.12 or newer)

### Dependencies (installed automatically)

- [ScalFMM](https://gitlab.inria.fr/solverstack/ScalFMM): C++ library that performs N-body simulation using kernel independent Fast Multipole Method

## Installing

Since this is a header only project, it does not require any installation. However, we recommend the following procedure in order to set up all dependencies correctly.

### Installing Dependencies and ViennaEMC

The following steps are recommended if you want to install the library and the dependencies:

    git clone github.com/ViennaTools/ViennaEMC.git
    cd ViennaEMC
    mkdir build && cd build
    cmake .. -DCMAKE_INSTALL_PREFIX=/path/to/your/custom/install/
    make buildDependencies
    make install

This will install the necessary headers and CMake files to the specified path. If DCMAKE_INSTALL_PREFIX is not specified, it will be installed to the standard path for your system, usually /usr/local/ .

### Installing ViennaEMC with pre-installed dependencies

If you want to use your own install of scalFMM, just specify the directory in CMake:

    git clone github.com/ViennaTools/ViennaEMC.git
    cd ViennaEMC
    mkdir build && cd build
    cmake .. -DCMAKE_INSTALL_PREFIX=/path/to/your/custom/install/ -Dscalfmm_DIR=/path/to/scalfmm/install
    make install

## Integration in CMake projects

In order to use this library in your CMake project, add the following lines to the CMakeLists.txt of your project:

    set(ViennaEMC_DIR "/path/to/your/custom/install/lib/cmake/ViennaEMC/")
    find_package(ViennaEMC REQUIRED)
    add_executable(...)
    target_include_directories(${PROJECT_NAME} PUBLIC ${VIENNAEMC_INCLUDE_DIRS})
    target_link_libraries(${PROJECT_NAME} ${VIENNAEMC_LIBRARIES})

## Building examples

The examples can be built using cmake

    mkdir build && cd build
    cmake .. -DVIENNAEMC_BUILD_EXAMPLES=ON
    make

The examples can be divided into device simulations and bulk simulations. While the classes of the library are written for device simulation, for bulk simulations some custom classes were written.

The examples that perform bulk simulations are:

- [Bulk Simulation of Silicon:](examples/bulkSimulation/) This example simulates electron transport in bulk silicon. The band structure of silicon is approximated within this simulation with the analytical band approach, considering only the minima close to the X-points of the first Brillouin zone. The considered scatter mechanisms include acoustic and optical phonon scattering.
- [Bulk Simulation of Silicon with the usage of real-space particle-particle interactions:](examples/bulkSimulationFMM/) This example simulates electron transport in bulk silicon as described for the previous example. Additionally, real-space particle-particle interactions with the help of [ScalFMM](https://gitlab.inria.fr/solverstack/ScalFMM) are included. Moreover, also donors are considered for the particle-particle interactions. 
- [Bulk Simulation of a Single Layer of Molybdenum Disulfide (MoS2):](examples/singleLayerMoS2/) This example simulates the characteristics of a monolayer of Mos2, in which the transport only happens in two dimensions. The simulation can be based on different papers, depending on the paper different valleys and characteristics are used. The applied scatter mechanisms are again acoustic and optical phonons.

The given examples that perform device simulations are:

- [2D MOSFET Simulation:](examples/mosfet2D/) This examples simulates a n-channel MOSFET built with silicon and is used to compare the results of this library with the results of the classical ensemble Monte Carlo simulator (CEMC) of [ViennaWD](https://www.iue.tuwien.ac.at/software/viennawd/). The device structure is described in more details within the example.
- [2D Resistor simulation:](examples/resistor2D/) This example simulates a resistor built from silicon. The device structure is described in more details within the example.

More information on each example can be found in the corresponding folder.

## Running Tests

The tests can be built using cmake:

    mkdir build && cd build
    cmake .. -DVIENNAEMC_BUILD_TESTS=ON
    make buildTests
    make test

## Documentation

The documentation is done with [Doxygen](https://www.doxygen.nl/index.html) and can be built locally:

    cd docs/doxygen
    ./make_doxygen.sh
    cd html/
    firefox index.html


## Plotting the Results

The results can be plotted using the Python - Files in the folder `helper/emcPlottingFiles`. `emcPlottingFiles` can be installed with pip by opening a terminal in the folder `helper` and using: `pip install emcPlottingFiles`, needed python libraries for that installation are: `numpy`, `matplotlib`, `pandas` and `scipy`.

For each given example in `examples` a plotting script is
provided that shows how to use the python files for plotting the current results. **Note:** to use these python scripts `emcPlottingFiles` has to be installed.

## Authors
List of former contributors: Laura Gollner, Robin Steiner and Anna Benzer

This library was developed under the aegis of the [Institute for Microelectronics](http://www.iue.tuwien.ac.at/) at [TU Wien](https://www.tuwien.at/).

## License

See file LICENSE in the base directory.