/*
=============================================================================
   Copyright (c) 2013, Institute for Microelectronics, TU Wien
   http://www.iue.tuwien.ac.at
                             -----------------
              ViennaWD - The Vienna Wigner Decoherence Algorithms
                         Ensemble Monte Carlo Simulator
                             -----------------

   authors:    Marek Pobjecky                          
               Mihail Nedjalkov                  nedjalkov@iue.tuwien.ac.at

   license:    see file LICENSE in the base directory
============================================================================= 
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "emc.h"
#include <omp.h>

/*=== main program ===*/
int main(int argc, char *argv[]) {
    const_t constants;
    geometry_t *geometry;
    scatpar_t *scatpar;
    phys_quant_t *phys_quantities;
    output_t *outputdata;
    int threads;

    FILE *input1, *input2;

    if (argc != 3) {
        printf("Error: Wrong execution - Usage: %s device_dimensions.in  simulation_parameters.in\n", argv[0]);
        printf("Aborting ...\n");
        return -1;
    }

    /*=== allocate memory for data ===*/
    geometry = (geometry_t *) calloc(1, sizeof(geometry_t));
    scatpar = (scatpar_t *) calloc(1, sizeof(scatpar_t));
    phys_quantities = (phys_quant_t *) calloc(1, sizeof(phys_quant_t));

    /*        INITIALIZATION        */
    /* ============================ */
    printf("===========================================================\n");
    /*=== initialize material parameters ===*/
    oooMatParInitialization(&constants);
    /*=== initialize scattering parameters ===*/
    oooScatParInitialization(scatpar);
    printf("===========================================================\n");
    /*=== read geometry data and doping concentrations from file ===*/
    input1 = fopen(argv[1], "r");
    if (input1 != NULL){		                                                                       /* default:   */
        if (fscanf(input1, "%*s %lf", &geometry->deviceLength) == 0) geometry->deviceLength = 600e-9;  /*  600 nm    */
        if (fscanf(input1, "%*s %lf", &geometry->deviceHeight) == 0) geometry->deviceHeight = 100e-9;  /*  100 nm    */
        if (fscanf(input1, "%*s %lf", &geometry->deviceWidth) == 0) geometry->deviceWidth = 1e-6;      /*  1 um      */
        if (fscanf(input1, "%*s %lf", &geometry->meshSize) == 0) geometry->meshSize = 1e-9;            /*  1 nm      */
        if (fscanf(input1, "%*s %d", &scatpar->kTLevels) == 0) scatpar->kTLevels = 2;                  /*  2         */
        if (fscanf(input1, "%*s %d", &scatpar->precision) == 0) scatpar->precision = 10; /* 10 steps for e- energies */
        if (fscanf(input1, "%*s %d", &threads) == 0) threads = 2;
        scatpar->dopingCon = oooIntegrate(constants, scatpar->kTLevels);                               /*  5e25 m^-3 */
        constants.Ni = 1.0 * scatpar->dopingCon;
        constants.debyeLength /= sqrt(constants.Ni);
        omp_set_num_threads(threads);
        fclose(input1);
    } else {
        printf("ERROR: File '%s' not found!\n", argv[1]);
        free(geometry);
        free(scatpar);
        free(phys_quantities);
        exit(1);
    }
    printf("Doping read from file '%s':\n", argv[1]);
    printf("Doping = %.5e m^-3\n", scatpar->dopingCon);
    printf("Threads = %d\n", threads);
    printf("===========================================================\n");
    /*=== initialize device structure ===*/
    oooDeviceStructureInitialization(constants, geometry, scatpar);
    printf("===========================================================\n");
    /*=== read simulation times and voltages from file ===*/
    input2 = fopen(argv[2], "r");
    if (input2 != NULL)
    {                                                                                                   /* default:  */
        if (fscanf(input2, "%*s %lf", &scatpar->dt) == 0) scatpar->dt = 1.5e-16;                        /* 1.5e-16 s */
        if (fscanf(input2, "%*s %lf", &scatpar->totalTime) == 0) scatpar->totalTime = 10e-12;           /* 10 ps     */
        if (fscanf(input2, "%*s %lf", &scatpar->transientTime) == 0) scatpar->transientTime = 3e-12;    /* 3 ps      */
        if (fscanf(input2, "%*s %lf", &scatpar->averTime) == 0) scatpar->averTime = 2.0e-13;            /* 2.0e-13 s */
        if (fscanf(input2, "%*s %lf", &phys_quantities->inputVs) == 0) phys_quantities->inputVs = 0.0;  /* 0.0 V     */
        if (fscanf(input2, "%*s %lf", &phys_quantities->inputVd) == 0) phys_quantities->inputVd = 0.0;  /* 0.0 V     */
        if (fscanf(input2, "%*s %lf", &constants.conv_omega) == 0) constants.conv_omega = 1.8;          /* 1.8       */
        if (fscanf(input2, "%*s %lf", &constants.conv_tolerance) == 0) constants.conv_tolerance = 1e-4; /* 1e-4      */
        fclose(input2);
    } else {
        printf("ERROR: File '%s' not found!\n", argv[2]);
        free(geometry);
        free(scatpar);
        free(phys_quantities);
        exit(1);
    }
    printf("Times read from file '%s':\n", argv[2]);
    printf("dt\t\t= %.1e s\nTotal time\t= %.1e s\nTransient time\t= %.1e s\nAveraging time\t= %.1e s\n",
           scatpar->dt, scatpar->totalTime, scatpar->transientTime, scatpar->averTime);
    printf("===========================================================\n");
    printf("Voltages read from file '%s':\n", argv[2]);
    printf("Vsource\t= %.5f V\nVdrain\t= %.5f V\n",
           phys_quantities->inputVs, phys_quantities->inputVd);
    printf("===========================================================\n");

    /*=== execute EMC simulator ===*/
    outputdata = EMC(constants, geometry, scatpar, phys_quantities);

    /*=== write data output to files ===*/
    oooWriteOutput(outputdata, geometry, phys_quantities);
    printf("Output written.\n");

    /*=== free memory ===*/
    oooPoissonSOR(constants, geometry, phys_quantities, 2);  /* free memory for poisson coefficients */
    free(geometry);
    free(scatpar);
    free(phys_quantities);

    free(outputdata->current_cumulative);
    free(outputdata->current_from_charge);
    free(outputdata);

    return 0;
}