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
int main(int argc, char* argv[])
{
  const_t constants;
  geometry_t *geometry;
  scatpar_t *scatpar;
  phys_quant_t *phys_quantities;
  output_t *outputdata;
  int threads;

  FILE *input1, *input2;

  if (argc != 3)
  {
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
  if (input1 != NULL)
  {										/* default:   */
    if (fscanf(input1,"%lf",&geometry->deviceLength) == 0) geometry->deviceLength = 600e-9;	/*  1 um      */
    if (fscanf(input1,"%lf",&geometry->deviceHeight) == 0) geometry->deviceHeight = 100e-9;
    if (fscanf(input1,"%lf",&geometry->deviceWidth) == 0) geometry->deviceWidth = 1e-6;		/*  1 um      */
    if (fscanf(input1,"%lf",&geometry->meshSize) == 0) geometry->meshSize = 1e-9;		/*  1 nm      */
    if (fscanf(input1,"%d",&scatpar->kTLevels) == 0) scatpar->kTLevels = 2;		/*  1 nm      */
 //   if (fscanf(input1,"%lf",&scatpar->dopingCon) == 0) scatpar->dopingCon =  5e25;            /*  5e25 m^-3 */
//    if (fscanf(input1,"%lf",&scatpar->dopingCon) < 10 ) 
        scatpar->dopingCon =  oooIntegrate(constants, scatpar->kTLevels); 	/*  5e25 m^-3 */
                                                        constants.Ni = 1.0 * scatpar->dopingCon;
                                                        constants.debyeLength /= sqrt(constants.Ni);
    if (fscanf(input1,"%d",&scatpar->precision) == 0) scatpar->precision = 10;                 /* 10 steps for e- energies */
    if (fscanf(input1,"%d",&threads) == 0) threads = 2;                 /* 10 steps for e- energies */
    omp_set_num_threads(threads);
    
    printf("\n\nee = %.5e\n\n", constants.tau_ee);
    
//    FILE *myDopingData;
//    myDopingData = fopen("myDopingData.txt", "w");
//    for (int i = 1; i < 30; i++)
//    {
//        fprintf(myDopingData, "%d, %fl\n", i, oooIntegrate(constants, i));
//    }
//    fclose(myDopingData);

//    double delta_e = constants.kt * scatpar->kTLevels;
//    scatpar->energy_levels = calloc(scatpar->kTLevels, sizeof(double));
//    scatpar->cumulative = calloc(scatpar->kTLevels, sizeof(double));
//      printf("scatpar->precision = %d\n", scatpar->precision);
//    for (int i=0;i<scatpar->precision;i++) {
//      scatpar->energy_levels[i] = constants.ef - (0.5*delta_e + i * delta_e / scatpar->precision) / constants.q;
//      scatpar->cumulative[i] = scatpar->energy_levels[i] - scatpar->energy_levels[0]
//              + log((exp(scatpar->energy_levels[i])+1) / (exp(scatpar->energy_levels[0])+1));
//      printf("energy_levels = %f, cumulative = %f\n", scatpar->energy_levels[i], scatpar->cumulative[i]);
//    }

//                                                         scatpar->dopingCon = 1e22;
//                                                         constants.Ni = 1e16;
//                                                         constants.debyeLength /= sqrt(constants.Ni);

    fclose(input1);
  }
  else
  {
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
  input2 = fopen(argv[2],"r");
  if (input2 != NULL)
  {											        	/* default:  */
    if (fscanf(input2,"%lf",&scatpar->dt) == 0) scatpar->dt = 1.5e-16;         		       	 	/* 1.5e-16 s */
    if (fscanf(input2,"%lf",&scatpar->totalTime) == 0) scatpar->totalTime = 10e-12;              	/* 10 ps     */
    if (fscanf(input2,"%lf",&scatpar->transientTime) == 0) scatpar->transientTime = 3e-12;     		/* 3 ps      */
    if (fscanf(input2,"%lf",&scatpar->averTime) == 0) scatpar->averTime = 2.0e-13;             	 	/* 2.0e-13 s */
    if (fscanf(input2,"%lf",&phys_quantities->inputVs) == 0) phys_quantities->inputVs = 0.0;            /* 0.0 V     */
    if (fscanf(input2,"%lf",&phys_quantities->inputVd) == 0) phys_quantities->inputVd = 0.0;//1.0;      /* 1.0 V     */
    if (fscanf(input2,"%lf",&constants.conv_omega) == 0) constants.conv_omega = 1.8;	  		/* 1.8       */
    if (fscanf(input2,"%lf",&constants.conv_tolerance) == 0) constants.conv_tolerance = 1e-4;   	/* 1e-4      */
    fclose(input2);
  }
  else
  {
    printf("ERROR: File '%s' not found!\n", argv[2]);
    free(geometry);
    free(scatpar);
    free(phys_quantities);
    exit(1);
  }
  printf("Times read from file '%s':\n", argv[2]);
  printf("dt             = %.1e s\nTotal time     = %.1e s\nTransient time = %.1e s\nAveraging time = %.1e s\n",
	 scatpar->dt, scatpar->totalTime, scatpar->transientTime, scatpar->averTime);
  printf("===========================================================\n");
  printf("Voltages read from file '%s':\n", argv[2]);
  printf("Vsource    = %.5f V\nVdrain     = %.5f V\n",
	 phys_quantities->inputVs, phys_quantities->inputVd);
  printf("===========================================================\n");

  /*=== execute EMC simulator ===*/
  outputdata = EMC(constants, geometry, scatpar, phys_quantities);

  /*=== write data output to files ===*/
  oooWriteOutput(outputdata, geometry, phys_quantities);
//  oooWriteOutput(outputdata);
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