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

/********************************************************************/
/* Writes output files pre-processing, averages and post-processing */
/********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "emc.h"


/* this function calculates the time averaged potential */
int oooAveragePotential(geometry_t *geometry, phys_quant_t *phys_quantities, int mode, double timesteps)
{
  int i,j;
  static double avgpot[MAXNX][MAXNY];

  if (mode == 0) /* initialization */
  {
    for (i = 0; i < MAXNX; i++)
      for (j = 0; j < MAXNY; j++)
        avgpot[i][j] = 0.0;
    return 0;
  }
  if (mode == 1) /* sum up */
  {
    for (i = 0; i <= geometry->nxmax; i++)
      for (j = 0; j <= geometry->nymax; j++)
        avgpot[i][j] += phys_quantities->potential[i][j];
    return 0;
  }
  if (mode == 2) /* divide & write back */
  {
    for (i = 0; i <= geometry->nxmax; i++)
      for (j = 0; j <= geometry->nymax; j++)
      {
        avgpot[i][j] /= timesteps;
        phys_quantities->potential[i][j] = avgpot[i][j]; 
      }
    return 0;
  }
  return 0;
}

/* this function calculates the time averaged electron density */
int oooAverageDensity(geometry_t *geometry, phys_quant_t *phys_quantities, int mode, double timesteps)
{
  int i,j;
  static double dens[MAXNX][MAXNY];

  if (mode == 0) /* initialization */
  {
    for (i = 0; i < MAXNX; i++)
      for (j = 0; j < MAXNY; j++)
         dens[i][j] = 0.0;
    return 0;
  }
  if (mode == 1) /* sum up */
  {
    for (i = 0; i <= geometry->nxmax; i++)
      for (j = 0; j <= geometry->nymax; j++)
        dens[i][j] += phys_quantities->elecDensity[i][j];
    return 0;
  }
  if (mode == 2) /* divide & write back */
  {
    for (i = 0; i <= geometry->nxmax; i++)
      for (j = 0; j <= geometry->nymax; j++)
      {
        dens[i][j] /= timesteps;
        phys_quantities->elecDensity[i][j] = dens[i][j]; 
      }
    return 0;
  }
  return 0;
}

int oooWritePotential(const_t constpar, geometry_t *geometry, phys_quant_t *phys_quantities, int flagSolution, int iter)
{
  static int i, j;
  static double denn, r1, sum;

  FILE *output1, *output2, *output3;

  char file1[100];
  sprintf(file1, "init_%d_potential.csv", iter);
  char file2[100];
  sprintf(file2, "init_%d_electron_density.csv", iter);
  char file3[100];
  sprintf(file3, "init_%d_sheet_density_x.csv", iter);

  output1 = fopen(file1, "w");
  output2 = fopen(file2, "w");
  output3 = fopen(file3, "w");

//  output1 = fopen("init_potential.csv", "w");
//  output2 = fopen("init_electron_density.csv", "w");
//  output3 = fopen("init_sheet_density_x.csv", "w");

  for (i = 0; i <= geometry->nxmax; ++i)
  {
      sum = 0.0;
      for (j = 0; j <= geometry->nymax; ++j)   
      {
          /*=== Electron Density Calculation ===*/
	  if (flagSolution == 0)
	      denn = exp(phys_quantities->potential[i][j]); 
	  else if (flagSolution == 1)
          {
	      denn = phys_quantities->elecDensity[i][j] * constpar.Ni;
              if (phys_quantities->elecDensity[i][j] == 0.0) denn = constpar.Ni;     
	  }
	  fprintf(output2,"%d %d %f \n",i, j, denn); 

	  /*=== Potential Calculation ===*/
//          r1 = (constpar.delta_Ec - phys_quantities->potential[i][j]) * constpar.vt;
          r1 = phys_quantities->potential[i][j] * constpar.vt;
	  fprintf(output1,"%d %d %f \n",i, j, r1); 

          /*=== Sheet Density Calculation ===*/
          sum += phys_quantities->elecDensity[i][j] * constpar.Ni * geometry->gridy[j + 1] * constpar.debyeLength;
//          sum += phys_quantities->elecDensity[i][j] * geometry->gridy[j + 1] * constpar.debyeLength;
      }
      fprintf(output3,"%.6e \n", sum * 1e-4);
  }

  fclose(output1);
  fclose(output2);
  fclose(output3); 

  return 0;
}


/********************************************************************/
/* calculate all output quantities                                  */ 
/* 		                                  (Dev.independent) */
/********************************************************************/
output_t *oooCalculateOutput(const_t constpar, geometry_t *geometry, phys_quant_t *phys_quantities, int nTimeStepsAv, int simtime, int iterTotal, currents_t *current_cumulative, curr_from_charge_t *current_from_charge)
//output_t *oooCalculateOutput(const_t constpar, geometry_t *geometry, phys_quant_t *phys_quantities, int nTimeStepsAv, int simtime, int iterTotal)
{
  int i,j;
  double sum;
  double xpos = 0.0; 
  double ypos = 0.0;
  output_t *outputdata;

  /*=== allocate memory for output data ===*/
  outputdata = (output_t *) calloc(1, sizeof(output_t));

  outputdata->iterTotal = iterTotal;
  outputdata->current_cumulative = current_cumulative;
  outputdata->current_from_charge = current_from_charge;

  /*=== save max mesh coordinates ===*/
  outputdata->x_max = geometry->nxmax;
  outputdata->y_max = geometry->nymax;

  outputdata->totaltime = simtime;

  for (i = 0; i <= geometry->nxmax; i++)
  {
    sum = 0.0;

    /*=== x-axis calculation ===*/
    outputdata->x_axis[i] = xpos;
    xpos += geometry->gridx[i + 1] * constpar.debyeLength;

    for (j = 0; j <= geometry->nymax; j++)
    {
      /*=== y-axis calculation ===*/
      outputdata->y_axis[j] = ypos;
      ypos += geometry->gridy[j + 1] * constpar.debyeLength;

      /*=== Potential Calculation ===*/
//      outputdata->potential[i][j] = (constpar.delta_Ec - phys_quantities->potential[i][j]) * constpar.vt;
      outputdata->potential[i][j] = phys_quantities->potential[i][j] * constpar.vt;
//      outputdata->potential[i][j] = (constpar.delta_Ec - phys_quantities->potential[i][j]) * constpar.vt;
//      outputdata->potential[i][j] = (phys_quantities->potential[i][j]) * constpar.vt;
          
      /*=== Electron Density Calculation ===*/
      if (phys_quantities->elecDensity[i][j] == 0.0) outputdata->electrondensity[i][j] = constpar.Ni;
        else outputdata->electrondensity[i][j] = phys_quantities->elecDensity[i][j] * constpar.Ni;

      /*=== Sheet Density Calculation ===*/
      sum += phys_quantities->elecDensity[i][j] * constpar.Ni * geometry->gridy[j + 1] * constpar.debyeLength;
//      sum += phys_quantities->elecDensity[i][j] * geometry->gridy[j + 1] * constpar.debyeLength;
  
      /*=== Electric Field Calculation ===*/
      if (i == geometry->nxmax)
      {
//	outputdata->fieldXY_x[i][j] = 0.0;
       outputdata->fieldXY_x[i][j] = -(phys_quantities->potential[i-1][j] - phys_quantities->potential[i][j]) 
	   / geometry->meshSize * constpar.vt;
	if (j == geometry->nymax) outputdata->fieldXY_y[i][j] = 0.0;
	else
	  outputdata->fieldXY_y[i][j] = -(phys_quantities->potential[i][j] - phys_quantities->potential[i][j+1]) 
	     / geometry->meshSize * constpar.vt;	      
      }
      else if (j == geometry->nymax)
      {
	outputdata->fieldXY_x[i][j] = -(phys_quantities->potential[i][j] - phys_quantities->potential[i+1][j]) 
	   / geometry->meshSize * constpar.vt;
	outputdata->fieldXY_y[i][j] = 0.0;
      }
      else
      {
	outputdata->fieldXY_x[i][j] = -(phys_quantities->potential[i][j] - phys_quantities->potential[i+1][j]) 
	   / geometry->meshSize * constpar.vt;
	outputdata->fieldXY_y[i][j] = -(phys_quantities->potential[i][j] - phys_quantities->potential[i][j+1]) 
	   / geometry->meshSize * constpar.vt; 
      }

      /*=== Current Density Calculation ===*/
      outputdata->currentdensityX[i][j] = phys_quantities->curSumX[i][j] / ((double)(nTimeStepsAv) * geometry->deviceWidth);
      outputdata->currentdensityY[i][j] = phys_quantities->curSumY[i][j] / ((double)(nTimeStepsAv) * geometry->deviceWidth);
      outputdata->currentdensity[i][j]  = sqrt(phys_quantities->curSumX[i][j] * phys_quantities->curSumX[i][j] + 
		     		              phys_quantities->curSumY[i][j] * phys_quantities->curSumY[i][j]) / ((double)(nTimeStepsAv) * geometry->deviceWidth);	

    }
    outputdata->sheetdensity[i] = sum * 1e-4; 

    outputdata->fieldX[i]  = (phys_quantities->fieldX_avg[i] / (double)(nTimeStepsAv)) / 1e2;
    outputdata->fieldY[i]  = (phys_quantities->fieldY_avg[i] / (double)(nTimeStepsAv)) / 1e2;
    outputdata->fieldX2[i] = (phys_quantities->fieldX_avg2[i] / (double)(nTimeStepsAv)) / 1e2;
    outputdata->fieldY2[i] = (phys_quantities->fieldY_avg2[i] / (double)(nTimeStepsAv)) / 1e2;
    outputdata->density[i] = (phys_quantities->densityX_avg[i] / (double)(nTimeStepsAv)) / 1e4;

    outputdata->velocityX[i] = phys_quantities->velxSum[i] / (double)(nTimeStepsAv);
    outputdata->velocityY[i] = phys_quantities->velySum[i] / (double)(nTimeStepsAv);
    outputdata->energy[i]    = phys_quantities->enerSum[i] / (double)(nTimeStepsAv);
  }
  
  return outputdata;
}


/********************************************************************/
/* write all output quantities into files                           */ 
/* 		                                  (Dev.independent) */
/********************************************************************/
//int oooWriteOutput(output_t *outputdata)
int oooWriteOutput(output_t *outputdata, geometry_t *geometry, phys_quant_t *phys_quantities)
{
  int i,j;

  FILE *output1, *output2, *output3, *output4, *output5, *output6, *output7, 
       *output8, *output9, *output10, *output11, *output12, *output13, *output14; 

  output1 = fopen("x_axis.csv", "w");
  output2 = fopen("y_axis.csv", "w");
  output3 = fopen("potential.csv", "w");
  output4 = fopen("electron_density.csv", "w");
  output5 = fopen("sheet_density_x.csv", "w");
  output6 = fopen("fields_density_x.csv", "w");
  output7 = fopen("fieldXY.csv", "w");
  output8 = fopen("v_e_aver.csv", "w");
  output9 = fopen("current_densityX.csv", "w");
  output10 = fopen("current_densityY.csv", "w");
  output11 = fopen("current_density.csv", "w");
  output12 = fopen("total_simulation_time.csv", "w");

  output13 = fopen("cur_from_charge_SD.csv", "w");
  output14 = fopen("current_cumulative.csv", "w");

  fprintf(output12,"%d\n", outputdata->totaltime);
  fclose(output12);

  for (i = 0; i <= outputdata->x_max; ++i)
  {
    fprintf(output1, "%e \n", outputdata->x_axis[i]);

    fprintf(output5,"%e \n", outputdata->sheetdensity[i]);
    fprintf(output6,"%e \t %e \t %e \t %e \t %e \n", outputdata->fieldX[i], outputdata->fieldY[i], outputdata->fieldX2[i], outputdata->fieldY2[i], outputdata->density[i]); 
    fprintf(output8, "%f \t %f \t %f \n", outputdata->velocityX[i], outputdata->velocityY[i], outputdata->energy[i]);

    for (j = 0; j <= outputdata->y_max; ++j)
    {
      if (i == 0) fprintf(output2, "%e \n", outputdata->y_axis[j]);

      fprintf(output3,"%d %d %.10f \n",i, j, outputdata->potential[i][j]);
      fprintf(output4,"%d %d %f \n",i, j, outputdata->electrondensity[i][j]);

      fprintf(output7, "%d %d %e %e \n", i, j, outputdata->fieldXY_x[i][j], outputdata->fieldXY_y[i][j]);

      fprintf(output9, "%d %d %f \n",i, j, outputdata->currentdensityX[i][j]);
      fprintf(output10,"%d %d %f \n",i, j, outputdata->currentdensityY[i][j]);
      fprintf(output11,"%d %d %f \n",i, j, outputdata->currentdensity[i][j]);
    }
  }

  for (i = 0; i < outputdata->iterTotal; i++)
  {
    fprintf(output14, "%f \t %f \t %f \t %f \t %f\n", outputdata->current_from_charge[i].Time, outputdata->current_cumulative[i].Is_cumul,       outputdata->current_cumulative[i].Id_cumul, 
                                                                                               outputdata->current_cumulative[i].Is_momentary,   outputdata->current_cumulative[i].Id_momentary);
    fprintf(output13, "%f \t %f \t %f \t %f \t %f \n", outputdata->current_from_charge[i].Time, outputdata->current_from_charge[i].sourceCurr,   outputdata->current_from_charge[i].drainCurr, 
                                                                                                outputdata->current_from_charge[i].sourceFactor, outputdata->current_from_charge[i].drainFactor);
  }
  
  printf("Current = %.8e, Resistivity = %f, Rho = %.8e\n", 
          outputdata->current_from_charge[outputdata->iterTotal - 1].sourceCurr,
          (phys_quantities->inputVd - phys_quantities->inputVs) / outputdata->current_from_charge[outputdata->iterTotal - 1].sourceCurr,
          (phys_quantities->inputVd - phys_quantities->inputVs) * geometry->deviceHeight * geometry->deviceWidth / 
            (geometry->deviceLength * outputdata->current_from_charge[outputdata->iterTotal - 1].sourceCurr));

  fclose(output1);
  fclose(output2);
  fclose(output3);
  fclose(output4);
  fclose(output5);
  fclose(output6);
  fclose(output7);
  fclose(output8);
  fclose(output9);
  fclose(output10);
  fclose(output11);

  fclose(output13);
  fclose(output14);

  return 0;
}