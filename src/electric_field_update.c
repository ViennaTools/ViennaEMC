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
#include "emc.h"


/********************************************************************/
/* Calculate the electric field at the half-mesh point distance     */
/* fx_field, fy_field are the electric fields used for free flight. */ 
/* (Dev.independent) 						    */
/********************************************************************/
int oooElectricFieldUpdate(const_t constpar, geometry_t *geometry, scatpar_t *scatpar, phys_quant_t *phys_quantities, double *Time)
{
    static int i, j;
    static double sum1, sum2, sum3, sum4, sumn;
    static double func1[MAXNX][MAXNY], func2[MAXNX][MAXNY], func3[MAXNX][MAXNY], func4[MAXNX][MAXNY], n[MAXNX][MAXNY];

    /*=== calculate electric field from potential ===*/
    for (i = 0; i <= geometry->nxmax; ++i)  
      for (j = 0; j <= geometry->nymax; ++j)
      {
	  if (i == geometry->nxmax)
//       	  if ((i == geometry->nxmax)||(i == 0))
	  {
//              phys_quantities->fxField[i][j] = 0;
//              phys_quantities->fyField[i][j] = 0;
              phys_quantities->fxField[i][j] = (phys_quantities->potential[i-1][j] - phys_quantities->potential[i][j])
                                               / geometry->meshSize * constpar.vt;
	      if (j == geometry->nymax)
		  phys_quantities->fyField[i][j] = (phys_quantities->potential[i][j-1] - phys_quantities->potential[i][j]) 
                                                   / geometry->meshSize * constpar.vt;
	      else
		  phys_quantities->fyField[i][j] = (phys_quantities->potential[i][j] - phys_quantities->potential[i][j+1]) 
		                                   / geometry->meshSize * constpar.vt;
	  }
	  else if (j == geometry->nymax)
	  {
	      phys_quantities->fxField[i][j] = (phys_quantities->potential[i][j] - phys_quantities->potential[i+1][j]) 
                                               / geometry->meshSize * constpar.vt;
	      phys_quantities->fyField[i][j] = (phys_quantities->potential[i][j-1] - phys_quantities->potential[i][j]) 
                                               / geometry->meshSize * constpar.vt;
	  }
	  else
	  {
//              if (i==0 && j==0) printf("\n\nphys_quantities->potential[i][j] = %f\n\n", phys_quantities->potential[i][j]);
              
	      phys_quantities->fxField[i][j] = (phys_quantities->potential[i][j] - phys_quantities->potential[i+1][j]) 
                                               / geometry->meshSize * constpar.vt;
	      phys_quantities->fyField[i][j] = (phys_quantities->potential[i][j] - phys_quantities->potential[i][j+1]) 
                                               / geometry->meshSize * constpar.vt;
	  }
      }
    
    /*=== calculate local electric field ===*/
    for (i = 0; i <= geometry->nxmax; ++i)
      for (j = 0; j <= geometry->nymax; ++j)
/*      {
        if (i == geometry->nxmax)
            phys_quantities->localfieldx[i][j] = (phys_quantities->fxField[i-1][j] + phys_quantities->fxField[i][j]) * 0.5;
        else
            phys_quantities->localfieldx[i][j] = (phys_quantities->fxField[i][j] + phys_quantities->fxField[i+1][j]) * 0.5;
        
        if (j == geometry->nymax)
            phys_quantities->localfieldy[i][j] = (phys_quantities->fyField[i][j-1] + phys_quantities->fyField[i][j]) * 0.5;
        else
            phys_quantities->localfieldy[i][j] = (phys_quantities->fyField[i][j] + phys_quantities->fyField[i][j+1]) * 0.5;
      }
*/
      {
        if (i == geometry->nxmax) 
            phys_quantities->localfieldx[i][j] = (phys_quantities->fxField[i][j-1] + phys_quantities->fxField[i][j]) * 0.5;
        else 
            phys_quantities->localfieldx[i][j] = (phys_quantities->fxField[i][j] + phys_quantities->fxField[i][j+1]) * 0.5;
        
        if (j == geometry->nymax) 
            phys_quantities->localfieldy[i][j] = (phys_quantities->fyField[i-1][j] + phys_quantities->fyField[i][j]) * 0.5;
        else 
            phys_quantities->localfieldy[i][j] = (phys_quantities->fyField[i][j] + phys_quantities->fyField[i+1][j]) * 0.5;
      }

    /*=== Calculate instantaneous Fx(i), Fy(i) and SheetDensity(i) ===*/
    /*=== Sum up only after the initial transient ===*/
    if (*Time > scatpar->transientTime) 
    {
      for (i = 0; i <= geometry->nxmax; ++i)   
	for (j = 0; j <= geometry->nymax; ++j) 
	{
	  n[i][j] = phys_quantities->elecDensity[i][j] * constpar.Ni;
	  
	  func1[i][j] = phys_quantities->fxField[i][j] * n[i][j];	
	  func2[i][j] = phys_quantities->fyField[i][j] * n[i][j];
	  func3[i][j] = phys_quantities->fxField[i][j];
	  func4[i][j] = phys_quantities->fyField[i][j];
	}

      for (i = 0; i <= geometry->nxmax; ++i) 
      {
	sum1 = 0.0;
	sum2 = 0.0;
	sum3 = 0.0;
	sum4 = 0.0;
	sumn = 0.0;

	for (j = 0; j <= geometry->nymax; ++j)
	{
	  sum1 += func1[i][j] * geometry->meshSize;
	  sum2 += func2[i][j] * geometry->meshSize;
	  
	  sum3 += func3[i][j] * geometry->meshSize;
	  sum4 += func4[i][j] * geometry->meshSize;

          sumn +=     n[i][j] * geometry->meshSize;
	}

	if (sumn != 0.0)
	{
	  phys_quantities->fieldX_avg[i] += sum1 / sumn; 
	  phys_quantities->fieldY_avg[i] += sum2 / sumn; 
	}
	phys_quantities->fieldX_avg2[i] += sum3; 
	phys_quantities->fieldY_avg2[i] += sum4;
	phys_quantities->densityX_avg[i] += sumn; 
      }
    }

    return 0;
} 
