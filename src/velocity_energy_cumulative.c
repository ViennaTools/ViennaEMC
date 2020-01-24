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
#include <math.h>
#include "emc.h"
#include <omp.h>


/********************************************************************/
/*  Calculate cumulative velocities and energies of the particles   */
/*  				             (device independent)   */
/********************************************************************/
currents_t oooVelocityEnergyCumulative(const_t constpar, geometry_t *geometry, scatpar_t *scatpar, el_data_t * particles, phys_quant_t *phys_quantities, double *Time, int *nTimeStepsAv)
{
  static int i, j, n, iv, weight;
  static double velx, vely, ee, denom, factor, qfactor, momentaryCur[MAXNX];
  static double velxSum[MAXNX], velySum[MAXNY], eNumber[MAXNX], enerSum[MAXNX], currentSum[MAXNX];	/* 1D sums */
  static double velxSum_2D[MAXNX][MAXNY], velySum_2D[MAXNX][MAXNY], eNumber_2D[MAXNX][MAXNY];  		/* 2D sums */

  static double counter = 0; /* counter = kcount1.kcounter */

  currents_t currents;

  /*=== Calculate instantaneous v(i), e(i) and currentdensity(i) ===*/
  ++counter;

  for (i = 0; i <= geometry->nxmax; ++i)
  {
    eNumber[i] = 0.0; 
    velxSum[i] = 0.0; 
    velySum[i] = 0.0; 
    enerSum[i] = 0.0;
    for (j = 0; j <= geometry->nymax; ++j)  
    {
	eNumber_2D[i][j] = 0;
	velxSum_2D[i][j] = 0.0;
	velySum_2D[i][j] = 0.0;
    }
  }

//  #pragma omp parallel
//  #pragma omp for
  #pragma omp for schedule(static, 1)
  for (n = 0; n <= scatpar->n_used; ++n) 
  {
//    i = (rintf) (particles[n].p[5] / geometry->meshSize);
    i = (int) (particles[n].p[5] / geometry->meshSize);
    j = (int) (particles[n].p[6] / geometry->meshSize);
    if      (i < 0)                i = 0;
    else if (i > geometry->nxmax)  i = geometry->nxmax;

    weight = particles[n].p[7];
    iv     = particles[n].ip;
    ee     = particles[n].energy;

    denom  = 1.0 / (constpar.af2 * ee + 1.0);
    if (iv == 1)
    {
	velx = constpar.hm[0] * particles[n].p[1] * denom;
        vely = constpar.hm[1] * particles[n].p[2] * denom;
    }  
    else if (iv == 2)
    {
	velx = constpar.hm[1] * particles[n].p[1] * denom;
	vely = constpar.hm[0] * particles[n].p[2] * denom;
    }
    else if (iv == 3)
    {
	velx = constpar.hm[2] * particles[n].p[1] * denom;
	vely = constpar.hm[1] * particles[n].p[2] * denom;
    }
    else 
    {
	printf("invalid iv = %d at particle index n = %d (n used = %d)\n", iv, n, scatpar->n_used);
    }

    velxSum[i] += velx * weight; 
    velySum[i] += vely * weight; 
    enerSum[i] += ee   * weight; 
    eNumber[i] +=        weight; 
    velxSum_2D[i][j] += velx * weight;
    velySum_2D[i][j] += vely * weight;
    eNumber_2D[i][j] += weight;
  }

  qfactor = constpar.q / geometry->meshSize;

  if (*Time < scatpar->transientTime)
  {
    for (i = 0; i <= geometry->nxmax; ++i) 
    {
      currentSum[i]  += velxSum[i] * qfactor;
      momentaryCur[i] = velxSum[i] * qfactor;

    }
  }
  else if (*Time == scatpar->transientTime)
  {
    counter = 1;
    for (i = 0; i <= geometry->nxmax; ++i) 
    {
      currentSum[i] = 0.0;
      phys_quantities->curSumX[i][j] = 0.0;
      phys_quantities->curSumY[i][j] = 0.0;

      currentSum[i]   += velxSum[i] * qfactor;
      momentaryCur[i]  = velxSum[i] * qfactor;
      
      if (eNumber[i] != 0.0)
      {
	factor = 1.0 / eNumber[i]; 		               	/* number of electrons */
	phys_quantities->velxSum[i] += velxSum[i] * factor; 	/* mean velocity x     */
	phys_quantities->velySum[i] += velySum[i] * factor; 	/* mean velocity y     */
	phys_quantities->enerSum[i] += enerSum[i] * factor; 	/* mean energy	       */      
      }
    }
    for (j = 0; j <= geometry->nymax; ++j)   
      for (i = 0; i <= geometry->nxmax; ++i) 
      {
	  phys_quantities->elSum[i][j]    = eNumber_2D[i][j];
	  phys_quantities->curSumX[i][j] += velxSum_2D[i][j] * qfactor;
	  phys_quantities->curSumY[i][j] += velySum_2D[i][j] * qfactor;
      }
  }
  else /* if (*Time > scatpar->transientTime) */
  {
    for (i = 0; i <= geometry->nxmax; ++i) 
    {
      currentSum[i]  += velxSum[i] * qfactor;
      momentaryCur[i] = velxSum[i] * qfactor;
      
      if (eNumber[i] != 0.0)
      {
	factor = 1.0 / eNumber[i];
	phys_quantities->velxSum[i] += velxSum[i] * factor;
	phys_quantities->velySum[i] += velySum[i] * factor;
	phys_quantities->enerSum[i] += enerSum[i] * factor;
      }
    }
    for (j = 0; j <= geometry->nymax; ++j)   
      for (i = 0; i <= geometry->nxmax; ++i) 
      {
	  phys_quantities->elSum[i][j] = (phys_quantities->elSum[i][j] * (counter - 1)
				       + eNumber_2D[i][j]) / counter;
	  phys_quantities->curSumX[i][j] += velxSum_2D[i][j] * qfactor;
	  phys_quantities->curSumY[i][j] += velySum_2D[i][j] * qfactor;
      }
  }
  
  /*=== calculate cumulative and momentary source & drain currents ===*/ 
  denom = 1.0 / (counter * geometry->deviceWidth);

  currents.Is_cumul     = currentSum[geometry->nxmax / 2] * denom;
  currents.Id_cumul     = currentSum[geometry->nxmax / 2] * denom;
//  currents.Is_cumul     = currentSum[geometry->ns + 3] * denom;
//  currents.Id_cumul     = currentSum[geometry->nd - 3] * denom;
//  currents.Is_momentary = momentaryCur[geometry->ns + 3] / geometry->deviceWidth;
//  currents.Id_momentary = momentaryCur[geometry->nd - 3] / geometry->deviceWidth;

  return currents;
}
