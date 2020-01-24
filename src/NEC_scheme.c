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

#include <math.h>
#include "emc.h"
#include <omp.h>

/********************************************************************/
/* NEC method for the charge assignment at the node points          */
/********************************************************************/
int oooChargeAssignmentNEC(const_t constpar, geometry_t *geometry, scatpar_t *scatpar, el_data_t *particles, phys_quant_t *phys_quantities)
{
  static int i, j, n;
  static double denn, teglo, elecNumber[MAXNX][MAXNY], normFactor;

  /*=== Evaluate multiplication constant ===*/
  normFactor = 1.0 / (geometry->cellVolume * constpar.Ni);

  /*=== Reset the charge vector ===*/
  for (i = 0; i <= geometry->nxmax; ++i)   
    for (j = 0; j <= geometry->nymax; ++j) 
      elecNumber[i][j] = 0.0;

//  #pragma omp parallel
//  #pragma omp for
  #pragma omp for schedule(static, 1)
  /*=== Charge assignment part ===*/
  for (n = 0; n <= scatpar->n_used; ++n) 
  {
    i = (int) (particles[n].p[5] / geometry->meshSize); 
    j = (int) (particles[n].p[6] / geometry->meshSize); 

    teglo = particles[n].p[7] * 0.25; 
    elecNumber[i]  [j]   += teglo;
    elecNumber[i]  [j+1] += teglo;
    elecNumber[i+1][j]   += teglo;
    elecNumber[i+1][j+1] += teglo;
  }

  /*=== Calculate electron density ===*/
  for (i = 0; i <= geometry->nxmax; ++i)   
    for (j = 0; j <= geometry->nymax; ++j) 
    {
      denn = elecNumber[i][j];
      if (i == 0 || i == geometry->nxmax)  denn *= 2.0;
      if (j == 0 || j == geometry->nymax)  denn *= 2.0;
      phys_quantities->elecDensity[i][j] = denn * normFactor;
    }

  return 0;
}

