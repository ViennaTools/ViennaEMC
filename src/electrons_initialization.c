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


/********************************************************************/
/* Initialize number of electrons based on the equilibrium solution */
/* of the Poisson equation  (Uniform mesh assumed, Dev.independent) */
/********************************************************************/
int oooElectronsInitialization(const_t constpar, geometry_t *geometry, scatpar_t *scatpar, el_data_t *particles, phys_quant_t *phys_quantities)
{
  static int i, j, ne, n_ij;
  static double denn;

  ne = 0; 
  for (i = 0; i <= geometry->nxmax; ++i)  
    for (j = 0; j <= geometry->nymax; ++j)
    {
//      denn = exp(phys_quantities->potential[i][j] * constpar.vt) * constpar.Ni;
      denn = scatpar->dopingCon;
//      printf("denn = %f\n", denn);
      if (i == 0 || i == geometry->nxmax)  denn *= 0.5;
      if (j == 0 || j == geometry->nymax)  denn *= 0.5;
      if (i == geometry->nxmax) denn *= 1.01;
      n_ij = (rintf) (denn * geometry->cellVolume);
//      printf("denn * geometry->cellVolume = %f\n", denn * geometry->cellVolume);
//      printf("n_ij = %i\n", n_ij);
      while(n_ij > 0.0)
      {
	if (ne > MAXEN)
	  printf("Actual number of electrons = %i exceeds %i\n",ne,MAXEN);
	oooInitKspaceFD(constpar, geometry, scatpar, particles, &ne, &i, &j);
	oooInitRealspace(geometry, particles, &ne, &i, &j);
	n_ij -= particles[ne].p[7];
	++ne;
      }
    }
  
  /* set index of last used electron (scatpar->n_used) */
  scatpar->n_used = ne - 1;	
  
  printf("Number of electrons max allowed = %i\n",MAXEN);
  printf("Number of electrons initialized = %i",scatpar->n_used);
  if (scatpar->n_used == -1) printf(" / (%.3e)", denn*geometry->cellVolume);

  for (i = scatpar->n_used + 1; i <= MAXEN; ++i)
    particles[i].ip = 9;

  return 0;
} 
