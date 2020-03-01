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
#include <stdio.h>
#include "emc.h"
#include <omp.h>


/********************************************************************/
/* Calculate the number of carriers in the source and drain regions */
/* based on doping.                               (device specific) */
/********************************************************************/
SDcarriers_t oooSDcarrierNumber(const_t constpar, geometry_t *geometry, phys_quant_t *phys_quantities)
{
  int i, j;
  double denn;
  SDcarriers_t SDcarriers;

  /* source contact */
  i = 0;
  for (j = 0; j <= geometry->nymax; ++j)
  {
    denn = phys_quantities->doping * constpar.Ni * 0.5;//0.5;
//    denn = phys_quantities->doping[i][j] * 0.5;//0.5;
    if (j == 0 || j == geometry->nymax)  denn *= 0.5;

    SDcarriers.nScarriers[j] = (int) (denn * geometry->cellVolume + 0.5);
//    printf("Source Carriers at mesh point x = %d: %d\n", j, SDcarriers.nScarriers[j]);  /*debug*/
  }

  /* drain contact */
  i = geometry->nxmax;
  for (j = 0; j <= geometry->nymax; ++j)
  {
    denn = phys_quantities->doping * constpar.Ni * 0.5;
//    denn = phys_quantities->doping[i][j] * 0.5; 
    if (j == 0 || j == geometry->nymax)  denn *= 0.5;

    SDcarriers.nDcarriers[j] = (int) (denn * geometry->cellVolume + 0.5);
//    printf("Drain Carriers at mesh point x = %d: %d\n", j, SDcarriers.nDcarriers[j]);   /*debug*/
  }

  return SDcarriers;
}


/********************************************************************/
/* Check the charge neutrality of source and drain contacts (FET)   */
/********************************************************************/
int oooCheckSourceDrainContacts(const_t constpar, geometry_t *geometry, scatpar_t *scatpar, el_data_t *particles, phys_quant_t *phys_quantities)
{
  static int j, n, ne, ix, jy, nEleDiff, npts[MAXNY], nptd[MAXNY];
  static SDcarriers_t SDcarriers;
  static int SDcarriers_init = 0;

  /*=== Calculate number of carriers in source and drain regions (called only once) ===*/
  if (SDcarriers_init == 0)
  {
    SDcarriers = oooSDcarrierNumber(constpar, geometry, phys_quantities);
    SDcarriers_init = 1;
  }

  for (j = 0; j <= geometry->nymax; ++j)
  {
    npts[j] = 0;
    nptd[j] = 0;
  }

//  #pragma omp parallel
//  #pragma omp for
  #pragma omp for schedule(static, 1)
//printf("0. npt[%i] = %i\n", 10, npt[10]);
  for (n = 0; n <= scatpar->n_used; ++n)
    if (particles[n].ip != 9)
    {
      ix = (int) (particles[n].p[5] / geometry->meshSize + 0.5); /* x position */
      jy = (int) (particles[n].p[6] / geometry->meshSize + 0.5); /* y position */

      /*=== Delete extra carriers at the source and drain contacts ===*/

      /* source contact */
//      ix = 0;
      if ((ix == 0) && (jy >= 0 && jy <= geometry->nymax))
      {
//printf("1. scatpar->n_used = %i, SDcarriers.nScarriers[jy]=%i\n", scatpar->n_used, SDcarriers.nScarriers[jy]);
//          printf("npts[%i] = %i\n",jy,npts[jy]);
//          printf("nScarriers[%i] = %i\n",jy,SDcarriers.nScarriers[jy]);
	if (npts[jy] < SDcarriers.nScarriers[jy])
	  npts[jy] += particles[n].p[7];
	else
	{
//          printf("del - nScarriers[%i] = %i, ix = %i, npts[%i] = %i\n",jy, SDcarriers.nScarriers[jy], ix, jy,npts[jy]);
	  particles[n].ip = 9;
	  phys_quantities->Iss_eli += particles[n].p[7];
	}
      }

      /* drain contact */
//      ix == geometry->nxmax;
      if ((ix == geometry->nxmax) && (jy >= 0 && jy <= geometry->nymax))
      {
	if (nptd[jy] < SDcarriers.nDcarriers[jy])
	  nptd[jy] += particles[n].p[7];
	else
	{
//          printf("del - nDcarriers[%i] = %i, ix = %i, nptd[%i] = %i\n",jy, SDcarriers.nDcarriers[jy], ix, jy,nptd[jy]);
	  particles[n].ip = 9;
	  phys_quantities->Idd_eli += particles[n].p[7];
	}
      }
    }
//          printf("Iss_eli = %f\n", phys_quantities->Iss_eli);
//printf("1. npts[%i] = %i\n", 10, npts[10]);
//printf("1. nptd[%i] = %i\n", 10, npts[10]);
//printf("1. scatpar->n_used = %i\n", scatpar->n_used);


  /*=== Create carriers at the source and drain contacts ===*/

  /* source contact */
  ix = 0;
//printf("SC_a scatpar->n_used = %i\n", scatpar->n_used);
//        printf("\n\n");
  for (jy = 0; jy <= geometry->nymax; ++jy)
  {
      nEleDiff = SDcarriers.nScarriers[jy] - npts[jy];
      if (nEleDiff < 0)  continue;

      ne = scatpar->n_used;
      while(nEleDiff > 0)
      {
//      printf("nScarriers[%i] = %i, ix = %i, npts[%i] = %i\n",jy, SDcarriers.nScarriers[jy], ix, jy,npts[jy]);
        ++ne;
	oooInitKspaceFD(constpar, geometry, scatpar, particles, &ne, &ix, &jy);
	oooInitRealspace(geometry, particles, &ne, &ix, &jy);
	nEleDiff -= particles[ne].p[7];
        phys_quantities->Iss_cre += particles[ne].p[7];
      }
      scatpar->n_used = ne;

  }
//          printf("Iss_cre = %f\n", phys_quantities->Iss_eli);

  /* drain contact */
  ix = geometry->nxmax;
  for (jy = 0; jy <= geometry->nymax; ++jy)
  {
      nEleDiff = SDcarriers.nDcarriers[jy] - nptd[jy];
      if (nEleDiff < 0)  continue;

      ne = scatpar->n_used;
      while(nEleDiff > 0)
      {
//        printf("nDcarriers[%i] = %i, ix = %i, nptd[%i] = %i\n",jy, SDcarriers.nDcarriers[jy], ix, jy,nptd[jy]);
	++ne;
	oooInitKspaceFD(constpar, geometry, scatpar, particles, &ne, &ix, &jy);
	oooInitRealspace(geometry, particles, &ne, &ix, &jy);
	nEleDiff -= particles[ne].p[7];
        phys_quantities->Idd_cre += particles[ne].p[7];
      }
      scatpar->n_used = ne;
  }
  return 0;
}


/********************************************************************/
/* Delete extra particles if their ip-number iv=9 (Dev.independent) */
/********************************************************************/
int oooDeleteParticles(scatpar_t *scatpar, el_data_t *particles) {
    static int i, j, nFix, flagConv;

#pragma omp parallel for shared(scatpar, particles) private(flagConv, nFix)
    for (i = 0; i <= scatpar->n_used; ++i) {
    if (particles[i].ip == 9) {
        flagConv = 0;
        nFix = scatpar->n_used + 1;

        while (!flagConv) {
            if (particles[nFix].ip != 9) {
                particles[i].ip = particles[nFix].ip;

                for (j = 0; j <= 7; ++j)
                    particles[i].p[j] = particles[nFix].p[j];

                particles[i].energy = particles[nFix].energy;
                particles[nFix].ip = 9;
                --nFix;
                flagConv = 1;
            }

            --nFix;
            if (nFix < i) flagConv = 1;
        }
    }
}

  /* set n_used to last valid particle */
  while(1)
  {
//printf("scatpar->n_used = %i...\n", scatpar->n_used);
      if (particles[scatpar->n_used].ip != 9) break;
      --scatpar->n_used;
      if (scatpar->n_used < 0) break;
//printf("scatpar->n_used = %i...\n", scatpar->n_used);
  }

  return 0;
}
