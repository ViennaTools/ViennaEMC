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
/* DEVICE STRUCTURE INITIALIZATION    		                    */
/*  (1) defines source, drain and gate regions                      */
/*  (2) generates mesh and defines various domains with domainflag: */
/*      1->ohmic contact 2->semiconductor/oxide interface 3->bulk   */
/*  (3) defines doping regions 1->source 2->drain 3->gate 4->bulk   */
/*  (4) calculates coefficients for Poisson equation solution       */
/*  (5) defines initial doping and initial potential                */
/********************************************************************/
int oooDeviceStructureInitialization(const_t constpar, geometry_t *geometry, scatpar_t *scatpar)
{
  static int i, j;
  geometry->nxmax = (int) (geometry->deviceLength / geometry->meshSize);
  geometry->nymax = (int) (geometry->deviceHeight / geometry->meshSize );

  if (geometry->nxmax >= MAXNX)  printf("x-mesh size exceeds maximum size\n\n");
  if (geometry->nymax >= MAXNY)  printf("y-mesh size exceeds maximum size\n\n");
  printf("Mesh size: nxmax = %i  nymax = %i\n", geometry->nxmax, geometry->nymax);

  /*=== Volume calculation ===*/
  geometry->cellVolume = geometry->meshSize * geometry->meshSize * geometry->deviceWidth;

  /*=== Mesh generation (assumes uniform mesh) ===*/
  for (i = 0; i <= geometry->nxmax + 1; ++i)
    geometry->gridx[i] = geometry->meshSize / constpar.debyeLength;
  for (j = 0; j <= geometry->nymax + 1; ++j)
    geometry->gridy[j] = geometry->meshSize / constpar.debyeLength;
  printf("Mesh generated: deviceLength = %i  deviceWidth = %i\n", geometry->nxmax, geometry->nymax);

  /*=== Define domains for Poisson equation solution ===*/
  for (i = 1; i <= geometry->nxmax-1; ++i)
    for (j = 0; j <= geometry->nymax; ++j)             geometry->domainflag[i][j] = 3; /* Bulk */

  /*=== contacts ===*/
  i = 0;
  for (j = 0; j <= geometry->nymax; ++j)                      geometry->domainflag[i][j] = 1; /* Source contact */
  i = geometry->nxmax;
  for (j = 0; j <= geometry->nymax; ++j)		      geometry->domainflag[i][j] = 1; /* Drain contact */

  /*=== Define doping regions ===*/
  for (i = 0; i<= geometry->nxmax; ++i)
    for (j = 0; j <= geometry->nymax; ++j)                 geometry->regionflag[i][j] = 1; /* Doped strip */

  return 0;
}

