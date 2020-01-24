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
/* Initialize carrier energy and wavevector according to the        */
/* Fermi-Dirac statistics                         (Dev.independent) */
/********************************************************************/
int oooInitKspaceFD(const_t constpar, geometry_t *geometry, scatpar_t *scatpar, el_data_t *particles, int *ne, int *i, int *j)
{
  static int iv;
  int iRegion;
  static double e, k, rr, ct, tc, kx, ky, kz, st, fai;

  /*=== Set particle energy ===*/
// === Semiconductor (Si) === //
//  do { rr = oooRand(); } while (rr <= 1e-6);
//   e = -(constpar.vt * 1.5) * log(rr);

// === Copper, based on FD === //
  rr = oooRand();
  e =constpar.ef-log((1.0+exp(0.5*scatpar->kTLevels))/(exp(scatpar->kTLevels*0.5*rr))-1.0);
  
/*  for (int i = 0; i<11; i++){
      double ii = i/10.;
      rr = 1 / (exp(scatpar->kTLevels*(ii-0.5))+1);
      e = (constpar.ef + (rr-0.5)*scatpar->kTLevels);
      rr1 = (i==0)?1e-6:ii;
      e1 = -1.5 * log(rr1);
      printf("ef = %f, nkT = %d, i = %d, e = %f, e1 = %f \n", constpar.ef, scatpar->kTLevels, i, e, e1);
  }*/

  /*=== Set initial valley index and region ===*/
  rr = oooRand() * 3.0;
  if      (rr <= 1.0)  iv = 1;
  else if (rr <= 2.0)  iv = 2;
  else if (rr <= 3.0)  iv = 3;
  iRegion = geometry->regionflag[*i][*j];

//if (iRegion < 1) printf("iRegion = %i...\n", iRegion);
  k = constpar.smh * sqrt(e * (constpar.q * constpar.af * e + 1.0));
  //k = constpar.smh * sqrt(e * (constpar.af * e + 1.0));
  //printf("k = %f\n",k);

  /*=== Initial free-flight ===*/
  do { rr = oooRand(); } while (rr <= 1e-6);
  tc = -log(rr) * scatpar->taumax[iRegion - 1]; 
//  printf("scatpar->taumax[iRegion - 1] = %e \n", scatpar->taumax[iRegion - 1]);
//  printf("tc = %e \n", tc);

  fai = 2.0 * M_PI * oooRand();
  ct  = 1.0 - oooRand() * 2.0;
  st  = sqrt(1.0 - ct * ct);
  kx  = k * st * cos(fai);
  ky  = k * st * sin(fai);
  kz  = k * ct;

  /*=== Check boundaries ===*/
  if (*i == 0 && kx < 0.0)                kx = -kx;
  if (*i == geometry->nxmax && kx > 0.0)  kx = -kx;
  if (*j == 0 && ky < 0.0)            	  ky = -ky;
  if (*j == geometry->nymax && ky > 0.0)  ky = -ky;

  /*=== Map particle atributes ===*/
  particles[*ne].p[0]   = (double)iRegion;
  particles[*ne].p[1]   = kx;
  particles[*ne].p[2]   = ky;
  particles[*ne].p[3]   = kz;
  particles[*ne].p[4]   = tc;
  particles[*ne].p[7]   = 1.0;
  particles[*ne].ip     = iv;
  particles[*ne].energy = e;
//if (particles[*ne].p[0] < 1) printf("particles[%i].p[0] = %e...\n", *ne, particles[*ne].p[0]);
  return 0;
}

/********************************************************************/
/* Initialize carrier energy and wavevector according to the        */
/* Maxwell-Boltzmann statistics                   (Dev.independent) */
/********************************************************************/
int oooInitKspaceMW(const_t constpar, geometry_t *geometry, scatpar_t *scatpar, el_data_t *particles, int *ne, int *i, int *j)
{
  static int iv;
  int iRegion;
  static double e, k, rr, ct, tc, kx, ky, kz, st, fai;
  static double rr1;

  /*=== Set particle energy ===*/
  do { rr = oooRand(); } while (rr <= 1e-6);
   e = -(constpar.vt * 1.5) * log(rr);

//   do { rr1 = oooRand(); } while (rr <= 1e-6);
//   if (rr < 0.5) e = constpar.ef-e;
//   else e = constpar.ef+e;
   
  /*=== Set initial valley index and region ===*/
  rr = oooRand() * 3.0;
  if      (rr <= 1.0)  iv = 1;
  else if (rr <= 2.0)  iv = 2;
  else if (rr <= 3.0)  iv = 3;
  iRegion = geometry->regionflag[*i][*j];

//if (iRegion < 1) printf("iRegion = %i...\n", iRegion);
  k = constpar.smh * sqrt(e * (constpar.q * constpar.af * e + 1.0));
  //k = constpar.smh * sqrt(e * (constpar.af * e + 1.0));
  //printf("k = %f\n",k);

  /*=== Initial free-flight ===*/
  do { rr = oooRand(); } while (rr <= 1e-6);
  tc = -log(rr) * scatpar->taumax[iRegion - 1]; 
//  printf("scatpar->taumax[iRegion - 1] = %e \n", scatpar->taumax[iRegion - 1]);
//  printf("tc = %e \n", tc);

  fai = 2.0 * M_PI * oooRand();
  ct  = 1.0 - oooRand() * 2.0;
  st  = sqrt(1.0 - ct * ct);
  kx  = k * st * cos(fai);
  ky  = k * st * sin(fai);
  kz  = k * ct;

  /*=== Check boundaries ===*/
  if (*i == 0 && kx < 0.0)                kx = -kx;
  if (*i == geometry->nxmax && kx > 0.0)  kx = -kx;
  if (*j == 0 && ky < 0.0)            	  ky = -ky;
  if (*j == geometry->nymax && ky > 0.0)  ky = -ky;

  /*=== Map particle atributes ===*/
  particles[*ne].p[0]   = (double)iRegion;
  particles[*ne].p[1]   = kx;
  particles[*ne].p[2]   = ky;
  particles[*ne].p[3]   = kz;
  particles[*ne].p[4]   = tc;
  particles[*ne].p[7]   = 1.0;
  particles[*ne].ip     = iv;
  particles[*ne].energy = e;
//if (particles[*ne].p[0] < 1) printf("particles[%i].p[0] = %e...\n", *ne, particles[*ne].p[0]);
  return 0;
}


/********************************************************************/
/* Initialize carrier position on node location   (Dev.independent) */
/********************************************************************/
int oooInitRealspace(geometry_t *geometry, el_data_t *particles, int *ne, int *i, int *j)
{
  static double xPosition, yPosition;

  /*=== Define particle coordinates based on grid location ===*/
  if (*i != 0 && *i != geometry->nxmax)
    xPosition = *i + oooRand() - 0.5;
  else if (*i == 0)
    xPosition = oooRand() * 0.5;
  else if (*i == geometry->nxmax)
    xPosition = geometry->nxmax - oooRand() * 0.5;
  
  if (xPosition < 0.0)              xPosition = -xPosition;
  if (xPosition > geometry->nxmax)  xPosition = geometry->nxmax * 2.0 - xPosition;

  if (*j != 0 && *j != geometry->nymax)
    yPosition = *j + oooRand() - 0.5;
  else if (*j == 0)
      yPosition = oooRand() * 0.5;
  else if (*j == geometry->nymax)
    yPosition = geometry->nymax - oooRand() * 0.5;
  
  if (yPosition < 0.0)              yPosition = -yPosition;
  if (yPosition > geometry->nymax)  yPosition = geometry->nymax * 2.0 - yPosition;

  /*=== Map particle atributes ===*/
  particles[*ne].p[5] = xPosition * geometry->meshSize; 
  particles[*ne].p[6] = yPosition * geometry->meshSize; 
  
  return 0;
}