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
/* Check boundaries: Artificial boundaries + Source & Drain regions */
/********************************************************************/
particle_t oooCheckBoundary(geometry_t *geometry, phys_quant_t *phys_quantities, particle_t particle)
{
//  static int ix;
  static double deviceLength, deviceDepth;

  deviceLength = geometry->nxmax * geometry->meshSize;
  deviceDepth  = geometry->nymax * geometry->meshSize;

//  ix = (rintf) (particle.xPosition / geometry->meshSize);

  /*=== Check left and right boundary ===*/
  if (particle.xPosition < 0.0)
  {
    particle.iv = 9;
    phys_quantities->Iss_out += particle.weight;
//      printf("Iss_out = %f\n", phys_quantities->Iss_out);
//    particle.xPosition = -particle.xPosition;
//    particle.kx        = -particle.kx;
  } 
  else if (particle.xPosition > deviceLength)
  {
    particle.iv = 9;
    phys_quantities->Idd_out += particle.weight;
//    particle.xPosition =  deviceLength * 2.0 - particle.xPosition;
//    particle.kx        = -particle.kx;
  }

  /*=== Check bottom and top boundary ===*/
  if (particle.yPosition < 0)
  {
    particle.yPosition = -particle.yPosition;
    particle.ky        = -particle.ky;
  }
  else if (particle.yPosition > deviceDepth)
  {
    particle.yPosition =  deviceDepth * 2.0 - particle.yPosition;
    particle.ky        = -particle.ky;
  }

/*
  if (particle.yPosition > deviceDepth)
  {
    particle.yPosition =  deviceDepth * 2.0 - particle.yPosition;
    particle.ky        = -particle.ky;
  }
  else if (particle.yPosition < 0.0)
  {
    if (ix <= geometry->ns)
    {
      particle.iv = 9;
      phys_quantities->Iss_out += particle.weight;
    } 
    else if (ix >= geometry->nd)
    {
      particle.iv = 9;
      phys_quantities->Idd_out += particle.weight;
    } 
    else
    {
      particle.ky = -particle.ky;
      particle.yPosition = -particle.yPosition;
    }
  }
*/
  return particle;
} 


/********************************************************************/
/*  PERFORM THE K-SPACE AND REAL-SPACE MOTION OF THE CARRIERS       */
/********************************************************************/
particle_t oooDrift(const_t constpar, geometry_t *geometry, phys_quant_t *phys_quantities, particle_t particle, double *tau)
{
  static double fx,fy,qh1,dkx,dky, kxf,kyf,kzf,gk, avkx,avky,denom,dx,dy;

  /*=== Calculate Poisson electric field ===*/
//  fx = 1.0/(50e-9)*constpar.vt;//phys_quantities->localfieldx[particle.ixFix][particle.jyFix];
//  fy = 0.0;//phys_quantities->localfieldy[particle.ixFix][particle.jyFix];
  fx = phys_quantities->localfieldx[particle.ixFix][particle.jyFix];
  fy = phys_quantities->localfieldy[particle.ixFix][particle.jyFix];

  /*=== Update particle momentum and energy for X Valleys ===*/
  qh1 = -constpar.qh * *tau;
    
  if (particle.iv == 1)
  {
     dkx = qh1 * fx * constpar.tm[0]; 
     dky = qh1 * fy * constpar.tm[1]; 
  } 
  else if (particle.iv == 2)
  {
     dkx = qh1 * fx * constpar.tm[1]; 
     dky = qh1 * fy * constpar.tm[0]; 
  }
  else if (particle.iv == 3)
  {
     dkx = qh1 * fx * constpar.tm[2]; 
     dky = qh1 * fy * constpar.tm[1]; 
  }
  else /* added to prevent undefined states of dx, dy */
  {
     dkx = 0.0;
     dky = 0.0;
  }

  kxf = particle.kx + dkx;        
  kyf = particle.ky + dky;        
  kzf = particle.kz;              
  gk  = constpar.hhm * (kxf * kxf + kyf * kyf + kzf * kzf); 
  particle.e = gk * 2 / (sqrt(constpar.af4 * gk + 1.0) + 1.0);

  /*=== Update particle's position using the 'leap-frog' scheme ===*/
  avkx = (particle.kx + kxf) * 0.5; 
  avky = (particle.ky + kyf) * 0.5;  
  particle.kx = kxf;
  particle.ky = kyf;
  particle.kz = kzf;

  denom = *tau / (constpar.af2 * particle.e + 1.0);
  if (particle.iv == 1) 
  {
     dx = constpar.hm[0] * avkx * denom;
     dy = constpar.hm[1] * avky * denom;
  } 
  else if (particle.iv == 2)
  {
     dx = constpar.hm[1] * avkx * denom;
     dy = constpar.hm[0] * avky * denom;
  } 
  else if (particle.iv == 3)
  {
     dx = constpar.hm[2] * avkx * denom;
     dy = constpar.hm[1] * avky * denom;
  }
  else /* added to prevent undefined states of dx, dy */
  {
     dx = 0.0;
     dy = 0.0;
  }

  particle.xPosition += dx;
  particle.yPosition += dy;

  /*=== Check simulation boundaries ===*/
  //  printf("drift\n");
  particle = oooCheckBoundary(geometry, phys_quantities, particle);
//  printf("drift2\n");

  if (particle.iv != 9) 
  {
//  printf("drift3\n");
//  printf("particle.xPosition = %f, particle.xPosition = %f, geometry->meshSize = %f\n", particle.xPosition, particle.xPosition, geometry->meshSize);
    particle.ixFix   = (rintf) (particle.xPosition / geometry->meshSize);
//    particle.ixFix   = (int) (particle.xPosition / geometry->meshSize);
//  printf("drift4\n");
    particle.jyFix   = (int) (particle.yPosition / geometry->meshSize);
//  printf("drift5\n");
//  printf("particle.ixFix = %i, particle.jyFix = %i\n", particle.ixFix, particle.jyFix);
    particle.iRegion = geometry->regionflag[particle.ixFix][particle.jyFix];
//  printf("drift6\n");
  }
//  printf("drift7\n");
  return particle;
}


/********************************************************************/
/* Long period (>2×1018) random number generator of L'Ecuyer with   */
/* Bays-Durham shuffle and added safeguards. Returns uniform random */ 
/* deviate between 0.0 and 1.0 (exclusive of the endpoint values).  */
/* Call with iso a negative int to initialize; thereafter, do not   */
/* alter iso between successive deviates in a sequence. RNMX should */ 
/* approximate the largest floating value that is less than 1.      */
/********************************************************************/
double oooRand()
{
  /* Initialized data */
  static int iso  = -1345;
  static int iso2 = 123456789;
  static int ggiv[32] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
  static int iy = 0;
  static int j, k;
  
  /* System generated locals */
  int tmp;
  double ret_val, r2;

  tmp = iso;
  if (tmp <= 0)
  {
    tmp = (-tmp > 1 ? -tmp : 1); /* Computing MAX: iso = max(-iso,1) */
    iso2 = tmp;     
    for (j = 40; j >= 1; --j)
    {
      k    = tmp / 53668;
      tmp = (tmp - k * 53668) * 40014 - k * 12211;

      if (tmp < 0)
	tmp += 2147483563;
      if (j <= 32)
	ggiv[j - 1] = tmp;
    }
    iy = ggiv[0];
  }
  k = tmp / 53668;
  tmp = (tmp - k * 53668) * 40014 - k * 12211;
  if (tmp < 0)
    tmp += 2147483563;
  
  k = iso2 / 52774;
  iso2 = (iso2 - k * 52774) * 40692 - k * 3791;
  if (iso2 < 0)
    iso2 += 2147483399;
  
  j = iy / 67108862 + 1;
  iy = ggiv[j - 1] - iso2;
  ggiv[j - 1] = tmp;

  iso = tmp;
  if (iy < 1)
    iy += 2147483562;
  
  r2 = iy * 4.6566130573917691e-10; /* Computing MIN: rand = min(AM*iy,RNMX) */ 
  ret_val = (r2 < .99999987999999995 ? r2 : .99999987999999995);

  return ret_val;
} 
