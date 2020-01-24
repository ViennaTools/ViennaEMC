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
#include <stdio.h>


/********************************************************************/
/*  ISOTROPIC SCATTERING PROCESS                                    */
/*  uniform probability density for scattering in all directions    */
/********************************************************************/
particle_t oooIsotropic_g(const_t constpar, scatpar_t *scatpar, particle_t particle, int *iFix)
{
  static double rknew, fi, ct, st;

  particle.e += scatpar->w[*iFix - 1][particle.iRegion - 1];
  /*=== Update carrier wavevector ===*/
  rknew = constpar.smh * sqrt(particle.e * (constpar.af * particle.e + 1.0));

  fi = 2.0 * M_PI * oooRand();
  ct = 1.0 - oooRand() * 2.0;
  st = sqrt(1.0 - ct * ct);
  particle.kx = rknew * st * cos(fi);
  particle.ky = rknew * st * sin(fi);
  particle.kz = rknew * ct;

  return particle;
}


/********************************************************************/
/*  ISOTROPIC SCATTERING PROCESS                                    */
/*     uniform probability density for scattering in all directions */
/********************************************************************/
particle_t oooIsotropic_f(const_t constpar, scatpar_t *scatpar, particle_t particle, int *iFix, int *ivalley_type)
{
  static double rknew, fi, ct, rr, st;
  static int iv0;

  /*=== Update carrier energy ===*/
  particle.e += scatpar->w[*iFix - 1][particle.iRegion - 1];
  /*=== Update carrier wavevector ===*/
  rknew = constpar.smh * sqrt(particle.e * (constpar.af * particle.e + 1.0));

  fi = 2.0 * M_PI * oooRand();
  ct = 1.0 - oooRand() * 2.0;
  st = sqrt(1.0 - ct * ct);
  particle.kx = rknew * st * cos(fi);
  particle.ky = rknew * st * sin(fi);
  particle.kz = rknew * ct;

  rr = oooRand();
  if(*ivalley_type == 2) /* selectMech = 2; f type inter X->X */
  {
      if (particle.iv == 1)
      {
	  if (rr <= 0.5)  iv0 = 2;
	  else            iv0 = 3;
      }
      else if (particle.iv == 2)
      {
	  if (rr <= 0.5)  iv0 = 1;
	  else            iv0 = 3;
      }
      else if (particle.iv == 3)
      {
	  if (rr <= 0.5)  iv0 = 1;
	  else            iv0 = 2;
      }
  }

  particle.iv = iv0;
  return particle;
}


/********************************************************************/
/* COULOMB SCATTERING ANGLE: randomize the polar angle with the     */
/* assumption that the scattering process is elastic                */
/********************************************************************/
particle_t oooCoulombAngleBH(const_t constpar, scatpar_t *scatpar, particle_t particle)
{
  static double kxy, k, ct0, st0, cfi0, sfi0, ge, rr, ct, st, fi, kxp, kyp, kzp;

  /*=== Update carrier energy ===*//* enew = e + w(i_fix,i_region) */
  /*=== Calculate the rotation angles ===*/
  kxy  = sqrt(particle.kx * particle.kx + particle.ky * particle.ky);
  k    = sqrt(kxy * kxy + particle.kz * particle.kz);
  ct0  = particle.kz / k;
  st0  = kxy / k;
  cfi0 = particle.kx / kxy;
  sfi0 = particle.ky / kxy;

  /*=== Randomize momentum in the rotated coordinate system ===*/
  ge  = particle.e * (constpar.af * particle.e + 1.0);
  rr  = oooRand();
  ct  = 1.0 - rr * 2.0 / ((1 - rr) * 4.0 * ge / scatpar->debyeEnergy + 1.0);

  st  = sqrt(1.0 - ct * ct);
  fi  = 2.0 * M_PI * oooRand();
  kxp = k * st * cos(fi);
  kyp = k * st * sin(fi);
  kzp = k * ct;

  /*=== Return to the original coordinate system ===*/
  particle.kx =  kxp * cfi0 * ct0 - kyp * sfi0 + kzp * cfi0 * st0;
  particle.ky =  kxp * sfi0 * ct0 + kyp * cfi0 + kzp * sfi0 * st0;
  particle.kz = -kxp * st0 + kzp * ct0;   /* e = enew */

  return particle;
}

/********************************************************************/
/* COPPER ELECTRON ELECTRON SCATTERING: randomize the direction of  */
/* the momentum and energy                                          */
/********************************************************************/
particle_t oooCopperElecElec(const_t constpar, scatpar_t *scatpar, particle_t particle)
{
  int iRegion;
  static double e, k, rr, ct, tc, kx, ky, kz, st, fai;

  /*=== Set particle energy ===*/
  do { rr = oooRand(); } while (rr <= 1e-6);
   e = -(constpar.vt * 1.5) * log(rr);

//if (iRegion < 1) printf("iRegion = %i...\n", iRegion);
//  k = constpar.smh * sqrt(e * (constpar.af * e + 1.0));
  k = constpar.smh * sqrt(e * (constpar.q * constpar.af * e + 1.0));

  fai = 2.0 * M_PI * oooRand();
  ct  = 1.0 - oooRand() * 2.0;
  st  = sqrt(1.0 - ct * ct);
  kx  = k * st * cos(fai);
  ky  = k * st * sin(fai);
  kz  = k * ct;

  /*=== Map particle atributes ===*/
  particle.kx   = kx;
  particle.ky   = ky;
  particle.kz   = kz;
  particle.e    = e;

  return particle;
}


/********************************************************************/
/* SELECT SCATTERING MECHANISM AND PERFORM THE SCATTERING.  Change  */
/* energy and wavevectors of particles: kx,ky,kz,iv,energy,iRegion. */
/* Select by flagMech:                                              */
/*   1-> Isotropic scattering (acoustic, intervalley g-phonons)     */
/*   2-> Isotropic scattering (intervalley f-phonons)               */
/*   3-> Coulomb scattering ==> small-angle scattering              */
/*   4-> Electron-Electron scattering                               */
/********************************************************************/
particle_t oooScatterCarrier(const_t constpar, scatpar_t *scatpar, particle_t particle)
{
  static int i, iTop, iMech, loc, selectMech;
  static double rr, boundLower, boundUpper;

  /*=== Calculate index to the scattering table ===*/
  loc = (int) (particle.e / scatpar->de);
  if (loc == 0)    loc = 1;
  if (loc > NLEV)  loc = NLEV;

  /*=== Select scattering mechanism ===*/
  iTop = scatpar->maxScatMech[particle.iRegion - 1];
  rr = oooRand();
  if (rr >= scatpar->ScatTable[loc - 1][iTop - 1][particle.iRegion - 1])
      return particle; /* self-scattering */

  if (rr < scatpar->ScatTable[loc - 1][0][particle.iRegion - 1])
    iMech = 1;
  else if (iTop > 1)
    for (i = 1; i <= iTop - 1; ++i)
    {
      boundLower = scatpar->ScatTable[loc - 1][i - 1][particle.iRegion - 1];
      boundUpper = scatpar->ScatTable[loc - 1][i][particle.iRegion - 1];
      if (rr >= boundLower && rr < boundUpper)
      {
	iMech = i + 1;
	break;
      }
    }

  /*=== Perform scattering (change energy and randomize momentum) ===*/
  selectMech = scatpar->flagMech[iMech - 1][particle.iRegion - 1];

  if      (selectMech == 1)  particle = oooIsotropic_g(constpar, scatpar, particle, &iMech);
  else if (selectMech == 2)  particle = oooIsotropic_f(constpar, scatpar, particle, &iMech, &selectMech);
  else if (selectMech == 3)  particle = oooCoulombAngleBH(constpar, scatpar, particle);
  else if (selectMech == 4)  particle = oooCopperElecElec(constpar, scatpar, particle);
  /* else oooIsotropic_f(&iMech, &selectMech); */

  return particle;
}
