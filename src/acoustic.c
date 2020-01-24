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
/*  Calculation of acoustic phonon scattering rate                  */
/********************************************************************/
/*=== Acoustic scattering rate for the X Valleys ===*/
int oooRateAcoustic(const_t constpar, scatpar_t *scatpar, int *iCount, int *iReg)
{
  static int i;
  static double scatRate, Const, ee, fe, cL;
  FILE *output10;

  /* Calculate constants */
  cL = constpar.density * scatpar->vsound * scatpar->vsound;

  Const = sqrt(2.0 * constpar.amd * constpar.q) * constpar.amd * constpar.vt * constpar.qh *
          constpar.qh * constpar.qh / (M_PI * cL * constpar.hbar) * 
          scatpar->sigma * scatpar->sigma;

  output10 = fopen("rateAcoustic.csv", "w");

  /* Create scattering table */
  ++(*iCount);
  for (i = 1; i <= NLEV; ++i)
  {
    ee = i * scatpar->de;
    fe = ee * (constpar.af * ee + 1.0);
    scatRate = Const * sqrt(fe) * (constpar.af2 * ee + 1.0);
    scatpar->ScatTable[i - 1][*iCount - 1][*iReg - 1] = scatRate;
    fprintf(output10, "%f \t %e \n", ee, scatRate); 
  }
  fclose(output10);

  scatpar->flagMech[*iCount - 1][*iReg - 1] = 1;
  scatpar->w[*iCount - 1][*iReg - 1] = 0.0;

  return 0;
}
