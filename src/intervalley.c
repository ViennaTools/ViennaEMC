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
/*  Calculation of INTERVALLEY PHONONS scattering rate              */
/*  (absorption + emission)                                         */
/********************************************************************/
int oooRateIntervalley(const_t constpar, scatpar_t *scatpar, intervalley_param_t *intervalley, int *iCount, int *iReg, double *w0, char gORf)
{
  static int i;
  static double Const, ab, ee, ef, gf, em, absorption, emission, rnq;
  FILE *output10;
  FILE *output11;

  /*=== Calculate constants ===*/
  rnq   = 1.0 / (exp(*w0 / constpar.vt) - 1.0);
  Const = intervalley->finalValleys * intervalley->couplingConst * 
    intervalley->couplingConst * constpar.qh * sqrt(constpar.q) 
    / (sqrt(2.0) * M_PI * constpar.density * *w0 * constpar.hbar) * 
    constpar.amd * sqrt(constpar.amd);
  /* MN Note that phonon_GX is the energy and not the frequency of the phonon;
     Thus one hbar less in the denominator in comparison with the usual formulas. */

  /*=== Absorption ===*/
  ++(*iCount);
    
  if (gORf == 'g')
    output10=fopen("rateIntervalleyAbg.csv","w");
  else
    output10=fopen("rateIntervalleyAbf.csv","w");
  
  ab = rnq * Const;
  for (i = 1; i <= NLEV; ++i) 
  {
    ee = i * scatpar->de; 
    ef = ee + *w0 - intervalley->deltaFi;
    gf = ef * (constpar.af * ef + 1.0);

    if (ef <= 0.0) 
      absorption = 0.0;
    else
      absorption = ab * sqrt(gf) * (constpar.af2 * ef + 1.0);

    scatpar->ScatTable[i - 1][*iCount - 1] = absorption;
    fprintf(output10,"%f \t %e \n", ee, absorption); 
  }
  fclose(output10);
  
  scatpar->flagMech[*iCount - 1] = intervalley->i_mech;
  scatpar->w[*iCount - 1] =  *w0 -intervalley->deltaFi;
    
  /*=== Emission ===*/
  ++(*iCount);

  if (gORf == 'g')
    output11=fopen("rateIntervalleyEmg.csv","w");
  else
    output11=fopen("rateIntervalleyEmf.csv","w");

  em = (rnq + 1.0) * Const;
  for (i = 1; i <= NLEV; ++i)
  {
    ee = i * scatpar->de;
    ef = ee - *w0 - intervalley->deltaFi;
    gf = ef * (constpar.af * ef + 1.0);

    if (ef <= 0.0)  emission = 0.0;
    else            emission = em * sqrt(gf) * (constpar.af2 * ef + 1.0);
    
    scatpar->ScatTable[i - 1][*iCount - 1] = emission;
    fprintf(output11,"%f \t %e \n", ee, emission); 
  }
  fclose(output11);

  scatpar->flagMech[*iCount - 1] = intervalley->i_mech;
  scatpar->w[*iCount - 1] =  -(*w0) - intervalley->deltaFi;
  
  return 0;
}
