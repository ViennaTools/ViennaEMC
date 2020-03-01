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
/* Calculation of COULOMB SCATTERING rate (Brooks-Herring approach) */
/* Assumption ==> elastic scattering process                        */
/********************************************************************/
int oooRateCoulombBH(const_t constpar, scatpar_t *scatpar, int *iCount)
{
  static int i;
  static double debyeLengthSq, Const, ee, ge, factor, scatRate;
  FILE *output10;

  /*=== Calculate constants ===*/
//  debyeLengthSq = constpar.eps_sc * constpar.vt / constpar.q / fabs(scatpar->dopingCon);
//  scatpar->debyeEnergy = constpar.hhm / debyeLengthSq;
//  printf("scatpar->debyeEnergy = %e\n",scatpar->debyeEnergy);

//  factor = debyeLengthSq * constpar.qh * constpar.qh / constpar.eps_sc;
//  Const  = fabs(scatpar->dopingCon) * constpar.amd / M_PI *
//           sqrt(constpar.amd * 2.0 * constpar.q) * factor * factor;
// 
//  output10 = fopen("rateCoulomb.csv","w");

  /*=== Calculate scattering rate ===*/
  ++(*iCount);
  for (i = 1; i <= NLEV; ++i)
  {
//    ee = i * scatpar->de;
//    ge = ee * (constpar.af * ee + 1.0);
//    
//    factor = sqrt(ge) * (constpar.af2 * ee + 1.0) / (ge * 4.0 / scatpar->debyeEnergy + 1.0);
//    scatRate = Const * factor;

//    scatpar->ScatTable[i - 1][*iCount - 1][*iReg - 1] = scatRate;
    scatpar->ScatTable[i - 1][*iCount - 1] = 1.0/constpar.tau_ee;
//    fprintf(output10,"%f \t %e \n", ee, scatRate); 
  }
//  fclose(output10);
  
//  scatpar->flagMech[*iCount - 1][*iReg - 1] = 3;
  scatpar->flagMech[*iCount - 1] = 4;
  scatpar->w[*iCount - 1] = 0.0;

  return 0;
}