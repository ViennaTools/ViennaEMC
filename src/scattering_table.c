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
#include "emc.h"


/********************************************************************/
/* CREATE THE SCATTERING TABLE FOR EACH REGION                      */
/* Scattering mechanisms: Acoustic phonons, Intervalley g-phonons,  */
/*                        Intervalley f-phonons, Coulomb scattering */
/*      flag_mech = 1 ==> isotropic scattering process              */
/*      flag_mech = 2 ==> polar optical phonons                     */
/********************************************************************/
int oooScatteringTable(const_t constpar, scatpar_t *scatpar)
{
  static int i;
  static double w0;
  static int iCount, iRegion;
  intervalley_param_t intervalley;

  for (iRegion = 1; iRegion <= DOPREG; ++iRegion)
  {
    iCount = 0;   
//    /* Acoustic phonons scattering rate */
//    if (scatpar->acousticscattering == 1)
//      oooRateAcoustic(constpar, scatpar, &iCount, &iRegion);
//
    /* Coulomb scattering: Brooks-Herring approach */
    if (scatpar->coulombscattering == 1)
      oooRateCoulombBH(constpar, scatpar, &iCount, &iRegion);
//
//    /* Intervalley scattering: zero-order interaction (g-process) */
//    if (scatpar->intervalley0g == 1)
//    {
//      w0 = scatpar->phonon0g;
//      intervalley.couplingConst = scatpar->defpot0g;
//      intervalley.deltaFi = 0.0;
//      intervalley.finalValleys = 1.0;
//      intervalley.i_mech = 1;
//      oooRateIntervalley(constpar, scatpar, &intervalley, &iCount, &iRegion, &w0, 'g');
//    }
//    /* Intervalley scattering: zero-order interaction (f-process) */
//    if (scatpar->intervalley0f == 1)
//    {
//      w0 = scatpar->phonon0f;
//      intervalley.couplingConst = scatpar->defpot0f;
//      intervalley.deltaFi = 0.0;
//      intervalley.finalValleys = 4.0;
//      intervalley.i_mech = 2;
//      oooRateIntervalley(constpar, scatpar, &intervalley, &iCount, &iRegion, &w0, 'f');
//    }
//    /* Intervalley scattering: first-order interaction (g-process) */
//    if (scatpar->intervalley1g == 1)
//    {
//      w0 = scatpar->phonon1g;
//      intervalley.couplingConst = scatpar->defpot1g;
//      intervalley.deltaFi = 0.0;
//      intervalley.finalValleys = 1.0;
//      intervalley.i_mech = 1;
//      oooRateIntervalley(constpar, scatpar, &intervalley, &iCount, &iRegion, &w0, 'g');
//    }
//    /* Intervalley scattering: first-order interaction (f-process) */
//    if (scatpar->intervalley1f == 1)
//    {
//      w0 = scatpar->phonon1f;
//      intervalley.couplingConst = scatpar->defpot1f;
//      intervalley.deltaFi = 0.0;
//      intervalley.finalValleys = 4.0;
//      intervalley.i_mech = 2;
//      oooRateIntervalley(constpar, scatpar, &intervalley, &iCount, &iRegion, &w0, 'f');
//    }
//
//
      
//    iCount=1;
//    for (int i = 1; i <= NLEV; ++i)
//      scatpar->ScatTable[i - 1][iCount-1][iRegion - 1] = 1.0/3e-15;
//  
//    scatpar->flagMech[iCount-1][iRegion - 1] = 3;
//    scatpar->w[iCount-1][iRegion - 1] = 0.0;
//
//    scatpar->maxScatMech[iRegion - 1] = iCount;
//    scatpar->taumax[iRegion - 1] = 3e-15;
//    printf("Region = %i  taumax = %e  ", iRegion, scatpar->taumax[iRegion - 1]);

    /* Impact ionization: need to add */
    oooRenormalizeTable(scatpar, &iCount, &iRegion);

    printf("Mechanisms included = %i\n",iCount);
    if (iCount > MAXSC)
      printf("Number of scattering mechanisms exceeds maximum number!\n");
  }

  return 0;
}


/********************************************************************/
/*  Renormalizes the scattering table for X valleys                 */
/********************************************************************/
int oooRenormalizeTable(scatpar_t *scatpar, int *iCount, int *iReg)
{
  static int i, k, imax;
  static double tau;
/*
  FILE *table;
  table = fopen("scattable.csv","w");
*/
  scatpar->maxScatMech[*iReg - 1] = *iCount;
  imax = scatpar->maxScatMech[*iReg - 1];

  if (imax == 0) scatpar->taumax[*iReg - 1] = 2e-15;	/* added for the case that all scattering mechanisms are disabled */
  if (imax >= 1)
  {
   if (imax > 1)
      for (i = 2; i <= imax; ++i)
         for (k = 1; k <= NLEV; ++k)
            scatpar->ScatTable[k - 1][i - 1][*iReg - 1] += scatpar->ScatTable[k - 1][i - 2][*iReg - 1];
   tau = 0.0;
   for (i = 1; i <= NLEV; ++i)
      if (scatpar->ScatTable[i - 1][imax - 1][*iReg - 1] > tau)
         tau = scatpar->ScatTable[i - 1][imax - 1][*iReg - 1];

   for (i = 1; i <= imax; ++i)
      for (k = 1; k <= NLEV; ++k)
      {
         scatpar->ScatTable[k - 1][i - 1][*iReg - 1] /= tau;
      }

   scatpar->taumax[*iReg - 1] = 1.0 / tau;
   printf("Region = %i  taumax = %e  ", *iReg, scatpar->taumax[*iReg - 1]);

    /* fclose(table); */
  }

  return 0;
}
