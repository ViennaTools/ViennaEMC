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
int oooScatteringTable(scatpar_t *scatpar) {
    /* Impact ionization: need to add */
    oooRenormalizeTable(scatpar, 0);
}


/********************************************************************/
/*  Renormalizes the scattering table for X valleys                 */
/********************************************************************/
int oooRenormalizeTable(scatpar_t *scatpar, int iCount) {
    static int i, k, imax;
    static double tau;

    scatpar->maxScatMech = iCount;
    imax = scatpar->maxScatMech;

    if (imax == 0)
        scatpar->taumax = 2e-15;    /* added for the case that all scattering mechanisms are disabled */
    if (imax >= 1) {
        if (imax > 1)
            for (i = 2; i <= imax; ++i)
                for (k = 1; k <= NLEV; ++k)
                    scatpar->ScatTable[k - 1][i - 1] += scatpar->ScatTable[k - 1][i - 2];
        tau = 0.0;
        for (i = 1; i <= NLEV; ++i)
            if (scatpar->ScatTable[i - 1][imax - 1] > tau)
                tau = scatpar->ScatTable[i - 1][imax - 1];

        for (i = 1; i <= imax; ++i)
            for (k = 1; k <= NLEV; ++k) {
                scatpar->ScatTable[k - 1][i - 1] /= tau;
            }

        scatpar->taumax = 1.0 / tau;
        printf("taumax = %e  ", scatpar->taumax);
    }

    return 0;
}
