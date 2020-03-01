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
#include <stdlib.h>
#include <math.h>
#include "emc.h"


/********************************************************************/
/* CALCULATE DISCRETIZATION COEFFICIENTS FOR POISSON EQUATION       */
/* domain1.nr: 1->ohmic contact; 2->semiconductor/oxide interface   */
/*             3->bulk region                    		    */
/* (Dev.independent)						    */
/********************************************************************/
poisson_coef_t *oooPoissonCoefficients(const_t constpar, geometry_t *geometry) {
    int i, j;
    poisson_coef_t *poisson;

    poisson = (poisson_coef_t *) calloc(1, sizeof(poisson_coef_t));

    for (i = 0; i <= geometry->nxmax; ++i)
        for (j = 0; j <= geometry->nymax; ++j)
            if (geometry->domainflag[i][j] == 3)    /* thesis eq. 4.9 */
            {
                /*  A = (yy(j) + yy(j-1)) / xx(i) B = (yy(j) + yy(j-1)) / xx(i-1)  */
                /*  D = (xx(i) + xx(i-1)) / yy(j) E = (xx(i) + xx(i-1)) / yy(j-1)  */
                poisson->coefA[i][j] = (geometry->gridy[j + 1] + geometry->gridy[j]) / geometry->gridx[i + 1];
                poisson->coefB[i][j] = (geometry->gridy[j + 1] + geometry->gridy[j]) / geometry->gridx[i];
                poisson->coefD[i][j] = (geometry->gridx[i + 1] + geometry->gridx[i]) / geometry->gridy[j + 1];
                poisson->coefE[i][j] = (geometry->gridx[i + 1] + geometry->gridx[i]) / geometry->gridy[j];

                if (i == geometry->nxmax) {
                    poisson->coefA[i][j] = 0.0;
                    poisson->coefB[i][j] *= 2.0;
                } else if (i == 0) {
                    poisson->coefA[i][j] *= 2.0;
                    poisson->coefB[i][j] = 0.0;
                }
                if (j == geometry->nymax) {
                    poisson->coefD[i][j] = 0.0;
                    poisson->coefE[i][j] *= 2.0;
                } else if (j == 0) {
                    poisson->coefE[i][j] = 0.0;
                    poisson->coefD[i][j] *= 2.0;
                }

                poisson->coefC[i][j] = -(poisson->coefA[i][j] + poisson->coefB[i][j] +
                                         poisson->coefD[i][j] + poisson->coefE[i][j]);
                poisson->alpha[i][j] = (geometry->gridx[i + 1] + geometry->gridx[i]) * 0.5 *
                                       (geometry->gridy[j + 1] + geometry->gridy[j]);
            } else if (geometry->domainflag[i][j] == 1) {
                poisson->coefA[i][j] = 0.0;
                poisson->coefB[i][j] = 0.0;
                poisson->coefD[i][j] = 0.0;
                poisson->coefE[i][j] = 0.0;
                poisson->coefC[i][j] = 1.0;
            }
    return poisson;
}


/********************************************************************/
/*  INITIALIZE DOPING AND POTENTIAL IN VARIOUS REGIONS:             */
/*  Doping density is renormalized with intrinsic carrier density   */
/*  Potential is renormalized with the thermal voltage              */
/*  (Dev.independent) 						    */
/********************************************************************/
int oooInitializeDopingPotential(const_t constpar, geometry_t *geometry, scatpar_t *scatpar,
                                 phys_quant_t *phys_quantities) {
    static int i, j;//, nRegion;
    static double dopTerm, factor;

    dopTerm = scatpar->dopingCon / constpar.Ni;
    factor = dopTerm + sqrt(dopTerm * dopTerm + 1.0);

    for (i = 0; i <= geometry->nxmax; ++i)
        for (j = 0; j <= geometry->nymax; ++j) {
            phys_quantities->doping = dopTerm;
            phys_quantities->potential[i][j] = log(factor);
        }

    return 0;
}

int oooDistributePotential(const_t constpar, geometry_t *geometry, scatpar_t *scatpar, phys_quant_t *phys_quantities) {
    static int i, j;
    for (i = 0; i <= geometry->nxmax; ++i) {
        for (j = 0; j <= geometry->nymax; ++j) {
//        printf("i = %i\n",i);
//        printf("geometry->nxmax = %i\n",geometry->nxmax);
//        printf("phys_quantities->inputVoltage_Vd = %.8f\n",phys_quantities->inputVoltage_Vd);
//        printf("phys_quantities->inputVoltage_Vs = %.8f\n",phys_quantities->inputVoltage_Vs);
            phys_quantities->potential[i][j] = (phys_quantities->inputVoltage_Vs + i *
                                                                                   (phys_quantities->inputVoltage_Vd -
                                                                                    phys_quantities->inputVoltage_Vs) /
                                                                                   (geometry->nxmax +
                                                                                    1));// / constpar.vt;
        }
//    printf("phys_quantities->potential[%i][all] = %.8f\n",i,phys_quantities->potential[i][0]);
    }
    return 0;
}

/********************************************************************/
/* SOLVE 2D-POISSON EQUATION USING THE SOR-METHOD 		    */
/* flagSolution:  0->equilibrium case; 1->non-equilibrium case      */
/*                2->free memory				    */
/* (Dev.independent)						    */
/********************************************************************/
int oooPoissonSOR(const_t constants, geometry_t *geometry, phys_quant_t *phys_quantities, int flagSolution) {
    static int coefinit = 0;
    static int i, j;
    static double forcingFunc[MAXNX][MAXNY], c2[MAXNX][MAXNY], denn, dennInv,
            dennSum, error, newerror, x1, x2, tempTerm, termA, termB, termD, termE;

    static poisson_coef_t *poisson = NULL;

    /*=== free memory at the end of the program ===*/
    if (flagSolution == 2) {
        free(poisson);
        return 0;
    }

    /*=== Calculate discretization coefficients for Poisson equation (called only once) ===*/
    if (coefinit == 0) {
        poisson = oooPoissonCoefficients(constants, geometry);
        coefinit = 1;
    }

    for (i = 0; i <= geometry->nxmax; ++i)
        for (j = 0; j <= geometry->nymax; ++j)
            if (geometry->domainflag[i][j] == 1) {
                forcingFunc[i][j] = phys_quantities->potential[i][j];
                c2[i][j] = poisson->coefC[i][j];
            }

    do {
        /*=== Calculate the forcing function ===*/
        for (i = 0; i <= geometry->nxmax; ++i)
            for (j = 0; j <= geometry->nymax; ++j) {
                if (geometry->domainflag[i][j] == 1)
                    continue;

                if (flagSolution == 0)
                    denn = exp(phys_quantities->potential[i][j]);
                else
                    denn = phys_quantities->elecDensity[i][j];
                dennInv = exp(-phys_quantities->potential[i][j]);
                dennSum = denn + dennInv;

                c2[i][j] = poisson->coefC[i][j] - poisson->alpha[i][j] * dennSum;
                tempTerm = dennInv - denn + phys_quantities->doping + phys_quantities->potential[i][j] * dennSum;

                if (geometry->domainflag[i][j] == 3)
                    forcingFunc[i][j] = -poisson->alpha[i][j] * tempTerm;
            }

        /*=== Update potential ===*/
        error = 0.0;
        for (i = 0; i <= geometry->nxmax; ++i)
            for (j = 0; j <= geometry->nymax; ++j) {
                if (i == 0) termB = 0.0;
                else termB = poisson->coefB[i][j] * phys_quantities->potential[i - 1][j];
                if (i == geometry->nxmax) termA = 0.0;
                else termA = poisson->coefA[i][j] * phys_quantities->potential[i + 1][j];
                if (j == 0) termE = 0.0;
                else termE = poisson->coefE[i][j] * phys_quantities->potential[i][j - 1];
                if (j == geometry->nymax) termD = 0.0;
                else termD = poisson->coefD[i][j] * phys_quantities->potential[i][j + 1];

                x1 = (forcingFunc[i][j] - termA - termB - termD - termE) /
                     c2[i][j];                        /* thesis eq. 4.11 */
                x2 = phys_quantities->potential[i][j] +
                     constants.conv_omega * (x1 - phys_quantities->potential[i][j]);        /* thesis eq. 4.12 */
                newerror = fabs(x2 - phys_quantities->potential[i][j]);
                if (newerror > error)
                    error = newerror;
                phys_quantities->potential[i][j] = x2;

            }
    } while (error > constants.conv_tolerance);

    return 0;
}