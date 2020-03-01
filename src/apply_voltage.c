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


/********************************************************************/
/*  APPLY VOLTAGE AT THE SOURCE, DRAIN AND BACK CONTACT             */
/********************************************************************/
int oooApplyVoltage(geometry_t *geometry, phys_quant_t *phys_quantities) {
    int i, j;
    double dopTerm, factor;

    /*=== Source region ===*/
    i = 0;
    for (j = 0; j <= geometry->nymax; ++j) {
        dopTerm = phys_quantities->doping * 0.5;
        factor = dopTerm + sqrt(dopTerm * dopTerm + 1.0);
        phys_quantities->potential[i][j] = log(factor) + phys_quantities->inputVoltage_Vs;
    }

    /*=== Drain region ===*/
    i = geometry->nxmax;
    for (j = 0; j <= geometry->nymax; ++j) {
        dopTerm = phys_quantities->doping * 0.5;
        factor = dopTerm + sqrt(dopTerm * dopTerm + 1.0);
        phys_quantities->potential[i][j] = log(factor) + phys_quantities->inputVoltage_Vd;
    }

    return 0;
}