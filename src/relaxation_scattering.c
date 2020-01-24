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
double oooRateRelaxation()
{
  static double tau = 1e-16;
  return 1.0/tau;
}
