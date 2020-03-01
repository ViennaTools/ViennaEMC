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
/*  Material parameter initialization (Si)                          */
/********************************************************************/
int oooMatParInitialization(const_t *constpar) {
    /*=== Define fundamental constants and general constpareters ===*/
    const double am0 = 9.11e-31;        /* am0   = 9.11e-31 [kg] */
    const double eps0 = 8.85419e-12;    /* eps_0 = 8.85419e-12   */
//  const double kb   = 1.38066e-23; 	/* kb    = 1.38066e-23   */
    const double kb = 1.38064852e-23;    /* kb    = 1.38066e-23   */

    /*=== Set temperature and effective masses ===*/
    double TL = 300.0;
    double am_l = 1.0;//0.91;			/* longitudinal eff. mass */
    double am_t = 1.0;//0.19;			/* transversal eff. mass */

//  constpar->q    = 1.60219e-19; 		 /* q = 1.60219e-19 eV */
    constpar->q = 1.60217662e-19;         /* q = 1.60219e-19 eV */
    constpar->hbar = 1.05459e-34;             /* h = 1.05459e-34 Js/rad */
    constpar->qh = constpar->q / constpar->hbar; /* qh = q/h */

    constpar->pi = 3.14159265358979323846264338327;    /* pi */

    constpar->kt = kb * TL;
    constpar->vt = kb * TL / constpar->q;     /* Vt = kB * T/q = 0.02585199 eV*/
    constpar->h = 2. * constpar->hbar * constpar->pi;             /* h = 1.05459e-34 Js/rad */
    printf("Lattice temperature = %f K\n", TL);
    printf("Thermal voltage = %f V\n", constpar->vt);

    /*=== Define material parameters for Si and SiO2 ===*/
    //UPDATE FOR COPPER
    constpar->n = 8.49e28;                          /* electrons per m^3 */
    constpar->hbar2m = constpar->hbar / (2 * am0);
//  constpar->ef          = constpar->h * constpar->h /  ( 8 * am0 ) * pow( 3. *  constpar->n / constpar->pi, 2/3. );
    constpar->ef = constpar->hbar * constpar->hbar / (2 * am0 * constpar->kt)
                   * pow(3. * constpar->n * constpar->pi * constpar->pi, 2 / 3.);
    constpar->phi = 0.;//4.7 / constpar->vt;                              /* Work function of Cu */
    constpar->rho = 1.7E-08;                          /* rho_Cu [Ohm.m] */
    constpar->tau_ee =
            am0 / (constpar->n * constpar->rho * constpar->q * constpar->q);    /* scattering for bulk copper */
    printf("ee = %.5e\n", constpar->tau_ee);
    constpar->eps_sc = eps0;//11.8 * eps0;				/* epsilon_Si */
    constpar->Ni = 2.35e26;//2.137567e24;//1e22;//8.45e22;//1.45e16;                  /* intrinsic carrier density [1/m^3] */
    constpar->density = 8960;//2329.0;  				/* density [kg/m3] */
//  constpar->debyeLength = sqrt(constpar->eps_sc * constpar->vt / (constpar->q * constpar->Ni));
    constpar->debyeLength = sqrt(constpar->eps_sc * constpar->vt / (constpar->q));
    constpar->delta_Ec = 0.575 / constpar->vt;            /* delta_Ec = 0.575/Vt = 22.242 = EG/2 /Vt */

    /*=== Map parameters into internal variables used in the code ===*/
    constpar->amd = pow(am_l * am_t * am_t, 1.0 / 3.0);
    constpar->amc = 3.0 / (1.0 / am_l + 2.0 / am_t);
    constpar->tm[0] = sqrt(constpar->amc / am_l);
    constpar->tm[1] = sqrt(constpar->amc / am_t);
    constpar->tm[2] = sqrt(constpar->amc / am_t);
    constpar->amd = constpar->amd * am0;
    constpar->amc = constpar->amc * am0;
    constpar->hm[0] = constpar->hbar / constpar->amc * constpar->tm[0];
    constpar->hm[1] = constpar->hbar / constpar->amc * constpar->tm[1];
    constpar->hm[2] = constpar->hbar / constpar->amc * constpar->tm[2];
    constpar->af = 0.;//0.5; 						/* nonparabolicity factor */
    constpar->af2 = constpar->af * 2.0;
    constpar->af4 = constpar->af * 4.0;
    constpar->smh = sqrt(constpar->amc * 2.0 * constpar->q) / constpar->hbar;
    constpar->hhm = constpar->hbar / constpar->amc / constpar->qh / 2.0;

    return 0;
}


/********************************************************************/
/*  Scattering parameter initialization                             */
/********************************************************************/
int oooScatParInitialization(scatpar_t *scatpar) {
    scatpar->de = 4.0 / NLEV;

    scatpar->vsound = 9040.0;    /* [m/s]  */
    scatpar->sigma = 6.55;        /* [eV]   */
    scatpar->defpot0g = 5.23e10;    /* [eV/m] */
    scatpar->defpot0f = 5.23e10;    /* [eV/m ] */
    scatpar->defpot1g = 0.0;        /* [eV/m] */
    scatpar->defpot1f = 0.0;        /* [eV/m] */
    scatpar->phonon0g = 63.0e-3;    /* [eV]   */
    scatpar->phonon0f = 59.0e-3;    /* [eV]   */
    scatpar->phonon1g = 27.8e-3;    /* [eV]   */
    scatpar->phonon1f = 29.0e-3;    /* [eV]   */

    /* scattering mechanism flags (1..on, 0..off) */
    scatpar->coulombscattering = 0;  //was 1
    scatpar->acousticscattering = 0; //was 1
    scatpar->intervalley0g = 0;// was 1
    scatpar->intervalley0f = 0;// was 1
    scatpar->intervalley1g = 0;
    scatpar->intervalley1f = 0;

    scatpar->surfaceroughness = 0;
    return 0;
}
