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
#include <omp.h>


int oooMapAttributes(el_data_t *particles, particle_t particle, int *i, double *dtau)
{
  /*=== Map particle atributes ===*/
  particles[*i].p[0]   = (double) particle.iRegion; 
  particles[*i].p[1]   = particle.kx; 
  particles[*i].p[2]   = particle.ky; 		
  particles[*i].p[3]   = particle.kz; 		
  particles[*i].p[4]   = *dtau;       		
  particles[*i].p[5]   = particle.xPosition;
  particles[*i].p[6]   = particle.yPosition;

  particles[*i].ip     = particle.iv;
  particles[*i].energy = particle.e;
  
  return 0;
}


/********************************************************************/
/*  PERFORM THE FREE-FLIGHT & SCATTER PART WITHIN ONE TIME INTERVAL */
/********************************************************************/
int oooFreeFlightScatter(const_t constpar, geometry_t *geometry, scatpar_t *scatpar, el_data_t *particles, phys_quant_t *phys_quantities) 
{

  static int i;
  static double dtau, rr, dt2, dt3, dtp;
  particle_t particle;

  /*=== Reset the electron number ===*/
  phys_quantities->Idd_out = 0.0;  phys_quantities->Iss_out = 0.0;
  phys_quantities->Idd_eli = 0.0;  phys_quantities->Iss_eli = 0.0;
  phys_quantities->Idd_cre = 0.0;  phys_quantities->Iss_cre = 0.0;

  /*=== Loop over all used carriers ===*/
//  #pragma omp parallel
//  #pragma omp for
  #pragma omp for schedule(static, 1)
  for (i = 0; i <= scatpar->n_used; ++i)
  {

    /*=== Inverse mapping of particle atributes ===*/
    particle.iRegion   = (rintf)(particles[i].p[0]); 
    particle.kx        = particles[i].p[1]; 
    particle.ky        = particles[i].p[2]; 
    particle.kz        = particles[i].p[3]; 
    dtau               = particles[i].p[4];
    particle.xPosition = particles[i].p[5];
    particle.yPosition = particles[i].p[6];
    particle.weight    = particles[i].p[7];
    particle.iv        = particles[i].ip;
    particle.e         = particles[i].energy;
//    particle.ixFix     = (int) (particle.xPosition / geometry->meshSize);
    particle.ixFix     = (rintf) (particle.xPosition / geometry->meshSize);
    particle.jyFix     = (int) (particle.yPosition / geometry->meshSize);

    /*=== Initial free-flight of the carriers ===*/
    dt2 = (dtau < scatpar->dt ? dtau : scatpar->dt); 				/* min */
//    printf("dt2 = %.e\n", dt2);
    particle = oooDrift(constpar, geometry, phys_quantities, particle, &dt2); 	/* time to advance the trajectory */
    if (particle.iv == 9)
    {
      oooMapAttributes(particles, particle, &i, &dtau);  
      continue;
    }

    /*=== Free-flight and scatter part ===*/
    if (dtau <= scatpar->dt)
      do
      {
	if (scatpar->maxScatMech > 0)
            particle = oooScatterCarrier(constpar, scatpar, particle); //scatter takes place
	do { rr = oooRand(); } while (rr <= 1e-6);
	dt3 = -log(rr) * scatpar->taumax;
	dtp = scatpar->dt - dtau; 	/* remaining time to scatter in dt-interval */
	dt2 = (dt3 <= dtp ? dt3 : dtp); /* min */
	particle = oooDrift(constpar, geometry, phys_quantities, particle, &dt2);
	if (particle.iv == 9)
	{
	  oooMapAttributes(particles, particle, &i, &dtau);  
	  break;
	}

	dtau += dt3; /* update times */
      } while (dtau < scatpar->dt);
    dtau -= scatpar->dt; 

    oooMapAttributes(particles, particle, &i, &dtau);
  }

  return 0;
}
