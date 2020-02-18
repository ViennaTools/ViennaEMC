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
#include <time.h>
#include "emc.h"


/*=== core Ensemble Monte Carlo code ===*/
output_t *EMC(const_t constants, geometry_t *geometry, scatpar_t *scatpar, phys_quant_t *phys_quantities) {
    int i, j;
    int iterRefPrint, iterTotal, nTimes, nTimeStepsAv, nTimeStepsAvPrint;
    double sourceFactor, drainFactor, sourceSum, drainSum;
    double Time, avgtimesteps;
    double CFACT;

    el_data_t *particles;
    currents_t *current_cumulative;
    curr_from_charge_t *current_from_charge;

    time_t tBegin, tEnd;

    /*=== allocate memory for particle data ===*/
    particles = (el_data_t *) calloc((MAXEN + 1), sizeof(el_data_t));
    Time = 0.0;
    sourceSum = 0.0;
    drainSum = 0.0;

    /*=== specify doping and initial potential in various regions ===*/
    oooInitializeDopingPotential(constants, geometry, scatpar, phys_quantities);
    /*=== calculate the scattering table ===*/
    oooScatteringTable(constants, scatpar);
    printf("===========================================================\n");

    CFACT = constants.q / scatpar->dt;

    /*=== solve equilibrium Poisson equation for electron initialization ===*/
    phys_quantities->inputVoltage_Vs = 0.0;
    phys_quantities->inputVoltage_Vd = 0.0;//1 / constants.vt;
    constants.conv_tolerance /= constants.vt;

    /*=== initialize potential again and write it out ===*/
    oooApplyVoltage(constants, geometry, phys_quantities);
    oooPoissonSOR(constants, geometry, phys_quantities, 0);
    oooElectricFieldUpdate(constants, geometry, scatpar, phys_quantities, &Time);
    oooWritePotential(constants, geometry, phys_quantities, 1, 0);

    /*=== initialize time averaged potential and density and averaging time ===*/
    avgtimesteps = 1e-12 / scatpar->dt;
    oooAveragePotential(geometry, phys_quantities, 0, 0.0);
    oooAverageDensity(geometry, phys_quantities, 0, 0.0);

    /*=== initialize electrons, count them, and write them out ===*/
    oooElectronsInitialization(constants, geometry, scatpar, particles, phys_quantities);

    /*=== apply bias at the contacts ===*/
    phys_quantities->inputVoltage_Vs = phys_quantities->inputVs / constants.vt;
    phys_quantities->inputVoltage_Vd = phys_quantities->inputVd / constants.vt;
    oooApplyVoltage(constants, geometry, phys_quantities);

    /*=== start the Boltzmann Monte Carlo procedure ===*/
    /*=== synchronize the times ===*/
    tBegin = time((time_t *) 0);

    scatpar->averTime = (int) (scatpar->averTime / scatpar->dt + 0.5) * scatpar->dt;
    scatpar->transientTime = (int) (scatpar->transientTime / scatpar->averTime + 0.5) * scatpar->averTime;
    scatpar->totalTime = (int) (scatpar->totalTime / scatpar->averTime + 0.5) * scatpar->averTime;
    nTimeStepsAv = (int) ((scatpar->totalTime - scatpar->transientTime) / scatpar->dt + 0.5);
    nTimeStepsAvPrint = (int) (scatpar->averTime / scatpar->dt + 0.5);
    nTimes = (int) (scatpar->totalTime / scatpar->dt + 0.5) / nTimeStepsAvPrint;
    iterTotal = nTimeStepsAvPrint * nTimes;

    printf("\nNumber of time steps for printing  = %i\n", nTimeStepsAvPrint);
    printf("Number of time steps for averaging = %i\n", nTimeStepsAv);
    printf("Total number of iterations = %i\n", iterTotal);
    printf("===========================================================\n");

    /*=== allocate memory for current data ===*/
    current_cumulative = calloc(iterTotal, sizeof(currents_t));
    current_from_charge = calloc(iterTotal, sizeof(curr_from_charge_t));

    printf("STARTING MONTE CARLO PROCEDURE...\n");

    for (j = 1, i = 0; j <= iterTotal; ++j) {
        Time = j * scatpar->dt;
        iterRefPrint = j % (nTimeStepsAvPrint);

        oooFreeFlightScatter(constants, geometry, scatpar, particles, phys_quantities);
        oooCheckSourceDrainContacts(constants, geometry, scatpar, particles, phys_quantities);
        oooDeleteParticles(scatpar, particles);
        oooChargeAssignmentNEC(constants, geometry, scatpar, particles, phys_quantities);
        printf("iter = %d\n", j);
        oooPoissonSOR(constants, geometry, phys_quantities, 1);
        oooElectricFieldUpdate(constants, geometry, scatpar, phys_quantities, &Time);

        /*=== calculate cumulative velocities, energies & currents; write cumulative and momentary currents into file ===*/
        current_cumulative[j - 1] = oooVelocityEnergyCumulative(constants, geometry, scatpar, particles,
                                                                phys_quantities, &Time, &nTimeStepsAv);
        current_from_charge[j - 1].Time = Time * 1e12;

        /*=== sum up the potential & density during the last 1ps for average calculation ===*/
        if (Time >= (scatpar->totalTime - 1e-12)) {
            oooAveragePotential(geometry, phys_quantities, 1, 0.0);
            oooAverageDensity(geometry, phys_quantities, 1, 0.0);
        }

        if (Time >= scatpar->transientTime) {
            drainFactor = phys_quantities->Idd_out + phys_quantities->Idd_eli - phys_quantities->Idd_cre;
            sourceFactor = phys_quantities->Iss_out + phys_quantities->Iss_eli - phys_quantities->Iss_cre;
            ++i;
            drainSum += drainFactor;
            sourceSum += sourceFactor;

            current_from_charge[j - 1].sourceCurr = -sourceSum / i * CFACT;
            current_from_charge[j - 1].drainCurr = drainSum / i * CFACT;
            current_from_charge[j - 1].sourceFactor = sourceFactor;
            current_from_charge[j - 1].drainFactor = drainFactor;

            /*=== just for plotting purpose at every averaging time ===*/
            if (iterRefPrint == 0) {
                time_t tCur = time((time_t *) 0);
                printf("===========================================================\n");
                printf("Time = %.2e \t\t\t CPU time = %ld s\n", Time, tCur - tBegin);
                printf("===========================================================\n");
                printf("        out  eli  cre  sum\t Current\n");
                printf("Source  %.f + %.f - %.f = %.f\t %e   %f\n",
                       phys_quantities->Iss_out, phys_quantities->Iss_eli, phys_quantities->Iss_cre,
                       sourceFactor, -sourceSum / i * CFACT, current_cumulative[j - 1].Is_cumul);
                printf("Drain   %.f + %.f - %.f = %.f\t %e   %f\n",
                       phys_quantities->Idd_out, phys_quantities->Idd_eli, phys_quantities->Idd_cre,
                       drainFactor, drainSum / i * CFACT, current_cumulative[j - 1].Id_cumul);
                printf("===========================================================\n");
                printf("Number of electrons in use = %i\n", scatpar->n_used);
                oooWritePotential(constants, geometry, phys_quantities, 1, j);
            }
        }
    }

    /*=== calculate and print total simulation time ===*/
    tEnd = time((time_t *) 0);
    printf("===========================================================\n");
    printf("SIMULATION COMPLETED in %ld s\n", tEnd - tBegin);
    printf("===========================================================\n");

    /*=== calculate averaged potential and density ===*/
    oooAveragePotential(geometry, phys_quantities, 2, avgtimesteps);
    oooAverageDensity(geometry, phys_quantities, 2, avgtimesteps);

    free(particles);

    return oooCalculateOutput(constants, geometry, phys_quantities, nTimeStepsAv, (int) (tEnd - tBegin), iterTotal,
                              current_cumulative, current_from_charge);
}
