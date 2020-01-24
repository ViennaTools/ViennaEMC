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
#define min(X, Y) (((X) < (Y)) ? (X) : (Y))
#define max(X, Y) (((X) > (Y)) ? (X) : (Y))

#ifndef M_PI
  #define M_PI 3.14159265358979323846264338327
#endif

#ifndef VIENNAWD_ENSEMBLEMC_EMC_H
#define VIENNAWD_ENSEMBLEMC_EMC_H

#define MAXNX 1101	/* max. mesh dimensions */
#define MAXNY 201

#define DOPREG 4	/* max. number of doping regions */
#define MAXSC 10	/* max. number of scattering mechanisms */
#define NLEV 1000	/* number of energy levels in the scattering table */
#define MAXEN 50000000	/* max. number of particles */  //was 5000000


/*=== data structure holds the number of carriers in the source and drain contacts; returned by oooSDcarrierNumber(...) ===*/
typedef struct
{
  int nScarriers[MAXNX], nDcarriers[MAXNX];
} SDcarriers_t;


/*=== data structure for intervalley scattering parameters ===*/
typedef struct
{
  double couplingConst, deltaFi, finalValleys;
  int i_mech;
} intervalley_param_t;


/*=== data structure for poisson discretization coefficients; returned by oooPoissonCoefficients(...) ===*/
typedef struct
{
  double coefA[MAXNX][MAXNY],
         coefB[MAXNX][MAXNY],
         coefC[MAXNX][MAXNY],
	 coefD[MAXNX][MAXNY],
         coefE[MAXNX][MAXNY],
         alpha[MAXNX][MAXNY];
} poisson_coef_t;


/*=== data structure for current output; returned by oooVelocityEnergyCumulative(...) */
typedef struct
{
  double Is_cumul, Id_cumul;			/* cumulative source & drain currents */
  double Is_momentary, Id_momentary;		/* momentary source & drain currents */
} currents_t;


typedef struct
{
  double Time;
  double sourceCurr, drainCurr;
  double sourceFactor, drainFactor;
} curr_from_charge_t;


/*=== data structure for constants and parameters ===*/
typedef struct
{
  double q, hbar, h, qh;				/* elementary charge and hbar */
  double amd, amc, smh, hhm, tm[3], hm[3];
  double af, af2, af4;				/* nonparabolicity factors */
  double delta_Ec;				/* delta_Ec = 0.575/Vt = 22.242 = EG/2 /Vt */

  double pi;
  double kt;
  double vt; 					/* thermal voltage */
  double debyeLength; 				/* debye length */

  double Ni; 					/* intrinsic density for Si */
  double density; 				/* Si density */
  double eps_sc, eps_ox, gamma; 		/* dielectric constants */
  
  double rho, tau_ee;                           /* added to deal with calculating copper bulk resistivity */
  double phi, hbar2m;
  double ef, n;                                     /* free electron density*/
  double conv_omega, conv_tolerance;		/* convergence parameters for the poisson solver */
} const_t;


/*=== data structure for device geometry ===*/
typedef struct
{
//  double lengthSD, depthSD, lengthG, depthB;    /* MOS structure data in [m] */
  double deviceLength;
  double deviceWidth;
  double deviceHeight;
//  double oxideThickness;
  double cellVolume;

  double meshSize;				/* mesh spacing (in both directions equal) in [m] */
  int nxmax, nymax;				/* maximum mesh dimensions */
//  int ns, nd, ny;				/* ns ... x-index where source region ends; nd ... x-index where drain region starts; ny ... y-index where contact region doping end */

  int regionflag[MAXNX][MAXNY];			/* doping region flags */
  int domainflag[MAXNX][MAXNY];			/* region flags for the poisson solver */

  double gridx[MAXNX+1];			/* mesh grid = meshSize / debyeLength */
  double gridy[MAXNY+1];
} geometry_t;


/*=== data structure for scattering parameters ===*/
typedef struct
{
  double sigma, vsound;					/* acoustic scattering constants */
//  double dopingCon[DOPREG], debyeEnergy[DOPREG];	/* coulomb scattering constants */
  double dopingCon, debyeEnergy;

  double tau_Cu;                                        /* Scattering time for Cu */
  
  /*=== intervalley scattering parameters ===*/
  double defpot0g, defpot0f, defpot1g, defpot1f;
  double phonon0g, phonon0f, phonon1g, phonon1f;

  /*=== scattering table parameters ===*/
  double de, w[MAXSC][DOPREG], ScatTable[NLEV][MAXSC][DOPREG];
  int maxScatMech[DOPREG], flagMech[MAXSC][DOPREG];

  /*=== scattering mechanism selection flags ===*/
  int acousticscattering, coulombscattering;
  int intervalley0g, intervalley0f, intervalley1g, intervalley1f;

  double dt, totalTime, averTime, transientTime;	/* timestep, total simulation time, averaging time and transient time */
  double taumax[DOPREG];

  int surfaceroughness;
  int n_used;						/* number of particles currently in use */
  
  /*=== distributing particle energy with F-D ===*/
  int kTLevels;
  int precision;
  double *energy_levels;
  double *cumulative;
} scatpar_t;


/*=== data structure for handling a single particle ===*/
typedef struct
{
  double kx, ky, kz,
         xPosition, yPosition,
         weight,
         e;
  int    iv, iRegion,
         ixFix, jyFix;
} particle_t;


/*=== data structure for the particle array ===*/
typedef struct
{
  double p[8];
  double energy;
  int ip;
} el_data_t;


/*=== data structure for physical simulation quantities ===*/
typedef struct
{
  double potential[MAXNX][MAXNY];
  double doping[MAXNX][MAXNY];
  double elecDensity[MAXNX][MAXNY];

  double fxField[MAXNX][MAXNY], fyField[MAXNX][MAXNY];		/* electric field */
  double localfieldx[MAXNX][MAXNY], localfieldy[MAXNX][MAXNY];  /* local electric field */

  double fieldX_avg[MAXNX], fieldY_avg[MAXNX], 			/* average electric field cuts along x-direction */
         fieldX_avg2[MAXNX], fieldY_avg2[MAXNX];
  double densityX_avg[MAXNX];					/* average density cut along x-direction */

  double biasbeta[MAXNX];

  double Idd_out, Idd_eli, Idd_cre;				/* Drain currents  */
  double Iss_out, Iss_eli, Iss_cre;				/* Source currents */

  double velxSum[MAXNX], velySum[MAXNX], 			/* average x,y velocities along x-direction */
         enerSum[MAXNX], 					/* average energy along x-direction */
         curSumX[MAXNX][MAXNY], curSumY[MAXNX][MAXNY], 		/* average currents x,y velocities along x-direction */
         elSum[MAXNX][MAXNY]; 					/* electron count at (x,y) */

//  double inputVs, inputVd, inputVg, inputVsub;			/* Input voltages  */
  double inputVs, inputVd;			/* Input voltages  */
  double inputVoltage_Vs,					/* Input voltages normalized by V_t */
         inputVoltage_Vd;
//         inputVoltage_Vg,
//         inputVoltage_Vsub;

} phys_quant_t;


/*=== data structure for output ===*/
typedef struct
{
  int totaltime;
  int x_max, y_max;
  int iterTotal;
  double x_axis[MAXNX],
         y_axis[MAXNY];

  currents_t *current_cumulative;
  curr_from_charge_t *current_from_charge;

  /* 2D quantities */
  double potential[MAXNX][MAXNY];
  double electrondensity[MAXNX][MAXNY];

  double fieldXY_x[MAXNX][MAXNY],
         fieldXY_y[MAXNX][MAXNY];

  double currentdensityX[MAXNX][MAXNY],
   	 currentdensityY[MAXNX][MAXNY],
  	 currentdensity[MAXNX][MAXNY];

  /* quantities along x direction with y = 0 */
  double fieldX[MAXNX], fieldX2[MAXNX],
         fieldY[MAXNX], fieldY2[MAXNY];
  double density[MAXNX];
  double sheetdensity[MAXNX];

  double velocityX[MAXNX],
         velocityY[MAXNX];
  double energy[MAXNX];
} output_t;


/*=== Initialization functions ===*/
int oooMatParInitialization(const_t *);
int oooScatParInitialization(scatpar_t *);
int oooDeviceStructureInitialization(const_t, geometry_t *, scatpar_t *);
int oooInitializeDopingPotential(const_t, geometry_t *, scatpar_t *, phys_quant_t *);
int oooElectronsInitialization(const_t, geometry_t *, scatpar_t *, el_data_t *, phys_quant_t *);
int oooInitKspaceMW(const_t, geometry_t *, scatpar_t *, el_data_t *, int *, int *, int *);
int oooInitKspaceFD(const_t, geometry_t *, scatpar_t *, el_data_t *, int *, int *, int *);
int oooInitRealspace(geometry_t *, el_data_t *, int *, int *, int *);
SDcarriers_t oooSDcarrierNumber(const_t, geometry_t *, phys_quant_t *);

double oooFunction(double, const_t);
double oooIntegrateLeftrect(double, double, double, const_t, double (*)());
double oooIntegrateRightrect(double, double, double, const_t, double (*)());
double oooIntegrateMidrect(double, double, double, const_t, double (*)());
double oooIntegrateTrapezium(double, double, double, const_t, double (*)());
double oooIntegrateSimpson(double, double, double, const_t, double (*)());
double oooIntegrate(const_t, double);

/*=== core functions ===*/
output_t *EMC(const_t, geometry_t *, scatpar_t *, phys_quant_t *);
double oooRand(void);
int oooApplyVoltage(const_t, geometry_t *, phys_quant_t *);
int oooPoissonSOR(const_t, geometry_t *, phys_quant_t *, int);
int oooDistributePotential(const_t, geometry_t *, scatpar_t *, phys_quant_t *);
int oooElectricFieldUpdate(const_t, geometry_t *, scatpar_t *, phys_quant_t *, double *);
int oooCheckSourceDrainContacts(const_t, geometry_t *, scatpar_t *, el_data_t *, phys_quant_t *);
int oooDeleteParticles(scatpar_t *, el_data_t *);
int oooChargeAssignmentNEC(const_t, geometry_t *, scatpar_t *, el_data_t *, phys_quant_t *);

/*=== Drift and scattering functions ===*/
particle_t oooDrift(const_t, geometry_t *, phys_quant_t *, particle_t, double *);
int oooScatteringTable(const_t, scatpar_t *);
particle_t oooScatterCarrier(const_t, scatpar_t *, particle_t);
int oooFreeFlightScatter(const_t, geometry_t *, scatpar_t *, el_data_t *, phys_quant_t *);
int oooRateAcoustic(const_t, scatpar_t *, int *, int *);
int oooRateCoulombBH(const_t, scatpar_t *, int *, int *);
int oooRateIntervalley(const_t, scatpar_t *, intervalley_param_t *, int *, int *, double *, char);
int oooRenormalizeTable(scatpar_t *, int *, int *);
double oooRateRelaxation();

int oooAveragePotential(geometry_t *, phys_quant_t *, int, double);
int oooAverageDensity(geometry_t *, phys_quant_t *, int, double);

currents_t oooVelocityEnergyCumulative(const_t, geometry_t *, scatpar_t *, el_data_t *, phys_quant_t *, double *, int *);

/*=== Output functions ===*/
output_t *oooCalculateOutput(const_t, geometry_t *, phys_quant_t *, int, int, int, currents_t *, curr_from_charge_t *);
//output_t *oooCalculateOutput(const_t, geometry_t *, phys_quant_t *, int, int, int);
//int oooWriteOutput(output_t *);
int oooWriteOutput(output_t *, geometry_t *, phys_quant_t *);
//int oooWriteElectricField(const_t, geometry_t *, phys_quant_t *, int *);
//int oooWriteVelocityEnergy(geometry_t *, phys_quant_t *, int *);
//int oooWriteCurrentDensity(geometry_t *, phys_quant_t *, int *);
//int oooWriteMesh(const_t, geometry_t *);
//int oooWritePotential(const_t, geometry_t *, phys_quant_t *, int);
int oooWritePotential(const_t, geometry_t *, phys_quant_t *, int, int);
//int oooWriteCurrents(int, currents_t *, curr_from_charge_t *);

#endif /*VIENNAWD_ENSEMBLEMC_EMC_H*/
