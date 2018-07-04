/*** Allard Lab jun.allard@uci.edu                    ***/

#define TWISTER genrand_real3()
#define NFILMAX         10
#define NMAX            400
#define NTMAX           2e9
#define NTADAPT         20000
#define NTCHECK         200000
#define DCHIMIN         1e-4
#define NBINS           100
#define NBINSPOLYMER    3000
#define PI              3.14159265359
#define INF             1e14
#define DCHIINIT        0.1
#define KSCRITICAL      0.01
#define MEMBRANE        1
#define MULTIPLE        0
#define STIFFEN         0
#define ELECTRO         1
#define HARDWALL        0
#define BASEBOUND       0
#define CPMAX           1e8
#define TALKATIVE       1
#define TXTPARAM        1
#define VISUALIZE       0
#define CD3ZETA         1
#define BINDTRANSITION  0

#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include "twister.c"


/*******************************************************************************/
//  GLOBAL VARIABLES
/*******************************************************************************/

/* General Global Variables */
char listName[100];
FILE *fList;
//
char paramsFilename[100], iSiteFilename[100], bSiteFilename[100], basicSiteFilename[100];
FILE *paramsFile, *iSiteList, *bSiteList, *basicSiteList;

long NFil,N[NFILMAX], ntNextStationarityCheck, iBin;
long iSite[NFILMAX][NMAX], iSiteTotal[NFILMAX], iSiteCurrent, iy,ty, stericOcclusion[NFILMAX][NMAX];
long Ncurrent;
double c0, c1, irLigand;

double ree[NFILMAX], rM[NFILMAX], rM2[NFILMAX], rMiSite[NFILMAX][NMAX], rM2iSite[NFILMAX][NMAX], rH[NFILMAX], ksStatistic;

long iseed;

double phi[NFILMAX][NMAX], theta[NFILMAX][NMAX], psi[NFILMAX][NMAX];
double phiPropose[NMAX], thetaPropose[NMAX], psiPropose[NMAX];
double r[NFILMAX][NMAX][3],t[NFILMAX][NMAX][3], e1[NFILMAX][NMAX][3], e2[NFILMAX][NMAX][3],
       rBase[NFILMAX][3], tBase[NFILMAX][3], e1Base[NFILMAX][3], e2Base[NFILMAX][3],
       rPropose[NFILMAX][NMAX][3],tPropose[NFILMAX][NMAX][3], e1Propose[NFILMAX][NMAX][3], e2Propose[NFILMAX][NMAX][3];
double norm;
double iLigandCenter[NFILMAX][NMAX][3];

double RGlobal[3][3], RLocal[3][3];
double e1_dot_t, e2_dot_t, e2_dot_e1;

long st;
long nf, nf2, nfPropose;
long proposals[2], accepts[2], nt, iChi, i, iPropose, ix, iParam, ntNextStationarityCheck,i2, iStart;
double E, ENew, rate[2], dChi[2], dChiHere, Force;
long constraintProposalsTotal;

int iSiteInputMethod;
long commandiSites;
char *iSiteLocations;
char input[4*NMAX];
long j,m;

/* Convergence Global Variables */
int convergedTF, constraintSatisfiedTF, verboseTF;
long convergenceVariableCounts[NBINS], convergenceVariableCountsPrevious[NBINS];
long polymerLocationCounts[NFILMAX][NMAX][NBINSPOLYMER];

/* STIFFEN Global Variables */
double StiffenRange, StiffSites[NFILMAX][NMAX];
int stiffCase, totalStiff[NFILMAX];

char occupiedSites[4*NMAX],occupiedSitesNoSpace[NMAX];
double iSiteOccupied[NFILMAX][NMAX];

/* MULTIPLE Global Variables*/
int bSiteInputMethod;
double brLigand;
double bLigandCenter[NFILMAX][NMAX][3];
long bSite[NFILMAX][NMAX], bSiteTotal[NFILMAX], bSiteCurrent, ib, ib2;
long bSiteCounter;

double bLigandCenterPropose[NFILMAX][NMAX][3];

double deliveryDistance;
long stericOcclusionBase[NFILMAX];
long membraneOcclusion[NFILMAX][NMAX], membraneAndSegmentOcclusion[NFILMAX][NMAX];
double localConcCutoff;
int deliveryMethod;
long boundToBaseDeliver[NFILMAX][NMAX];

/* ELECTRO Global Variables */
double Eelectro, EelectroNew;
double Erepulsion, Zrepulsion;
double parabolaDepth, parabolaWidth, wallParabolaK;
double PhosphorylatedSites[NFILMAX][NMAX];
int PhosElectroRange;
long basicSite[NFILMAX][NMAX], BasicSitesYN[NFILMAX][NMAX], basicSiteTotal[NFILMAX], basicSiteCurrent, iBasic;

/* BASEBOUND Global Variables */
double baserLigand;
double baseLigandCenter[NFILMAX][3];
double baseCenter[3];

/* MULTIPLE FILAMENT Variables*/
double baseSepDistance;

/*******************************************************************************/
//  INCLUDES
/*******************************************************************************/

#include "outputControl.c"
#include "getParameters.c"
#include "getSites.c"
#include "initializeStiffSites.c"
#include "initializePhosphorylatedSites.c"
#include "getBasicSites.c"
#include "metropolisJoint.c"


/*******************************************************************************/
//  MAIN
/*******************************************************************************/

// arguments: listName,
int main( int argc, char *argv[] )
{
        if (!TXTPARAM)
        {
            printf("This program is starting.");
            
            if(argv[1]) // listName
                strcpy(listName, argv[1]);
            if (TALKATIVE) printf("This is argument 1: %s\n", listName);
            
            if(argv[2]) // N - should be 113 for this code - human CD3 iSites specified
                N = atoi(argv[2]);
            if (TALKATIVE) printf("This is argument 2: %ld\n", N);
            
            //- binding site. iSite=0 is the first joint away from the origin. iSite=N-1 is the furthest joint
            //if(argv[3]) // iSite
            //iSite = atoi(argv[3]);
            
            //if(iSite == -1)
            //iSite = floor(N/2);
            
            if(argv[3]) // irLigand - RATIO OF ligand radius to kuhn length
                irLigand = atof(argv[3]);
            if (TALKATIVE) printf("This is argument 3: %f\n", irLigand);
            
            if(argv[4]) // brLigand - RATIO OF ligand radius to kuhn length
                brLigand = atof(argv[4]);
            if (TALKATIVE) printf("This is argument 4: %f\n", brLigand);
            
            Force = 0;
            if(argv[5]) // Force - Units of kBT/[kuhn length]
                Force = atof(argv[5]);
            if (TALKATIVE) printf("This is argument 5: %f\n", Force);
            
            // IF verboseTF = 0, one line summarizing the run is written to the file listName.
            // IF verboseTF = 1, one line is written each iteration to the file listName. (Use for making histograms).
            verboseTF = 0;
            if(argv[6]) // Verbose Output
                verboseTF = atoi(argv[6]);
            if (TALKATIVE) printf("This is argument 6: %d\n", verboseTF);
            
            if(argv[7]) //iSiteInputMethod - switch for how iSites are input
                iSiteInputMethod = atoi(argv[7]);
            if (TALKATIVE) printf("This is argument 7: %d\n", iSiteInputMethod);
            
            if(argv[8])
            {
                if(atoi(argv[8])!=-1) //iSite Location from command line
                {
                    iSite[0]= atoi(argv[8]);
                    iSiteTotal=1;
                    if (TALKATIVE) printf("This is argument 8: %ld\n", iSite[0]);
                    iSiteInputMethod=1;
                }
            }
            
            if(argv[9])
            {
                if(atoi(argv[9])!=-1)
                {
                    bSite[0]=atoi(argv[9]);
                    bSiteTotal=1;
                    if (TALKATIVE) printf("This is argument 9: %ld\n", bSite[0]);
                    bSiteInputMethod = 1;
                }
                else
                {
                    bSiteInputMethod = 0;
                }
            }
            ///////////Stiffening Parameters/////////////////
            if(argv[10]) // Stiffness Range - 0 = stiffen only the iSite, -1 = no stiffening at all
                StiffenRange = atof(argv[10]);
            if (TALKATIVE) printf("This is argument 10: %f\n", StiffenRange);
            
            if(argv[11]) // Stiffness Case - 0 = CD3Zeta Mouse
                stiffCase = atoi(argv[11]);
            if (TALKATIVE) printf("This is argument 11: %d\n", stiffCase);
            
            if(argv[12]) // Occupied (phosphorylated) iSites
                strcpy(occupiedSites,argv[12]);
            if (TALKATIVE) printf("This is argument 12: %s\n", occupiedSites);

            if(argv[13]) // Occupied (phosphorylated) iSites
                strcpy(occupiedSitesNoSpace,argv[13]);
            if (TALKATIVE) printf("This is argument 13: %s \n", occupiedSitesNoSpace);
            
            if(argv[14]) //iSite file
                strcpy(iSiteFilename, argv[14]);
            if (TALKATIVE) printf("This is argument 14: %s \n", iSiteFilename);
            
            if(argv[15]) //bSite file
                strcpy(bSiteFilename, argv[15]);
            if (TALKATIVE) printf("This is argument 15: %s \n", bSiteFilename);
            
            
            if(argv[16]) //bSiteInputMethod - switch for how bSites are input
                bSiteInputMethod = atoi(argv[16]);
            if (TALKATIVE) printf("This is argument 16: %d \n", bSiteInputMethod);
            
            if(argv[17])
            {
                parabolaDepth = atof(argv[17]);
                if (TALKATIVE) printf("This is the parabola depth: %f\n", parabolaDepth);
            }
            
            if(argv[18])
            {
                parabolaWidth = atof(argv[18]);
                if (TALKATIVE) printf("This is the parabola width: %f\n", parabolaWidth);
            }
            
            if(argv[19])
            {
                wallParabolaK = atof(argv[19]);
                if (TALKATIVE) printf("This is the wall parabola K: %f\n", wallParabolaK);
            }
            
            if(argv[20])
            {
                Erepulsion = atof(argv[20]);
                if (TALKATIVE) printf("This is the Erepulsion: %f\n", Erepulsion);
            }
            
            if(argv[21])
            {
                Zrepulsion = atof(argv[21]);
                if (TALKATIVE) printf("This is the Zrepulsion: %f\n", Zrepulsion);
            }
            
            if(argv[22]) //PhosElectroRange
                PhosElectroRange = atof(argv[22]);
            if (TALKATIVE) printf("This is argument 22: %d \n", PhosElectroRange);

            
        }
        if (TXTPARAM) // majority of parameters are in text file
        {
            if(argv[1]) //get name of parameters file
                strcpy(paramsFilename,argv[1]);
            if (TALKATIVE) printf("This is the parameter filename: %s\n", paramsFilename);
            
            // use text file to set parameters
            getParameters();
            
            // use rest of command line input to overwrite parameters as necessary
            if(argv[2])
            {
                if(atoi(argv[2])!=-1)
                    strcpy(listName,argv[2]);
                if (TALKATIVE) printf("This is the output filename: %s\n", listName);
            }
            
            if(argv[3])
            {
                if(atoi(argv[3])!=-1)
                    verboseTF = atoi(argv[3]);
                if (TALKATIVE) printf("This will print verbose: %ld\n", verboseTF);
            }
            
            if(argv[4])
            {
                if(atoi(argv[4])!=-1)
                    NFil = atoi(argv[4]);
                if (TALKATIVE) printf("This is the number of filaments: %ld\n", NFil);
            }
            
            if(argv[5])
            {
                if(atoi(argv[5])!=-1)
                    baseSepDistance = atof(argv[5]);
                if (TALKATIVE) printf("This is the base separation distance: %lf\n", baseSepDistance);
            }
            
            if(argv[6])
            {
                if(atoi(argv[6])!=-1)
                    Force = atof(argv[6]);
                if (TALKATIVE) printf("This is the force: %lf\n", Force);
            }

        }

	iseed = RanInitReturnIseed(0);
	
	metropolisJoint();

	return 0;
	
} // finished main
