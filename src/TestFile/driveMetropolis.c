/*** Allard Lab jun.allard@uci.edu                    ***/

#define TWISTER genrand_real3()
#define NMAX     400
#define NTMAX    1e9
#define NTADAPT  20000
#define NTCHECK  200000
#define DCHIMIN  1e-4
#define NBINS    100
#define PI       3.14159265359
#define INF      1e14
#define DCHIINIT 0.1
#define KSCRITICAL 0.005
#define MEMBRANE 0
#define MULTIPLE 1

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

char listName[100];
FILE *fList;

char listName2[100] = "TestOutput";
FILE *fList2;

long N, ntNextStationarityCheck, iBin;

long iSite[NMAX], iSiteTot, iSiteCurrent, iy,ty, stericOcclusion[NMAX];
double c0, c1, rLigand;
double ree, rM, rH, ksStatistic;
long rMCounts[NBINS], rMCountsPrevious[NBINS];
long iseed;

double phi[NMAX], theta[NMAX], psi[NMAX];
double phiPropose[NMAX], thetaPropose[NMAX], psiPropose[NMAX];
double r[NMAX][3],t[NMAX][3], e1[NMAX][3], e2[NMAX][3],
       rBase[3], tBase[3], e1Base[3], e2Base[3],
       rPropose[NMAX][3],tPropose[NMAX][3], e1Propose[NMAX][3], e2Propose[NMAX][3];
double norm;
double rLigandCenter[NMAX][3];
double RGlobal[3][3], RLocal[3][3];
double e1_dot_t, e2_dot_t, e2_dot_e1;

double StiffenRange,phosiSites[NMAX], Stiff[NMAX];
char phosphorylatediSites[4*NMAX],phosphorylatediSitesNoSpace[NMAX];

double rLigandCenterBound[NMAX][3];
long iSiteBound[NMAX], boundTotal, currentBoundSite, ib, ib2;

double rLigandCenterBase[3], deliveryDistance;
long stericOcclusionBase;
int deliveryMethod;

long boundToBaseDeliver[NMAX];

long constraintProposals, constraintProposalsTotal;

//long commandiSites;
//char *iSiteLocations;
//char input[4*NMAX];

long st;


long proposals[2], accepts[2], nt, iChi, i, iPropose, ix, iParam, ntNextStationarityCheck,i2, iStart;

double E, ENew, rate[2], dChi[2], dChiHere, ksStatistic, Force;

int convergedTF, constraintSatisfiedTF, verboseTF, testRun;

/*******************************************************************************/
//  INCLUDES
/*******************************************************************************/

#include "outputControl.c"
#include "metropolisJoint.c"

/*******************************************************************************/
//  MAIN
/*******************************************************************************/

// arguments: listName, dimension, nBar, N, dt
int main( int argc, char *argv[] )
{
    
    printf("This program is starting.");
    
	if(argv[1]) // listName
		strcpy(listName, argv[1]);
        printf("This is argument 1: %s\n", listName);

	if(argv[2]) // N - should be 143 for this code - human CD3 iSites specified
		N = atoi(argv[2]);
    printf("This is argument 2: %ld\n", N);
	
    //- binding site. iSite=0 is the first joint away from the origin. iSite=N-1 is the furthest joint
	//if(argv[3]) // iSite
		//iSite = atoi(argv[3]);
    
    //if(iSite == -1)
        //iSite = floor(N/2);
		
	if(argv[3]) // rLigand - RATIO OF ligand radius to kuhn length
		rLigand = atof(argv[3]);
    printf("This is argument 3: %f\n", rLigand);
    
    Force = 0;
    if(argv[4]) // Force - Units of kBT/[kuhn length]
        Force = atof(argv[4]);
    printf("This is argument 4: %f\n", Force);
    
    if(argv[5]) // Occupied (phosphorylated) iSites
        strcpy(phosphorylatediSites,argv[5]);
        printf("This is argument 5: %s\n", phosphorylatediSites);
    
    if(argv[6]) // Stiffness Range - 0 = stiffen only the iSite, -1 = no stiffening at all
        StiffenRange = atof(argv[6]);
    printf("This is argument 6: %f\n", StiffenRange);
    
    
    // IF verboseTF = 0, one line summarizing the run is written to the file listName.
    // IF verboseTF = 1, one line is written each iteration to the file listName. (Use for making histograms).
    verboseTF = 0;
    if(argv[7]) // Verbose Output
        verboseTF = atoi(argv[7]);
    printf("This is argument 7: %d\n", verboseTF);
    
    if(argv[8]) //Test Run - yes=1, no=0
        testRun = atoi(argv[8]);
    printf("This is argument 8: %d\n", testRun);
    
    if(argv[9]) // Occupied (phosphorylated) iSites
        strcpy(phosphorylatediSitesNoSpace,argv[9]);
    printf("This is argument 9: %s\n", phosphorylatediSitesNoSpace);
    
    if(argv[10]) //Delivery distance - how close to base it needs to be
        deliveryDistance = atof(argv[10]);
    printf("This is argument 10: %f\n", deliveryDistance);
    
    if(argv[11]) //Delivery method - 0 = within Base ligand site, 1 = within deliveryDistance
        deliveryMethod = atoi(argv[11]);
    printf("This is argument 11: %d\n", deliveryMethod);
    
//    if(argv[12]) //hardcoded vs command line iSites
//        commandiSites = atoi(argv[12]);
//    printf("This is argument 12: %d/n", commandiSites);
//
//    if (commandiSites)
//    {
//        if(argv[13])
//            iSiteTot=atoi(argv[13]);
//        printf("This is argument 13: %d/n", iSiteTot);
//        
//        if(argv[14])
//            strcpy(input,argv[14]);
//            strcpy(iSiteLocations,argv[14]);
//            printf("This is argument 14: %s and %s", iSiteLocations, input);
//    }
    
    
	iseed = RanInitReturnIseed(0);
	
	metropolisJoint();

	return 0;
	
} // finished main
