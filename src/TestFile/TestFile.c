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

long N, ntNextStationarityCheck, iBin;

long iSite[NMAX], iSiteTot, iSiteCurrent, iy, stericOcclusion[NMAX];
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
double rLigandCenter[3];
double RGlobal[3][3], RLocal[3][3];
double e1_dot_t, e2_dot_t, e2_dot_e1;

long proposals[2], accepts[2], nt, iChi, i, iPropose, ix, iParam, ntNextStationarityCheck,i2, iStart;

double E, ENew, rate[2], dChi[2], dChiHere, ksStatistic, Force;

int convergedTF, constraintSatisfiedTF, verboseTF;

/*******************************************************************************/
//  INCLUDES
/*******************************************************************************/

//#include "outputControl.c"
//#include "metropolisJoint.c"

/*******************************************************************************/
//  MAIN
/*******************************************************************************/

// arguments: listName, dimension, nBar, N, dt
int main( int argc, char *argv[] )
{

	if(argv[1]) // listName
		strcpy(listName, argv[1]);

	if(argv[2]) // N
		N = atoi(argv[2]);
	
    //- binding site. iSite=0 is the first joint away from the origin. iSite=N-1 is the furthest joint
	//if(argv[3]) // iSite
		//iSite = atoi(argv[3]);
    
    //if(iSite == -1)
        //iSite = floor(N/2);
		
	if(argv[3]) // rLigand - RATIO OF ligand radius to kuhn length
		rLigand = atof(argv[3]);
    
    Force = 0;
    if(argv[4]) // Force - Units of kBT/[kuhn length]
        Force = atof(argv[4]);
    
    // IF verboseTF = 0, one line summarizing the run is written to the file listName.
    // IF verboseTF = 1, one line is written each iteration to the file listName. (Use for making histograms).
    verboseTF = 0;
    if(argv[5]) // Verbose Output
        verboseTF = atoi(argv[5]);
    
	iseed = RanInitReturnIseed(0);
	
	metropolisJoint();

	return 0;
	
} // finished main
