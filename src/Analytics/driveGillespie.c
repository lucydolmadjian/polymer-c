/*** Allard Lab lclemens@uci.edu                   ***/

#define TWISTER genrand_real3()
#define NMAX       400
#define NTMAX      1e9
#define NTADAPT    20000
#define NTCHECK    200000
#define DCHIMIN    1e-4
#define NBINS      100
#define PI         3.14159265359
#define INF        1e14
#define ISITEMAX   9
#define STATEMAX   10000000000
#define ITMAX      1e9

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

char matrixName[1000];
FILE *ratesFile;
long iseed;

char outputName[1000];
FILE *outputFile;


double timeTotal,randTime[ISITEMAX],timeStep,timeSum;
double pathArray[STATEMAX][2];
int currentState,iy,it,iterations, stepCount,path;


double stateMatrix[STATEMAX][ISITEMAX];
int i,j;
int iSiteTotal,newState;

int sizeOfStateMatrix;
int verbose;

/*******************************************************************************/
//  INCLUDES
/*******************************************************************************/

#include "outputGillespie.c"
#include "runGillespie.c"


/*******************************************************************************/
//  MAIN
/*******************************************************************************/

// arguments: listName, dimension, nBar, N, dt
int main( int argc, char *argv[] )
{

    printf("This program is starting.\n");
    
    if(argv[1]) // matrixName
        strcpy(matrixName, argv[1]);
        printf("This is argument 1: %s\n", matrixName);
    
    if(argv[2]) //iSiteTotal
        iSiteTotal = atoi(argv[2]);
        printf("This is argument 2: %d\n", iSiteTotal);
    
    if(argv[3]) //total number of iterations
        iterations = atoi(argv[3]);
        printf("This is argument 3: %d\n", iterations);
    
    if(argv[4]) //output file name
        strcpy(outputName, argv[4]);
        printf("This is argument 4: %s\n", outputName);
    
    if(argv[5]) //verbose - 1 is verbose, 0 is not verbose
        verbose = atoi(argv[5])
        printf("This is argument 5: %d\n", verbose);

    
	iseed = RanInitReturnIseed(0);
	
	runGillespie();

	return 0;
	
} // finished main
