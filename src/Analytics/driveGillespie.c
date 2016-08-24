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


//

double transitionMatrix[64][64];

int i,j;


/*******************************************************************************/
//  INCLUDES
/*******************************************************************************/

//#include "outputGillespie.c"
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

    
	iseed = RanInitReturnIseed(0);
	
	runGillespie();

	return 0;
	
} // finished main
