/*** Allard Group jun.allard@uci.edu                    ***/

void outputGillespie();

/*******************************************************************************/
//  GLOBAL VARIABLES for output control
/*******************************************************************************/

double reeBar_sum, ree2Bar_sum, rMBar_sum;
long POcclude_sum[NMAX], Prvec0_sum[NMAX], POccludeBase_sum, PDeliver_sum[NMAX];
double topPaths[20][5],pathArrayShort[STATEMAX][3];
int factorial;

/********************************************************************************************************/
void initializeSummary()
{
    // summary variables
    reeBar_sum   = 0;
    ree2Bar_sum  = 0;


}

/********************************************************************************************************/
void outputGillespie()
{
    
    //find the top twenty fastest paths
    //are fastest and most probable the same??
    //presumably no chance of having a tie?
    
    //this seems like a bad method
    
    /********** Factorial ***********/
    
    //creates smaller matrix of all paths
    i=0;
    for (j=0;j<STATEMAX;j++)
    {
        if (pathArray[j][0] != 0)
        {
            pathArrayShort[i][0] = j; //path sequence
            pathArrayShort[i][1] = pathArray[j][0]; //total times path was taken
            pathArrayShort[i][2] = pathArray[j][1]; //sum of all times when path was taken
        }
    }
    
    for (i=1;i<=iSiteTotal;i++)
    {
        factorial *= i;
    }
    
    for (i=0;i<factorial;i++)
    {
        if (pathArrayShort[i][1]>freqPath)
    }
    
    topPaths[i][j]
    
    
    if (!verbose)
    {
        outputFile = fopen(outputName, "a");
        
        
        fprintf(outputFile, "%ld %f %f %f %ld %f %f %f %e",
                N,           // 1
                irLigand,    // 2
                brLigand,    // 3
                Force,       // 4
                nt,          // 5
                ksStatistic, // 6
                reeBar,      // 7
                ree2Bar,     // 8
                rMBar);      // 9
        
        for (iy=0;iy<iSiteTotal;iy++)
        {
            fprintf(fList, " %ld %e %e %e",
                iSite[iy], //10 + 4*iBind
                POcclude[iy], //11 + 4*iBind
                1-POcclude[iy], //12 + 4*iBind
                Prvec0[iy]); //13 + 4*iBind
        
        }
        
        fprintf(fList, " %d %e %e", -1, POccludeBase, 1-POccludeBase);
        
        for (ib=0;ib<bSiteTotal;ib++)
        {
            fprintf(fList, " %ld",
                bSite[ib]);
                    //PDeliver[ib]);
        }
        
        fprintf(fList, "\n");
        fclose(fList);
    }
    
    
    if (verbose)
    {
    
    }
    
}

	
}

