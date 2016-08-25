/*** Allard Group jun.allard@uci.edu                    ***/

void outputGillespie();

/*******************************************************************************/
//  GLOBAL VARIABLES for output control
/*******************************************************************************/

double topPaths[20][5],pathArrayShort[STATEMAX][5];
int factorial,frequency,topLocation,topPathsLocation[20],topFrequency;
double MFTP;

/********************************************************************************************************/
void outputGillespie()
{
    
    /************* MFTP ******************/
    MFTP = timeSum/iterations; //find mean first passage time - avg time it takes to move from 000000 to 111111
    
    
    
    /******************** PATHS ************************/
    
    //find the top twenty fastest paths
    //are fastest and most probable the same??
    //presumably no chance of having a tie?
    
    //this seems like a bad method
    
    //creates smaller matrix of all paths
    i=0;
    for (j=0;j<STATEMAX;j++)
    {
        if (pathArray[j][0] != 0)
        {
            pathArrayShort[i][0] = j; //path sequence
            pathArrayShort[i][1] = pathArray[j][0]; //total times path was taken
            pathArrayShort[i][2] = pathArray[j][1]; //sum of all times when path was taken
            pathArrayShort[i][3] = pathArray[j][0]/iterations; //probability of path being taken
            pathArrayShort[i][4] = pathArray[j][1]/pathArray[j][0]; //average time
        }
    }
    
    //find factorial
    for (i=1;i<=iSiteTotal;i++)
    {
        factorial *= i;
    }
    
    //find top twenty most frequent paths
    topFrequency=0;
    for (k=0;k<20;k++)
    {
        frequency=0;
        
        for (i=0;i<factorial;i++)
        {
            if (pathArrayShort[i][1]>frequency)
            {
                if (pathArrayShort[i][1]<topFrequency)
                {
                    frequency = pathArrayShort[i][1];
                    topPathsLocation[k] = i;
                }
            }
        }
        topLocation = topPathsLocation[k];
        topFrequency = pathArrayShort[topLocation][1];
    }
    
    //create list of top twenty paths with path, frequency, total time, probability, average time
    for (k=0;k<20;k++)
    {
        topLocation = topPathsLocation[k];
        for (i=0;i<5;i++)
        {
            topPaths[k][i] = pathArrayShort[topLocation][i];
        }
    }
    
    /***************** OUTPUT ********************/
    
    //print top twenty data
    if (!verbose)
    {
        outputFile = fopen(outputName, "a");
        
        
        fprintf(outputFile, "%f\n", MFTP);

        
        for (i=0;i<20;i++)
        {
            for (j=0;j<5;j++)
            {
                fprintf(outputFile, "%f ", topPaths[i][j]);
            }
            
            fprintf(outputFile, "\n");
        }
        
        fclose(outputFile);
    }
    
    //print all data
    if (verbose)
    {
        outputFile = fopen(outputName, "a");
        
        
        fprintf(outputFile, "%f\n", MFTP);
        
        //print all 720 path data
        for (i=0;i<factorial;i++)
        {
            for (j=0;j<5;j++)
            {
                fprintf(outputFile, "%f ", pathArrayShort[i][j]);
            }
            
            fprintf(outputFile, "\n");
        }
        
        fclose(outputFile);
        
        
    }
    
}


