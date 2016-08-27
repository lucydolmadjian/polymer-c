/*** Allard Group jun.allard@uci.edu                    ***/

void outputGillespie();

/*******************************************************************************/
//  GLOBAL VARIABLES for output control
/*******************************************************************************/

double topPaths[NTOPPATHS][5],pathArrayShort[STATEMAX][5],maxPath[5],leastPath[5];
int factorial,frequency,topLocation,topPathsLocation[NTOPPATHS],topFrequency, maxFreq, maxLocation,leastLocation;
double MFTP,leastFreq;
int pass;

/********************************************************************************************************/
void outputGillespie()
{
    
    //find factorial
    factorial=1;
    for (i=1;i<=iSiteTotal;i++)
    {
        factorial *= i;
    }
    
//    printf("This is %d factorial: %d\n", iSiteTotal,factorial);
//    
//    for (i=0;i<factorial;i++)
//    {
//        for (j=0;j<5;j++)
//        {
//            pathArrayShort[i][j]=0;
//        }
//    }
    
    
    /************* MFTP ******************/
    MFTP = timeSum/iterations; //find mean first passage time - avg time it takes to move from 000000 to 111111
    
    
    
    /******************** PATHS ************************/
    
    //find the top twenty fastest paths
    //are fastest and most probable the same??
    //presumably no chance of having a tie?
    
    //this seems like a bad method
    
    //creates smaller matrix of all paths
    i=0;
    for (j=0;j<pow(10,iSiteTotal+1);j++)
    {
        if (pathArray[j][0] != 0)
        {
            pathArrayShort[i][0] = j; //path sequence
            pathArrayShort[i][1] = pathArray[j][0]; //total times path was taken
            pathArrayShort[i][2] = pathArray[j][1]; //sum of all times when path was taken
            pathArrayShort[i][3] = pathArray[j][0]/iterations; //probability of path being taken
            pathArrayShort[i][4] = pathArray[j][1]/pathArray[j][0]; //average time
            i++;
        }
        
    }
    
    
    //find most frequent path
    maxFreq = 0;
    for (j=0;j<factorial;j++)
    {
        if (pathArrayShort[j][1] > maxFreq)
        {
            maxFreq = pathArrayShort[j][1];
            maxLocation = j;
        }
    }
    
    for (k=0;k<5;k++)
    {
        maxPath[k] = pathArrayShort[maxLocation][k];
    }
    
    
    //find least frequent path
    leastFreq = INF;
    for (j=0;j<factorial;j++)
    {
        if (pathArrayShort[j][1] < leastFreq)
        {
            leastFreq = pathArrayShort[j][1];
            leastLocation = j;
        }
    }
    
    for (k=0;k<5;k++)
    {
        leastPath[k] = pathArrayShort[leastLocation][k];
    }
  
    
    //find top twenty most frequent paths
    topFrequency=0;
    for (k=0;k<NTOPPATHS;k++)
    {
        frequency=0;
        
        for (i=0;i<factorial;i++)
        {
            if (pathArrayShort[i][1]>frequency)
            {
                j=0;
                pass = 1;
                while (j<NTOPPATHS && pass)
                {
                    if (i != topPathsLocation[j])
                    {
                        j++;
                    }
                    else
                    {
                        pass = 0;
                    }
                }
                if (pass)
                {
                    frequency = pathArrayShort[i][1];
                    topPathsLocation[k] = i;
                }
            }
        }
        //topLocation = topPathsLocation[k];
        //topFrequency = pathArrayShort[topLocation][1];
    }

    
    //create list of top twenty paths with path, frequency, total time, probability, average time
    for (k=0;k<NTOPPATHS;k++)
    {
        topLocation = topPathsLocation[k];
        for (i=0;i<5;i++)
        {
            topPaths[k][i] = pathArrayShort[topLocation][i];
        }
    }
    
    /***************** OUTPUT ********************/
    
    //print summary data

    summaryOutputFile = fopen(summaryOutputName, "a");
    
    
    fprintf(summaryOutputFile, "%f\n", MFTP);
    
    for (i=0;i<5;i++)
    {
        fprintf(summaryOutputFile, "%f ", maxPath[i]);
    }
    
    fprintf(summaryOutputFile, "\n");
    
    for (i=0;i<5;i++)
    {
        fprintf(summaryOutputFile, "%f ", leastPath[i]);
    }

    fprintf(summaryOutputFile, "\n");
    
    for (i=0;i<NTOPPATHS;i++)
    {
        for (j=0;j<5;j++)
        {
            fprintf(summaryOutputFile, "%f ", topPaths[i][j]);
        }
        
        fprintf(summaryOutputFile, "\n");
    }
    
    fclose(summaryOutputFile);

    
    //print all data


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


