/*** Allard Group jun.allard@uci.edu                    ***/

void runGillespie();

void runGillespie()
{
    
    /****************************** Create State Matrix from File ***********************/

    sizeOfRateMatrix = (int) pow(2,iSiteTotal);
    
    // import first rate matrix
    char line1[200];
    i=0;
    j=0;
    
    ratesFile = fopen(matrixName1,"r");
    
    while (fgets(line1, sizeof(line1), ratesFile))
    {
        rateMatrix1[i][j]=atof(line1);
        j++;
        if ((j%iSiteTotal) == 0)
        {
            i++;
            j=0;
        }
            
    }
    
    fclose(ratesFile);
    
    // import second rate matrix
    char line2[200];
    i=0;
    j=0;
    
    ratesFile = fopen(matrixName2,"r");
    
    while (fgets(line2, sizeof(line2), ratesFile))
    {
        // import reverse values, multiply by variable constant
        rateMatrix2[i][j]=reverseRate*atof(line2);
        j++;
        if ((j%iSiteTotal) == 0)
        {
            i++;
            j=0;
        }
        
    }
    
    fclose(ratesFile);
    
    // concatenate matrices into full rate matrix
    
    for (i=0;i<sizeOfRateMatrix;i++)
    {
        for (j=0;j<iSiteTotal;j++)
        {
            //efficient, but not as safe - include if/else/reject statements instead?
            rateMatrix[i][j] = rateMatrix1[i][j] + rateMatrix2[i][j];
        }
    }
            
    
    //debugging
    if(1)
    {
        
        
        // print rate matrix
        for (i=0;i<sizeOfRateMatrix;i++)
        {
            for (j=0;j<iSiteTotal;j++)
            {
                printf("%lf ", rateMatrix[i][j]);
            
            }
        
            printf("\n");
        }
    }
    
    /******************************* Binary Conversion ******************************************/
    
    binaryConversion();
    
    /******************************* Gillespie ******************************************/
    
    
    timeSum=0;
    it=0;

    timeTotal=0;
    
    currentState=0; //start at fully dephosphorylated
    stepCount = 0;
    
    // while less than number of desired steps or less than max steps
    while (it < iterations && it < ITMAX)
    {
        
            //initialize random time array and time step
            for (iy=0;iy<iSiteTotal;iy++)
            {
                randTime[iy]=0;
            }
            
            timeStep = INF;
            
            //Gillespie step
            
            for (iy=0;iy<iSiteTotal;iy++)
            {
                if (rateMatrix[currentState][iy]!=0)
                {
                    randTime[iy] = - log(TWISTER)/rateMatrix[currentState][iy]; //exponentially distributed random variable based on transition rate
                }
                else //this should be redundant since replaced zeros with reverse rates
                {
                    randTime[iy] = 0; //use 0 instead of infinity - then just remove these cases later
                }
            }
            
            //pick smallest of random times
            for (iy=0;iy<iSiteTotal;iy++)
            {
                if (randTime[iy]!= 0)  // 0 time is not an option, should be redundant
                {
                    if (randTime[iy]<timeStep)
                    {
                        timeStep = randTime[iy];
                        newState = iy;
                    }
                }
            }
            
            if (0)
            printf("This is the path chosen: %d\n", newState);
            
            //update time
            timeTotal += timeStep;
            stepCount++;
            
            // use this update for "forwards" transitionMatrix (i.e. forwards binary, backwards phosphorylation)
            currentState += pow(2,newState);
            
            // use this update for "backwards" transitionMatrix
            //currentState += pow(2,iSiteTotal-newState-1);
            
            if (0)
            printf("Current State is: %d \n", currentState);
            
            
        //}
        
        
        //for MFPT
        timeSum += timeTotal;
        
        
        if (TIME)
        {
            timeArray[it] = timeTotal;
        }
        
        it++;
    
    }
    
    outputGillespie();
    
    
}



