/*** Allard Group jun.allard@uci.edu                    ***/

void runGillespie();

void runGillespie()
{
    
    /****************************** Create State Matrix from File ***********************/

    sizeOfStateMatrix = (int) pow(2,iSiteTotal);
    
    
    char line[200];
    i=0;
    j=0;
    
    ratesFile = fopen(matrixName,"r");
    
    while (fgets(line, sizeof(line), ratesFile))
    {
        stateMatrix[i][j]=atof(line);
        j++;
        if ((j%iSiteTotal) == 0)
        {
            i++;
            j=0;
        }
            
    }
    
    fclose(ratesFile);
    
    // replace 0 rates with reverseRate variable
    
    for (i=0;i<sizeOfStateMatrix;i++)
    {
        for (j=0;j<iSiteTotal;j++)
        {
            if (stateMatrix[i][j]==0)
            {
                stateMatrix[i][j] = reverseRate;
            }
        }
    }
            
    
    //debugging
    if(1)
    {
        
        // need to improve this debugging part - won't work if reverse = forward in any case
//        // count instances of reverseRate
//        for (i=0;i<sizeOfStateMatrix;i++)
//        {
//            for (j=0;j<iSiteTotal;j++)
//            {
//                if (stateMatrix[i][j]==reverseRate)
//                {
//                    reverseRateTotalInstances++;
//                }
//            }
//        }
        
        printf("Total reverse: %f\n",reverseRateTotalInstances);
        
        // print state matrix
        for (i=0;i<sizeOfStateMatrix;i++)
        {
            for (j=0;j<iSiteTotal;j++)
            {
                printf("%lf ", stateMatrix[i][j]);
            
            }
        
            printf("\n");
        }
    }
    
    

    
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
                if (stateMatrix[currentState][iy]!=0)
                {
                    randTime[iy] = - log(TWISTER)/stateMatrix[currentState][iy]; //exponentially distributed random variable based on transition rate
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



