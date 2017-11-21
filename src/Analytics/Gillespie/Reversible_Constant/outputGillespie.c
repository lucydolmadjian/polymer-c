/*** Allard Group jun.allard@uci.edu                    ***/

void outputGillespie();

/*******************************************************************************/
//  GLOBAL VARIABLES for output control
/*******************************************************************************/

double finalBinaryState,averagePhosphorylationFraction,averageBound;
double sum_phosphorylationFraction,sum_time,sum_bound;
int finalPhosphorylationState;

/********************************************************************************************************/
void outputGillespie()
{
    // initializations
    sum_bound = 0;
    sum_phosphorylationFraction = 0;
    sum_time = 0;
    

    // find average phosphorylation state over last 100 states
//    for(i=0;i<numberStatesStored;i++)
//    {
//        sum_phosphorylationState += stateStorage[i];
//    }
//    averagePhosphorylationState = (double) sum_phosphorylationState/ (double) numberStatesStored;
//
    
    
    // find average phosphorylation state over last 100 states, taking time spent in state into account
    for(i=0;i<numberStatesStored;i++)
    {
        sum_bound                   += (double) timeStorage[i] * (double) stateStorage[i];
        sum_phosphorylationFraction += (double) timeStorage[i] * ((double) stateStorage[i] / (double) iSiteTotal);
        
        sum_time                    += (double) timeStorage[i];
    }
    averageBound                    = (double) sum_bound / (double) sum_time;
    averagePhosphorylationFraction  = (double) sum_phosphorylationFraction/ (double) sum_time;
    
    
    // look at final state
    finalPhosphorylationState = totalBound[finalState];
    finalBinaryState          = binaryState[finalState];
    
    /***************** OUTPUT ********************/

    outputFile = fopen(outputName, "a");

    
    fprintf(outputFile, "%e ", reverseRate);                        // 1
    fprintf(outputFile, "%d ", iSiteTotal);                         // 2
    fprintf(outputFile, "%f ", timeEnd);                            // 3
    fprintf(outputFile, "%f ", averagePhosphorylationFraction);     // 4
    fprintf(outputFile, "%f ", averageBound);                       // 5
    fprintf(outputFile, "%d ", finalState);                         // 6
    fprintf(outputFile, "%d ", finalPhosphorylationState);          // 7
    fprintf(outputFile, "%f ", finalBinaryState);                   // 8
    fprintf(outputFile, "%f ", finalTotalTime);                     // 9
    fprintf(outputFile, "\n");

    fclose(outputFile);

    
}


