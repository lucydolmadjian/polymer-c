/*** Allard Group jun.allard@uci.edu                    ***/

void outputGillespie();

/*******************************************************************************/
//  GLOBAL VARIABLES for output control
/*******************************************************************************/

double finalBinaryState,averagePhosphorylationState;
double sum_phosphorylationState,sum_time;
int finalPhosphorylationState;

/********************************************************************************************************/
void outputGillespie()
{
    // initializations
    sum_phosphorylationState = 0;
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
        sum_phosphorylationState += (double) timeStorage[i] * ((double) stateStorage[i] / (double) iSiteTotal);
        sum_time                 += (double) timeStorage[i];
    }
    averagePhosphorylationState = (double) sum_phosphorylationState/ (double) sum_time;
    
    
    // look at final state
    finalPhosphorylationState = totalBound[finalState];
    finalBinaryState = binaryState[finalState];
    
    /***************** OUTPUT ********************/

    outputFile = fopen(outputName, "a");

    
    fprintf(outputFile, "%e ", reverseRate);
    fprintf(outputFile, "%d ", iSiteTotal);
    fprintf(outputFile, "%f ", timeEnd);
    fprintf(outputFile, "%f ", averagePhosphorylationState);
    fprintf(outputFile, "%d ", finalState);
    fprintf(outputFile, "%d ", finalPhosphorylationState);
    fprintf(outputFile, "%f ", finalBinaryState);
    fprintf(outputFile, "%f ", finalTotalTime);
    fprintf(outputFile, "\n");

    fclose(outputFile);

    
}


