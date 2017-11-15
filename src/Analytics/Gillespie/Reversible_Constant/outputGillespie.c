/*** Allard Group jun.allard@uci.edu                    ***/

void outputGillespie();

/*******************************************************************************/
//  GLOBAL VARIABLES for output control
/*******************************************************************************/

double finalBinaryState,averagePhosphorylationState;
int sum_phosphorylationState,finalPhosphorylationState;

/********************************************************************************************************/
void outputGillespie()
{
    finalPhosphorylationState = totalBound[finalState];
    finalBinaryState = binaryState[finalState];
    
    // find average phosphorylation state over last 100 states
    for(i=0;i<numberStatesStored;i++)
    {
        sum_phosphorylationState += stateStorage[i];
    }
    averagePhosphorylationState = (double) sum_phosphorylationState/ (double) numberStatesStored;
    
    
    
    /***************** OUTPUT ********************/

    outputFile = fopen(outputName, "a");

    
    fprintf(outputFile, "%f ", reverseRate);
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


