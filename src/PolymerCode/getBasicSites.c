/*** Allard Group jun.allard@uci.edu                    ***/

void getBasicSites();


/*******************************************************************************/
//  GLOBAL VARIABLES for output control
/*******************************************************************************/


/********************************************************************************************************/
void getBasicSites()
{
    /********* INITIALIZE BASIC AND TYROSINE SITES *******************/
    

    basicSiteList = fopen(basicSiteFilename, "r");
    char line[200];
    iBasic=0;
    
    while (fgets(line, sizeof(line), basicSiteList))
    {
        basicSite[iBasic]=atoi(line);
        iBasic++;
    }
    
    fclose(iSiteList);
    
    basicSiteTotal=iBasic;

    
    for (iBasic=0; iBasic<basicSiteTotal;iBasic++)
    {
        if (basicSite[iBasic] >= N)
        {
            printf("Warning! Site is located past end of polymer!");
            fflush(stdout);
        }
    }
    
    //for debugging - prints a list of the iSites
    
    if (TALKATIVE)
    {
        for (iBasic=0;iBasic<basicSiteTotal;iBasic++)
        {
            printf("basicSite: %ld\n", basicSite[iBasic]);
            fflush(stdout);
        }
        
        printf("basicSiteTotal: %ld\n", basicSiteTotal);
        fflush(stdout);
    }


}

/********************************************************************************************************/

