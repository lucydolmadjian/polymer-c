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
    
    //initializes basic residues to 0 
    for(i=0;i<N;i++)
    {
        BasicSitesYN[i] =0;
    }
    
    //Full residue list of basic yes or no
    for (i=0;i<N;i++)
    {
        for(iBasic=0;iBasic<basicSiteTotal;iBasic++)
        {
            if(i==basicSite[iBasic])
            {

                BasicSitesYN[i]=1; //set that joint to "basic"

            }
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

