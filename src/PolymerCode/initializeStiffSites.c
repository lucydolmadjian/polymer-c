/*** Allard Group jun.allard@uci.edu                    ***/

void initializeStiffSites();

/*******************************************************************************/
//  GLOBAL VARIABLES for initializing stiff sites
/*******************************************************************************/
int stiffEnd, stiffStart;
double stiffiSites[NCHAINMAX][NMAX];
//
/********************************************************************************************************/
void initializeStiffSites()
{
    //initializes stiffiSites to 0 (none phosphorylated)
    for(ty=0;ty<iSiteTotal;ty++)
    {
        stiffiSites[ty]=0;
    }

    if (TALKATIVE) printf("These are the occupied sites: %s\n", occupiedSites);
    
    //read string and assign to double vector
    // 1 is occupied iSite (phosphorylated), 0 is unoccupied
    switch (stiffCase)
    {
        case 0:
            
            sscanf(occupiedSites,"%lf_%lf_%lf_%lf_%lf_%lf", &stiffiSites[0],&stiffiSites[1],&stiffiSites[2],&stiffiSites[3], &stiffiSites[4],&stiffiSites[5]);
            break;
            
            // could include a case to read occupiedSites from a file of either locations or of 0,1s
            // for now, this is only useful for CD3Zeta parameters (anything with six iSites)
            // could get rid of switch altogether for now
    }
    
    // for debugging, print which iSites are declared stiff
    if (TALKATIVE)
        for (iy=0;iy<iSiteTotal;iy++)
        {
            printf("stiffiSites[ %ld ] =  %f\n",iy, stiffiSites[iy]);
            fflush(stdout);
        }
    
    //initializes stiffened rods to 0 (none stiff)
    for(i=0;i<N;i++)
    {
        StiffSites[i] =0;
    }

    /********************************************************/
    /******************* STIFFEN SEGMENTS *******************/
    /********************************************************/
    for(ty=0;ty<iSiteTotal;ty++)
    {
        if(stiffiSites[ty]==1) //might want to check the truth value on this - equals for double?
        {
            // set beginning of stiffening range
            if(iSite[ty]-StiffenRange >= 0)
            {
                stiffStart=iSite[ty]-StiffenRange;
            }
            else // if stiffenrange goes below 0, start stiffening at 0
            {
                stiffStart=0;
            }
            
            //set end of stiffening range
            // if stiffenrange goes above N, end at N
            if(iSite[ty]+StiffenRange+1 >= N)
            {
                stiffEnd=N;
            }
            else
            {
                stiffEnd=iSite[ty]+StiffenRange+1;
            }
            
            // declare which segments are stiff, exclusive of right endpoint (because included +1 above)
            for(i=stiffStart;i<stiffEnd;i++)
            {
                StiffSites[i]=1; //set that joint to "stiff"
            }
        }
    }

    // for debugging, count and print the total number of stiff segments
    if (TALKATIVE)
    {
        totalStiff = 0;
        for (i=0;i<N;i++)
        {
            if (StiffSites[i]==1)
            {
                //printf("Stiffen[ %ld ] =  %f\n",i, StiffSites[i]);
                //fflush(stdout);
                totalStiff++;
            }
        }
        printf("Total Stiff: %d\n", totalStiff);
        fflush(stdout);
        
        if (totalStiff >= N)
        {
            // Include error for completely stiff
            // May cause convergence problems?
            printf("Error! Completely stiff!\n");
            fflush(stdout);
            
            exit(0);
        }
    }
}

/********************************************************************************************************/

