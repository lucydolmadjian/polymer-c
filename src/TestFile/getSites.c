/*** Allard Group jun.allard@uci.edu                    ***/

void getSites();


/*******************************************************************************/
//  GLOBAL VARIABLES for output control
/*******************************************************************************/


/********************************************************************************************************/
void getSites()
{
    /********* INITIALIZE ISITES *******************/
    
   
    switch (testRun) //switch in Batch Script to specify which set of iSites you want to use
    {
        case 0:  // iSites initialized for human CD3Zeta-Chain
            
            iSiteTotal = 7;
            
            
            for (iy=0;iy<iSiteTotal;iy++) //initializes iSite array
            {
                iSite[iy]=0;
            }
            
            //specify iSite locations by rod number (located at n-1) (i.e. if you have a polymer of N=50, and want an iSite at the 36th rod, input iSite[]=35)
            
            iSite[0]=42;
            iSite[1]=50;
            iSite[2]=61;
            iSite[3]=89;
            iSite[4]=101;
            iSite[5]=120;
            iSite[6]=131;
            break;
            
        case 1: // iSites for formin //for testing - N=10
            
            iSiteTotal = 3;
            
            for (iy=0;iy<iSiteTotal;iy++)
            {
                iSite[iy]=0;
            }
            
            iSite[0]=0;
            iSite[1]=3;
            iSite[2]=7;
            break;

            
        case 2: //iSites for testing
            
            iSiteTotal = 4;
            
            for(iy=0;iy<iSiteTotal;iy++)
            {
                iSite[iy]=0;
            }
            
            iSite[0]=10;
            iSite[1]=11;
            iSite[2]=25;
            iSite[3]=49;

            break;

    }
    
    
    //for debugging - prints a list of the iSites
    
//    for (iy=0;iy<iSiteTotal;iy++)
//    {
//        printf("iSite: %ld\n", iSite[iy]);
//        fflush(stdout);
//    }
    
    /********* INITIALIZE BOUND ISITES *******************/
    
    if (MULTIPLE) //if looking at multiple binding (i.e. MULTIPLE set to 1 in driveM)
    {
        
        bSiteTotal = 1; //total number of iSites bound
        
        bSite[0]=25; //specify bSites same way as iSites

    }


}

/********************************************************************************************************/

