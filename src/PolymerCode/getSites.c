/*** Allard Group jun.allard@uci.edu                    ***/

void getSites();


/*******************************************************************************/
//  GLOBAL VARIABLES for output control
/*******************************************************************************/


/********************************************************************************************************/
void getSites()
{
    /********* INITIALIZE ISITES *******************/
    
   
    switch (iSiteInputMethod) //switch in Batch Script to specify which set of iSites you want to use
    {
            
        case 0:  // iSites initialized for human CD3Zeta-Chain

            iSiteTotal = 6;

            for (iy=0;iy<iSiteTotal;iy++) //initializes iSite array
            {
                iSite[iy]=0;
            }

            //specify iSite locations by rod number (located at n-1) (i.e. if you have a polymer of N=50, and want an iSite at the 36th rod, input iSite[]=35)

            iSite[0]=20;
            iSite[1]=31;
            iSite[2]=59;
            iSite[3]=71;
            iSite[4]=90;
            iSite[5]=101;
            break;
            
        case 1: //do nothing, use command line input
            
            break;
            
        case 2: //input iSites from file
            
            iSiteList = fopen(iSiteFilename, "r");
            char line[200];
            iy=0;
            
            while (fgets(line, sizeof(line), iSiteList))
            {
                iSite[iy]=atoi(line);
                iy++;
            }
            
            fclose(iSiteList);
            
            iSiteTotal=iy;

            break;
        
        case 3: // use last site as only iSite
        
            iSiteTotal = 1;
            for(iy=0;iy<iSiteTotal;iy++)
            {
                iSite[iy]=0;
            }
            iSite[0] = N-1;
        
            break;

    }
    
    //Warning for possible user error
    for (iy=0; iy<iSiteTotal;iy++)
    {
        if (iSite[iy] >= N)
        {
            printf("Warning! Site is located past end of polymer!");
            fflush(stdout);
        }
    }
    
    //for debugging - prints a list of the iSites
    if (TALKATIVE)
    {
        for (iy=0;iy<iSiteTotal;iy++)
        {
            printf("iSite: %ld\n", iSite[iy]);
            fflush(stdout);
        }
        
        printf("iSiteTotal: %ld\n", iSiteTotal);
        fflush(stdout);
    }
    
    /*****************************************************/
    /********* INITIALIZE BOUND ISITES *******************/
    /*****************************************************/
    
    if (MULTIPLE) //if looking at multiple binding (i.e. MULTIPLE set to 1 in driveM)
    {
        switch (bSiteInputMethod)
        {
            case 0:
        
                bSiteTotal = 4; //total number of iSites bound
                
                bSite[0]=10;
                bSite[1]=10;
                bSite[2]=40;
                bSite[3]=40;
                break;

            case 1: //do nothing, use command line input
                
                break;
                
            case 2: //bSites for multiple binding of ZAP-70 to CD3 Zeta mouse
                
                for (iy=0; iy<iSiteTotal; iy++)
                {
                    iSiteOccupied[iy]=0;
                }
                
                sscanf(occupiedSites,"%lf_%lf_%lf_%lf_%lf_%lf", &iSiteOccupied[0],&iSiteOccupied[1],&iSiteOccupied[2],&iSiteOccupied[3], &iSiteOccupied[4],&iSiteOccupied[5]);
                
                bSiteCounter=0;
                for (iy=0;iy<iSiteTotal;iy++)
                {
                    if (iSiteOccupied[iy]==1)
                    {
                        bSite[bSiteCounter]=iSite[iy];
                        bSiteCounter++;
                    }
                }
                
                bSiteTotal=bSiteCounter;
                
                break;
                
            case 3: // read bound sites from file
                
                bSiteList = fopen(bSiteFilename, "r");
                char line[200];
                iy=0;
                
                while (fgets(line, sizeof(line), bSiteList))
                {
                    bSite[iy]=atoi(line);
                    iy++;
                }
                
                fclose(bSiteList);
                
                bSiteTotal=iy;
                
                break;
                
            case 4: // use last site as only bSite
                
                bSiteTotal = 1;
                for(ib=0;ib<bSiteTotal;ib++)
                {
                    bSite[ib]=0;
                }
                bSite[0] = N-1;
                
                break;
                
        }
        
        //Warning for possible user error
        for (iy=0; iy<bSiteTotal;iy++)
        {
            if (bSite[iy] >= N)
            {
                printf("Warning! Bound site is located past end of polymer!");
                fflush(stdout);
            }
        }
        
        //for debugging - prints a list of the bSites
        if (TALKATIVE)
        {
            for (ib=0; ib<bSiteTotal; ib++)
            {
                printf("bSite[%ld] is %ld \n", ib, bSite[ib]);
                fflush(stdout);
            }
            
            
            printf("bSiteTotal = %ld \n", bSiteTotal);
            fflush(stdout);
        }

    }
    
    
    


}

/********************************************************************************************************/

