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
//        case 0:  // iSites initialized for human CD3Zeta-Chain
//            
//            iSiteTotal = 7;
//            
//            
//            for (iy=0;iy<iSiteTotal;iy++) //initializes iSite array
//            {
//                iSite[iy]=0;
//            }
//            
//            //specify iSite locations by rod number (located at n-1) (i.e. if you have a polymer of N=50, and want an iSite at the 36th rod, input iSite[]=35)
//            
//            iSite[0]=42;
//            iSite[1]=50;
//            iSite[2]=61;
//            iSite[3]=89;
//            iSite[4]=101;
//            iSite[5]=120;
//            iSite[6]=131;
//            break;
            
        case 0:  // iSites initialized for human CD3Zeta-Chain

            iSiteTotal = 6;


            for (iy=0;iy<iSiteTotal;iy++) //initializes iSite array
            {
                iSite[iy]=0;
            }

            //specify iSite locations by rod number (located at n-1) (i.e. if you have a polymer of N=50, and want an iSite at the 36th rod, input iSite[]=35)

            iSite[0]=50;
            iSite[1]=61;
            iSite[2]=89;
            iSite[3]=101;
            iSite[4]=120;
            iSite[5]=131;
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
            
            iSiteTotal = 3;
            
            for(iy=0;iy<iSiteTotal;iy++)
            {
                iSite[iy]=0;
            }
            
            iSite[0]=0;
            iSite[1]=1;
            iSite[2]=2;
//            iSite[3]=3;
//            iSite[4]=4;
//            iSite[5]=5;
//            iSite[6]=6;
//            iSite[7]=7;
//            iSite[8]=8;
//            iSite[9]=9;
//            iSite[10]=10;
//            iSite[11]=11;
//            iSite[12]=12;
//            iSite[13]=13;


            break;
            
        case 3: //do nothing, use command line input
            
            break;
            
        case 4: //input iSites from file
            
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
    
    /********* INITIALIZE BOUND ISITES *******************/
    
    if (MULTIPLE) //if looking at multiple binding (i.e. MULTIPLE set to 1 in driveM)
    {
        switch (bSiteCommand)
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
                
                sscanf(occupiedSites,"%lf %lf %lf %lf %lf %lf", &iSiteOccupied[0],&iSiteOccupied[1],&iSiteOccupied[2],&iSiteOccupied[3], &iSiteOccupied[4],&iSiteOccupied[5]);
                
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
                
                //for debugging
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
                
                break;
                
            case 3:
                
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
                
        }
        
        //for debugging - prints a list of the iSites
        
        if (TALKATIVE)
        {
            for (iy=0;iy<bSiteTotal;iy++)
            {
                printf("bSite: %ld\n", bSite[iy]);
                fflush(stdout);
            }
            
            printf("bSiteTotal: %ld\n", bSiteTotal);
            fflush(stdout);
        }

    }
    
    
    


}

/********************************************************************************************************/

