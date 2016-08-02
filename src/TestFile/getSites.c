/*** Allard Group jun.allard@uci.edu                    ***/

void getSites();


/*******************************************************************************/
//  GLOBAL VARIABLES for output control
/*******************************************************************************/


/********************************************************************************************************/
void getSites()
{
    /********* INITIALIZE ISITES *******************/
    
    
        for (iy=0;iy<iSiteTotal;iy++)
        {
            iSite[iy]=0;
        }
    
    switch (commandiSites)
    {
    case 0:
    
    switch (testRun)
    {
        case 0:  // iSites initialized for human CD3Zeta-Chain
            
            iSiteTotal = 7;
            
            
            for (iy=0;iy<iSiteTotal;iy++)
            {
                iSite[iy]=0;
            }
            
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
            
            
        case 2: //test case 2 - stiffen none, but test all iSites, make Ratio half of Ratio for case 1
            
            iSiteTotal = 7;
            
            for (iy=0;iy<iSiteTotal;iy++)
            {
                iSite[iy]=0;
            }
            
            iSite[0]=0;
            iSite[1]=1;
            iSite[2]=2;
            iSite[3]=3;
            iSite[4]=4;
            iSite[5]=5;
            iSite[6]=6;
            break;
            
        case 3:
            
            iSiteTotal = 4;
            
            for(iy=0;iy<iSiteTotal;iy++)
            {
                iSite[iy]=0;
            }
            
            iSite[0]=10;
            iSite[1]=11;
            iSite[2]=23;
            iSite[3]=40;

            break;

    
    
    
                case 4:
    
                    for (iy=0;iy<iSiteTotal;iy++)
                    {
                        iSite[iy]=0;
                    }
    
                printf("Total iSites: %ld", iSiteTotal);
    
                //for debugging
                for(iy=0;iy<iSiteTotal;iy++)
            {
                printf("iSite[%ld] = %ld", iy, iSite[iy]);
            }
    
                    //char input[] = iSiteLocations;
                    printf("I want to split this into tokens: %s", input);
                    char* strArray[NMAX];
                    char *token = strtok(input, " ");
    
                    //for(int j = 0; j<NMAX;j++)
                    //{
                    //    strArray[j] = new char[4];
                    //}
    
                    while(token != NULL)
                    {
                        strcpy(strArray[st],token);
                        printf("This is the next token: %s\n",token); //for debugging
                        token = strtok(NULL, " ");
                        st++;
                    }
    
                    //for debugging
    
                    if (iSiteTotal!=st)
                    {
                        printf("Warning! iSite Total is %ld but Number of iSites in String is %ld ", iSiteTotal, st);
                    }
    
                    //reassign strings as doubles
                    for(iy=0;iy<st;iy++)
                    {
                        iSite[iy]=atof(strArray[iy]);
                    }
    
                    //for debugging
                    for(iy=0;iy<iSiteTotal;iy++)
                    {
                        printf("iSite[%ld] = %ld", iy, iSite[iy]);
                    }
    
                break;
        }
    
    /********* INITIALIZE BOUND ISITES *******************/
    
    if (MULTIPLE)
    {
        
        
        //add switch to be able to input 101010 etc to specify which iSites are bound/unbound
        
        //switch () //add more cases later
        //{
        //case 0: // arbitrary subset are occupied
        
        bSiteTotal = 4; //total number of iSites bound
        
        bSite[0]=10;
        bSite[1]=10;
        bSite[2]=40;
        bSite[3]=40;

        

        
        
        //currently identifying by location
        //what is the best way to do this - identify by location or identify by iSite number?
        //pro for location - can specify locations other than iSites to be bound - but then might want to change name
        //}
    }


}

/********************************************************************************************************/

