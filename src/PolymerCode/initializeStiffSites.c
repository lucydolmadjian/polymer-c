/*** Allard Group jun.allard@uci.edu                    ***/

void initializeStiffSites();


/*******************************************************************************/
//  GLOBAL VARIABLES for output control
/*******************************************************************************/

//
/********************************************************************************************************/
void initializeStiffSites()
{
        
        
        //initializes stiffiSites to 0 (none phosphorylated)
        for(ty=0;ty<iSiteTotal;ty++)
        {
            stiffiSites[ty]=0;
        }
        
        
        printf("This is a string: %s\n", occupiedSites);
        
     
        if (TALKATIVE) printf("This is a string: %s\n", occupiedSites);
        
        //read string and assign to double vector
        
        // 1 is occupied iSite (phosphorylated), 0 is unoccupied
        
        switch (stiffCase) //do I want this to be testRun? or another variable?
        {
            case 0:
                
                sscanf(occupiedSites,"%lf_%lf_%lf_%lf_%lf_%lf", &stiffiSites[0],&stiffiSites[1],&stiffiSites[2],&stiffiSites[3], &stiffiSites[4],&stiffiSites[5]);
                break;
                
            case 1:
                
                sscanf(occupiedSites,"%lf_%lf_%lf_%lf_%lf_%lf_%lf_%lf_%lf_%lf_%lf_%lf_%lf_%lf", &stiffiSites[0],&stiffiSites[1],&stiffiSites[2],&stiffiSites[3], &stiffiSites[4],&stiffiSites[5],&stiffiSites[6],&stiffiSites[7],&stiffiSites[8],&stiffiSites[9],&stiffiSites[10],&stiffiSites[11],&stiffiSites[12],&stiffiSites[13]);
                break;
                
            case 2:
                
                sscanf(occupiedSites,"%lf_%lf_%lf_%lf_%lf_%lf_%lf_%lf_%lf_%lf_%lf_%lf_%lf_%lf", &stiffiSites[0],&stiffiSites[1],&stiffiSites[2],&stiffiSites[3], &stiffiSites[4],&stiffiSites[5],&stiffiSites[6],&stiffiSites[7],&stiffiSites[8],&stiffiSites[9],&stiffiSites[10],&stiffiSites[11],&stiffiSites[12],&stiffiSites[13]);
                break;

//                sscanf(occupiedSites,"%lf_%lf", &stiffiSites[0], &stiffiSites[1]);
                //break;
                
        }
        
        
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
        
        
        
        //Stiffen segments
        for(ty=0;ty<iSiteTotal;ty++)
        {
            if(stiffiSites[ty]==1) //might want to check the truth value on this - equals for double?
            {
                if(iSite[ty]-StiffenRange >= 0)
                {
                    stiffStart=iSite[ty]-StiffenRange;
                }
                else
                {
                    stiffStart=0;
                }
                
                if(iSite[ty]+StiffenRange+1 >= N)
                {
                    stiffEnd=N;
                }
                else
                {
                    stiffEnd=iSite[ty]+StiffenRange+1;
                }

                for(i=stiffStart;i<stiffEnd;i++)
                {
                    StiffSites[i]=1; //set that joint to "stiff"
                }
            }
        }
    
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
                printf("Error! Completely stiff!\n");
                fflush(stdout);
                
                exit(0);
                
            }
        }
    
    

}

/********************************************************************************************************/

