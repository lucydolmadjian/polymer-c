/*** Allard Group jun.allard@uci.edu                    ***/

void initializeStiffSites();


/*******************************************************************************/
//  GLOBAL VARIABLES for output control
/*******************************************************************************/


/********************************************************************************************************/
void initializeStiffSites()
{

        
        //initializes phosiSites to 0 (none phosphorylated)
        for(ty=0;ty<iSiteTotal;ty++)
        {
            phosiSites[ty]=0;
        }
        
        
        printf("This is a string: %s\n", occupiedSites);
        
     
        if (TALKATIVE) printf("This is a string: %s\n", occupiedSites);
        
        //read string and assign to double vector
        
        // 1 is occupied iSite (phosphorylated), 0 is unoccupied
        
        switch (stiffCase) //do I want this to be testRun? or another variable?
        {
            case 0:
                
                sscanf(occupiedSites,"%lf %lf %lf %lf %lf %lf %lf", &phosiSites[0],&phosiSites[1],&phosiSites[2],&phosiSites[3], &phosiSites[4],&phosiSites[5],&phosiSites[6]);
                break;
                
            case 1:
                
                sscanf(occupiedSites,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &phosiSites[0],&phosiSites[1],&phosiSites[2],&phosiSites[3], &phosiSites[4],&phosiSites[5],&phosiSites[6],&phosiSites[7],&phosiSites[8],&phosiSites[9],&phosiSites[10],&phosiSites[11],&phosiSites[12],&phosiSites[13]);
                break;
                
            case 2:
                
                sscanf(occupiedSites,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &phosiSites[0],&phosiSites[1],&phosiSites[2],&phosiSites[3], &phosiSites[4],&phosiSites[5],&phosiSites[6],&phosiSites[7],&phosiSites[8],&phosiSites[9],&phosiSites[10],&phosiSites[11],&phosiSites[12],&phosiSites[13]);
                break;

//                sscanf(occupiedSites,"%lf %lf", &phosiSites[0], &phosiSites[1]);
                //break;
                
        }
        
        
        if (TALKATIVE)
            for (iy=0;iy<iSiteTotal;iy++)
            {
                printf("phosiSites[ %ld ] =  %f\n",iy, phosiSites[iy]);
				fflush(stdout);
            }
        
        //initializes stiffened rods to 0 (none stiff)
        for(i=0;i<N;i++)
        {
            Stiff[i] =0;
        }
        
        
        
        //Stiffen segments
        for(ty=0;ty<iSiteTotal;ty++)
        {
            if(phosiSites[ty]==1) //might want to check the truth value on this - equals for double?
            {
                for(i=(iSite[ty]-StiffenRange);i<(iSite[ty]+StiffenRange+1);i++) //can I even put this stuff in a for loop?
                {
                    Stiff[i]=1; //set that joint to "stiff"
                }
            }
        }
    
    
        if (TALKATIVE)
        {
            for (i=0;i<N;i++)
            {
                if (Stiff[i]==1)
                {
                    printf("Stiffen[ %ld ] =  %f\n",i, Stiff[i]);
                    fflush(stdout);
                }
            }
        }
    

}

/********************************************************************************************************/

