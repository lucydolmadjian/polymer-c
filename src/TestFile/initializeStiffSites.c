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
        
        
        
        //initializes phosiSites to 0 (none phosphorylated)
        for(ty=0;ty<iSiteTotal;ty++)
        {
            phosiSites[ty]=0;
        }
        
        
        if (TALKATIVE) printf("This is a string: %s\n", phosphorylatediSites);
        
        //read string and assign to double vector
        
        // 1 is occupied iSite (phosphorylated), 0 is unoccupied
        
        switch (testRun)
        {
            case 0:
                
                sscanf(phosphorylatediSites,"%lf %lf %lf %lf %lf %lf %lf", &phosiSites[0],&phosiSites[1],&phosiSites[2],&phosiSites[3], &phosiSites[4],&phosiSites[5],&phosiSites[6]);
                break;
                
            case 1:
                
                sscanf(phosphorylatediSites,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &phosiSites[0],&phosiSites[1],&phosiSites[2],&phosiSites[3], &phosiSites[4],&phosiSites[5],&phosiSites[6],&phosiSites[7],&phosiSites[8],&phosiSites[9],&phosiSites[10],&phosiSites[11],&phosiSites[12],&phosiSites[13]);
                break;
                
            case 2:
                
//                sscanf(phosphorylatediSites,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &phosiSites[0],&phosiSites[1],&phosiSites[2],&phosiSites[3], &phosiSites[4],&phosiSites[5],&phosiSites[6],&phosiSites[7],&phosiSites[8],&phosiSites[9],&phosiSites[10],&phosiSites[11],&phosiSites[12],&phosiSites[13]);
//                break;

                sscanf(phosphorylatediSites,"%lf %lf", &phosiSites[0], &phosiSites[1]);
                break;
                
        }
        
        
        if (TALKATIVE)
            for (iy=0;iy<iSiteTotal;iy++)
            {
                printf("phosiSites[ %ld ] =  %f\n",iy, phosiSites[iy]);
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
    

}

/********************************************************************************************************/

