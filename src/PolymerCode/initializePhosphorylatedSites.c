/*** Allard Group jun.allard@uci.edu                    ***/

void initializePhosphorylatedSites();


/*******************************************************************************/
//  GLOBAL VARIABLES for output control
/*******************************************************************************/

int totalPhosphorylated;

double phosiSites[NMAX];


//
/********************************************************************************************************/
void initializePhosphorylatedSites()
{
        
        
        //initializes phoshorylatedSites to 0 (none phosphorylated)
        for(ty=0;ty<iSiteTotal;ty++)
        {
            phosiSites[ty]=0;
        }
     
        if (TALKATIVE) printf("This is a string: %s\n", occupiedSites);
        
        //read string and assign to double vector
        
        // 1 is occupied iSite (phosphorylated), 0 is unoccupied

        //read string or read and interpret which of CD3Zeta Mouse tyrosines are phosphorylated
    
        sscanf(occupiedSites,"%lf_%lf_%lf_%lf_%lf_%lf", &phosiSites[0],&phosiSites[1],&phosiSites[2],&phosiSites[3], &phosiSites[4],&phosiSites[5]);


        //print which of the iSites are phosphorylated
        if (TALKATIVE)
            for (iy=0;iy<iSiteTotal;iy++)
            {
                printf("phosiSites[ %ld ] =  %f\n",iy, phosiSites[iy]);
				fflush(stdout);
            }
        
        //initializes stiffened rods to 0 (none stiff)
        for(i=0;i<N;i++)
        {
            PhosphorylatedSites[i] =0;
        }
        
        
        
        //Phosphorylate sites
        for(ty=0;ty<iSiteTotal;ty++)
        {
            if(phosiSites[ty]==1) //might want to check the truth value on this - equals for double?
            {
                for(i=(iSite[ty]-PhosElectroRange);i<(iSite[ty]+PhosElectroRange+1);i++)
                {
                    PhosphorylatedSites[i]=1; //set that joint to "phosphorylated"
                }
            }
        }
    
        if (TALKATIVE)
        {
            totalPhosphorylated = 0;
            for (i=0;i<N;i++)
            {
                if (PhosphorylatedSites[i]==1)
                {
                    //printf("PhosphorylatedSites[ %ld ] =  %f\n",i, PhosphorylatedSites[i]);
                    //fflush(stdout);
                    totalPhosphorylated++;
                }
            }
            printf("Total Phosphorylated: %d\n", totalPhosphorylated);
            fflush(stdout);

        }
    
    

}

/********************************************************************************************************/

