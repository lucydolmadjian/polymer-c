/*** Allard Group jun.allard@uci.edu                    ***/

void initializeSummary();
void finalizeSummary();
void dataRecording();

/*******************************************************************************/
//  GLOBAL VARIABLES for output control
/*******************************************************************************/

double reeBar_sum, ree2Bar_sum, rMBar_sum;
long POcclude_sum[NMAX], Prvec0_sum[NMAX], POccludeBase_sum, PDeliver_sum[NMAX];
double reeBar, ree2Bar, POcclude[NMAX], POccludeBase, PDeliver[NMAX], Prvec0[NMAX], reeiSite[NMAX], rMBar;

/********************************************************************************************************/
void initializeSummary()
{
    // summary variables
    reeBar_sum   = 0;
    ree2Bar_sum  = 0;
    for(iy=0;iy<iSiteTotal;iy++)
    {
        POcclude_sum[iy] = 0;
        Prvec0_sum[iy]   = 0;
    }
    POccludeBase_sum = 0;
//    for(ib=0; ib<bSiteTotal; ib++)
//    {
//        PDeliver[ib]=0;
//    }
    rMBar_sum    = 0;

}

/********************************************************************************************************/
void finalizeSummary()
{
    // finalize summary statistics
    reeBar = reeBar_sum/(double)(nt-NTCHECK);
    ree2Bar = ree2Bar_sum/(double)(nt-NTCHECK);

	
    for(iy=0;iy<iSiteTotal;iy++)
    {
        POcclude[iy] = (double)POcclude_sum[iy]/(double)(nt-NTCHECK);
        Prvec0[iy] = (double)Prvec0_sum[iy]/(4/3*PI*pow((double)N/(double)NBINS,3))/(double)(nt-NTCHECK);
    }
    
    POccludeBase = (double)POccludeBase_sum/(double)(nt-NTCHECK);
    
//    for(ib=0;ib<bSiteTotal;ib++)
//    {
//        PDeliver[ib] = (double)PDeliver_sum[ib]/(double)(nt-NTCHECK);
//    }
    
    rMBar = rMBar_sum/(double)(nt-NTCHECK);
    
    if (!verboseTF)
    {
        fList = fopen(listName, "a");
        
        
        if (CD3ZETA)
        {
            fprintf(fList, "%s %s ",occupiedSites,
                        occupiedSitesNoSpace);
        }
        
        fprintf(fList, "%ld %f %f %f %ld %f %f %f %e",
                N,           // 1
                irLigand,    // 2
                brLigand,    // 3
                Force,       // 4
                nt,          // 5
                ksStatistic, // 6
                reeBar,      // 7
                ree2Bar,     // 8
                rMBar);      // 9
        
        for (iy=0;iy<iSiteTotal;iy++)
        {
            fprintf(fList, " %ld %e %e %e",
                iSite[iy], //10 + 4*iBind
                POcclude[iy], //11 + 4*iBind
                1-POcclude[iy], //12 + 4*iBind
                Prvec0[iy]); //13 + 4*iBind
        
        }
        
        fprintf(fList, " %d %e %e", -1, POccludeBase, 1-POccludeBase);
        
        for (ib=0;ib<bSiteTotal;ib++)
        {
            fprintf(fList, " %ld",
                bSite[ib]);
                    //PDeliver[ib]);
        }
        
        fprintf(fList, "\n");
        fclose(fList);
    }
}

/********************************************************************************************************/
// Prepare stuff and optionally write to file - this function is called each timestep
void dataRecording()
{
	
    // end-to-end distance
	ree      = sqrt(  r[N-1][0]*r[N-1][0]       + r[N-1][1]*r[N-1][1]     + r[N-1][2]*r[N-1][2]);
	
    // distance from base to iSite
    for(iy=0;iy<iSiteTotal;iy++)
    {
        iSiteCurrent = iSite[iy];
        reeiSite[iy] = sqrt(r[iSiteCurrent][0]*r[iSiteCurrent][0] + r[iSiteCurrent][1]*r[iSiteCurrent][1] + r[iSiteCurrent][2]*r[iSiteCurrent][2]);
	}

    // distance of tip to membrane
	rM = r[N-1][2];
	
    // height (max distance to membrane)
	if  (0)
    {
        rH = 0;
        for(i=0;i<N;i++)
            if (r[i][2]>rH)
                rH = r[i][2];
    }
	
    // Verbose output: One line is written each iteration.
    if (verboseTF)
    {
        
        if ( (nt % 100) == 0) //only output every 100 time steps
        {
        // output results to file
        fList = fopen(listName, "a");
        fprintf(fList, "%ld %f %f %f %f %f %f %f %f %f %ld",
                nt,               // 1
                ree,              // 2
                rM,               // 3
                rH,               // 4
                E,                // 5
                dChi[0],          // 6
                dChi[1],          // 7
                rate[0],          // 8
                rate[1],          // 9
                ksStatistic,      // 10
        		constraintProposalsTotal);// 11
            
             if (VISUALIZE)
             {
                for (i=0;i<N;i++)
                {
                    fprintf(fList, " %f %f %f", r[i][0],r[i][1],r[i][2]);
                    fprintf(fList, " %f %f %f", phi[i],theta[i],psi[i]);
                }
             }
        
            for(iy=0;iy<iSiteTotal;iy++)
            {
                fprintf(fList, " %ld",stericOcclusion[iy]);
            }
        
            fprintf(fList, " %ld", stericOcclusionBase);
        
//        for(ib=0;ib<bSiteTotal;ib++)
//        {
//            fprintf(fList, " %ld", boundToBaseDeliver[ib]);
//        }
        
            fprintf(fList, "\n");

            fclose(fList);
        }
    } // finished verbose output
    
    if (nt>NTCHECK)
    {
        // Summary statistics (these will be written at the end of the run)
        reeBar_sum   += ree;
        ree2Bar_sum  += ree*ree;
        
        for(iy=0;iy<iSiteTotal;iy++)
        {
			POcclude_sum[iy] += (long)(stericOcclusion[iy]>0);
			Prvec0_sum[iy]   += (long)(reeiSite[iy] < (double)N/(double)NBINS);
			
        }
        POccludeBase_sum += (long)(stericOcclusionBase>0);
        //PDeliver_sum[ib] += (long)(boundToBaseDeliver>0);
        rMBar_sum    += rM;

        // update bins for KS test
		rMCounts[(long)floor(NBINS*fabs(rM)/N)]++;
    }
	return;
	
}

