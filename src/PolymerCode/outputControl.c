/*** Allard Group jun.allard@uci.edu                    ***/

void initializeSummary();
void finalizeSummary();
void dataRecording();

/*******************************************************************************/
//  GLOBAL VARIABLES for output control
/*******************************************************************************/

double reeBar_sum[NFILMAX], ree2Bar_sum[NFILMAX], rMBar_sum[NFILMAX], rM2Bar_sum[NFILMAX], rMiSiteBar_sum[NFILMAX][NMAX], rM2iSiteBar_sum[NFILMAX][NMAX];

long POcclude_sum[NFILMAX][NMAX], Prvec0_sum[NFILMAX][NMAX], POccludeBase_sum[NFILMAX], PDeliver_sum[NFILMAX][NMAX], PMembraneOcclude_sum[NFILMAX][NMAX],PMembraneSegmentOcclude_sum[NFILMAX][NMAX];

double reeBar[NFILMAX], ree2Bar[NFILMAX];
double POcclude[NFILMAX][NMAX], POccludeBase[NFILMAX];
double PDeliver[NFILMAX][NMAX], Prvec0[NFILMAX][NMAX],PMembraneOcclude[NFILMAX][NMAX],PMembraneSegmentOcclude[NFILMAX][NMAX];
double reeiSite[NFILMAX][NMAX], rMBar[NFILMAX], rM2Bar[NFILMAX], rMiSiteBar[NFILMAX][NMAX], rM2iSiteBar[NFILMAX][NMAX];

double distiSiteToLigand[NFILMAX][NMAX][NMAX], selfBind[NFILMAX][NMAX][NMAX], selfBindFraction[NFILMAX][NMAX][NMAX], localConcentration[NFILMAX][NMAX][NMAX];

double occupied[NFILMAX][NMAX];
double binSize[NFILMAX];
long binCurrent;

/********************************************************************************************************/
void initializeSummary()
{
    // summary variables
    for(nf=0;nf<NFil;nf++)
    {
        reeBar_sum[nf]   = 0;
        ree2Bar_sum[nf]  = 0;
        for(iy=0;iy<iSiteTotal[nf];iy++)
        {
            POcclude_sum[nf][iy]                = 0;
            Prvec0_sum[nf][iy]                  = 0;
            PMembraneOcclude_sum[nf][iy]        = 0;
        }
        POccludeBase_sum[nf] = 0;
    //    for(ib=0; ib<bSiteTotal[nf]; ib++)
    //    {
    //        PDeliver[nf][ib]=0;
    //    }
        rMBar_sum[nf]    = 0;
        rM2Bar_sum[nf]   = 0;
        
        //distance of iSites to membrane
        for (iy=0;iy<iSiteTotal[nf];iy++)
        {
            rMiSiteBar_sum[nf][iy] = 0;
            rM2iSiteBar_sum[nf][iy] = 0;
        }
    }
    
    for(nf=0;nf<NFil;nf++)
    {
        binSize[nf] = (double)(2*N[nf]) / NBINSPOLYMER;
    }
    
    if(MULTIPLE)
    {
        for(nf=0;nf<NFil;nf++)
        {
            for(iy=0; iy<iSiteTotal[nf]; iy++)
            {
                for(ib=0; ib<bSiteTotal[nf]; ib++)
                {
                    distiSiteToLigand[nf][iy][ib]  = 0;
                    selfBind[nf][iy][ib]           = 0;
                    selfBindFraction[nf][iy][ib]   = 0;
                    localConcentration[nf][iy][ib] = 0;
                }
            }
        }
    }

}

/********************************************************************************************************/
void finalizeSummary()
{
    // finalize summary statistics
    for(nf=0;nf<NFil;nf++)
    {
        reeBar[nf]  = reeBar_sum[nf]/(double)(nt-NTCHECK);
        ree2Bar[nf] = ree2Bar_sum[nf]/(double)(nt-NTCHECK);

        for(iy=0;iy<iSiteTotal[nf];iy++)
        {
            POcclude[nf][iy]         = (double)POcclude_sum[nf][iy]/(double)(nt-NTCHECK);
            Prvec0[nf][iy]           = (double)Prvec0_sum[nf][iy]/(4/3*PI*pow((double)N[nf]/(double)NBINS,3))/(double)(nt-NTCHECK);
            PMembraneOcclude[nf][iy] = (double)PMembraneOcclude_sum[nf][iy]/(double)(nt-NTCHECK);
            PMembraneSegmentOcclude[nf][iy] = (double)PMembraneSegmentOcclude_sum[nf][iy]/(double)(nt-NTCHECK);

        }
        
        POccludeBase[nf] = (double)POccludeBase_sum[nf]/(double)(nt-NTCHECK);
        
    //    for(ib=0;ib<bSiteTotal[nf];ib++)
    //    {
    //        PDeliver[nf][ib] = (double)PDeliver_sum[nf][ib]/(double)(nt-NTCHECK);
    //    }
        
        rMBar[nf] = rMBar_sum[nf]/(double)(nt-NTCHECK);
        rM2Bar[nf] = rM2Bar_sum[nf]/(double)(nt-NTCHECK);
        
        //average distance of iSites to membrane
        for (iy=0;iy<iSiteTotal[nf];iy++)
        {
            rMiSiteBar[nf][iy] = rMiSiteBar_sum[nf][iy]/(double)(nt-NTCHECK);
            rM2iSiteBar[nf][iy] = rM2iSiteBar_sum[nf][iy]/(double)(nt-NTCHECK);
        }
    }
    
    if(MULTIPLE)
    {
        for(nf=0;nf<NFil;nf++)
        {
            // determine local concentration
            for(iy=0;iy<iSiteTotal[nf];iy++)
            {
                for(ib=0;ib<bSiteTotal[nf];ib++)
                {
                    selfBindFraction[nf][iy][ib] = (double)selfBind[nf][iy][ib]/(double)(nt-NTCHECK);
                    localConcentration[nf][iy][ib] =  selfBindFraction[nf][iy][ib] / ((double) 4/3*pow(localConcCutoff,3)*PI);
                }
            }
        }
    }
    
    if (!verboseTF)
    {
        fList = fopen(listName, "a");
        
        for(nf=0;nf<NFil;nf++)
        {
            
            if (CD3ZETA)
            {
                    sscanf(occupiedSites,"%lf_%lf_%lf_%lf_%lf_%lf", &occupied[nf][0],&occupied[nf][1],&occupied[nf][2],&occupied[nf][3], &occupied[nf][4],&occupied[nf][5]);
                
                    for (iy=0; iy<iSiteTotal[nf];iy++)
                    {
                        fprintf(fList, "%lf ", occupied[nf][iy]);
                    }
                
                // eventually want this to depend on filament
                    fprintf(fList, "%s ", occupiedSitesNoSpace);
            }
            
            fprintf(fList, "%ld %ld %f %f %f %ld %f %f %f %e",
                    NFil,            // 1
                    N[nf],           // 2
                    irLigand,       // 3
                    brLigand,       // 4
                    Force,          // 5
                    nt,             // 6
                    ksStatistic[nf],    // 7
                    reeBar[nf],      // 8
                    ree2Bar[nf],     // 9
                    rMBar[nf]);      // 10
            
            for (iy=0;iy<iSiteTotal[nf];iy++)
            {
                fprintf(fList, " %ld %e %e %e %e",
                    iSite[nf][iy], //11 + 4*iBind
                    POcclude[nf][iy], //12 + 4*iBind
                    1-POcclude[nf][iy], //13 + 4*iBind
                    PMembraneOcclude[nf][iy], //14 +4*iBind
                    Prvec0[nf][iy]); //15 + 4*iBind
            
            }

            
            fprintf(fList, " %d %e %e", -1, POccludeBase[nf], 1-POccludeBase[nf]);
            
            for (ib=0;ib<bSiteTotal[nf];ib++)
            {
                fprintf(fList, " %ld",
                    bSite[nf][ib]);
                        //PDeliver[nf][ib]);
            }
            
            for (iy=0; iy<iSiteTotal[nf]; iy++)
            {
                for (ib=0; ib<bSiteTotal[nf]; ib++)
                {
                    fprintf(fList, " %f", selfBindFraction[nf][iy][ib]);
                }
            }
            
            for (iy=0; iy<iSiteTotal[nf]; iy++)
            {
                for (ib=0; ib<bSiteTotal[nf]; ib++)
                {
                    fprintf(fList, " %f", localConcentration[nf][iy][ib]);
                }
            }
            
            for(iy=0;iy<iSiteTotal[nf];iy++)
            {
                fprintf(fList, " %f", PMembraneSegmentOcclude[nf][iy]);
            }
            
            
            // print distribution of iSite Locations
            if(1)
            {
                for (iy=0;iy<iSiteTotal[nf];iy++)
                {
                    iSiteCurrent=iSite[nf][iy];
                    for (j=0;j<NBINSPOLYMER;j++)
                    {
                        fprintf(fList, " %ld", polymerLocationCounts[nf][iSiteCurrent][j]);
                    }
                }
            }
            
            // print distributions of end segments locations
            if(1)
            {
                for (j=0;j<NBINSPOLYMER;j++)
                {
                    Ncurrent = N[nf];
                    fprintf(fList, " %ld", polymerLocationCounts[nf][Ncurrent-1][j]);
                }
            }
            
            fprintf(fList, " %f", rM2Bar[nf]);
            
            for (iy=0;iy<iSiteTotal[nf];iy++)
            {
                fprintf(fList, " %f", rMiSiteBar[nf][iy]);
            }
            
            for (iy=0;iy<iSiteTotal[nf];iy++)
            {
                fprintf(fList, " %f", rM2iSiteBar[nf][iy]);
            }
            
            
        } // end printing data for each filament
        fprintf(fList, "\n");
        fclose(fList);
    }
}

/********************************************************************************************************/
// Prepare stuff and optionally write to file - this function is called each timestep
void dataRecording()
{
	
    // end-to-end distance
    for(nf=0;nf<NFil;nf++)
    {
        Ncurrent = N[nf];
        ree[nf]  = sqrt(r[nf][Ncurrent-1][0]*r[nf][Ncurrent-1][0]+
                        r[nf][Ncurrent-1][1]*r[nf][Ncurrent-1][1]+
                        r[nf][Ncurrent-1][2]*r[nf][Ncurrent-1][2]);
    }
	
    // distance from base to iSite
    for(nf=0;nf<NFil;nf++)
    {
        for(iy=0;iy<iSiteTotal[nf];iy++)
        {
            iSiteCurrent = iSite[nf][iy];
            reeiSite[nf][iy] = sqrt(r[nf][iSiteCurrent][0]*r[nf][iSiteCurrent][0] +
                                    r[nf][iSiteCurrent][1]*r[nf][iSiteCurrent][1] +
                                    r[nf][iSiteCurrent][2]*r[nf][iSiteCurrent][2]);
        }
    }

    // distance of tip to membrane
    for(nf=0;nf<NFil;nf++)
    {
        Ncurrent = N[nf];
        rM[nf] = r[nf][Ncurrent-1][2];
        rM2[nf] = r[nf][Ncurrent-1][2]*r[nf][Ncurrent-1][2];
    }
    
    //distance of iSites to membrane
    for(nf=0;nf<NFil;nf++)
    {
        for (iy=0;iy<iSiteTotal[nf];iy++)
        {
            iSiteCurrent = iSite[nf][iy];
            rMiSite[nf][iy] = r[nf][iSiteCurrent][2];
            rM2iSite[nf][iy] = r[nf][iSiteCurrent][2]*r[nf][iSiteCurrent][2];
        }
    }
	
    // height (max distance to membrane)
	if  (0)
    {
        for(nf=0;nf<NFil;nf++)
        {
            rH[nf] = 0;
            for(i=0;i<N[nf];i++)
                if (r[nf][i][2]>rH[nf])
                    rH[nf] = r[nf][i][2];
        }
    }


    // Verbose output: One line is written each iteration.
    if (verboseTF)
    {
        
        if ( (nt > NTCHECK && nt <= NTCHECK+4000) ) //only output 4000 runs, after initial transient
        {
        // output results to file
        fList = fopen(listName, "a");
        for(nf=0;nf<NFil;nf++)
        {
            fprintf(fList, "%ld %f %f %f %f %f %f %f %f %f %ld",
                    nt,               // 1
                    ree[nf],          // 2
                    rM[nf],           // 3
                    rH[nf],           // 4
                    E,                // 5
                    dChi[0],          // 6
                    dChi[1],          // 7
                    rate[0],          // 8
                    rate[1],          // 9
                    ksStatistic[nf],      // 10
                    constraintProposalsTotal);// 11
            
                 if (VISUALIZE)
                 {
                     
                    fprintf(fList, " %f %f %f", rBase[nf][0],rBase[nf][1],rBase[nf][2]);
                     
                    for (i=0;i<N[nf];i++)
                    {
                        fprintf(fList, " %f %f %f", r[nf][i][0],r[nf][i][1],r[nf][i][2]);
                    }
                     
                    for (i=0;i<iSiteTotal[nf];i++)
                    {
                        fprintf(fList, " %f %f %f", iLigandCenter[nf][i][0], iLigandCenter[nf][i][1], iLigandCenter[nf][i][2]);
                    }
                    if(BASEBOUND)
                    {
                        for (i=0;i<3;i++)
                        {
                            fprintf(fList," %f", baseCenter[i]);
                        }
                    }
                     
                 }
            
                for(iy=0;iy<iSiteTotal[nf];iy++)
                {
                    fprintf(fList, " %ld",stericOcclusion[nf][iy]);
                }
            
                for(iy=0;iy<iSiteTotal[nf];iy++)
                {
                    fprintf(fList, " %ld",membraneOcclusion[nf][iy]);
                }
            
                for(iy=0;iy<iSiteTotal[nf];iy++)
                {
                    fprintf(fList, " %ld",membraneAndSegmentOcclusion[nf][iy]);
                }
            
                fprintf(fList, " %ld", stericOcclusionBase[nf]);
            
    //        for(ib=0;ib<bSiteTotal[nf];ib++)
    //        {
    //            fprintf(fList, " %ld", boundToBaseDeliver[nf][ib]);
    //        }

        } // end of filament loop
            
            fprintf(fList, "\n");

            fclose(fList);
        }
    } // finished verbose output
    
    // Note:  eliminates initial transient
    if (nt>NTCHECK)
    {
        // Summary statistics (these will be written at the end of the run)
        for(nf=0;nf<NFil;nf++)
        {
            reeBar_sum[nf]   += ree[nf];
            ree2Bar_sum[nf]  += ree[nf]*ree[nf];
        }
        
        for(nf=0;nf<NFil;nf++)
        {
            for(iy=0;iy<iSiteTotal[nf];iy++)
            {
                POcclude_sum[nf][iy]         += (long)(stericOcclusion[nf][iy]>0);
                Prvec0_sum[nf][iy]           += (long)(reeiSite[nf][iy] < (double)N[nf]/(double)NBINS);
                PMembraneOcclude_sum[nf][iy] += (long)(membraneOcclusion[nf][iy]>0);
                PMembraneSegmentOcclude_sum[nf][iy] += (long)(membraneAndSegmentOcclusion[nf][iy]>0);
                
            }
        }
        for(nf=0;nf<NFil;nf++)
        {
            POccludeBase_sum[nf]  += (long)(stericOcclusionBase[nf]>0);
            //PDeliver_sum[nf][ib] += (long)(boundToBaseDeliver[nf]>0);
            rMBar_sum[nf]    += rM[nf];
            rM2Bar_sum[nf]   += rM2[nf];
        }
        
        for(nf=0;nf<NFil;nf++)
        {
            for (iy=0;iy<iSiteTotal[nf];iy++)
            {
                rMiSiteBar_sum[nf][iy]    += rMiSite[nf][iy];
                rM2iSiteBar_sum[nf][iy]   += rM2iSite[nf][iy];
            }
        }

        // update bins for KS test (fabs(rM)+ree will never be larger than 2N, so use 2N to normalize)
        for(nf=0;nf<NFil;nf++)
        {
            convergenceVariableCounts[nf][(long)floor(NBINS*(fabs(rM[nf])+ree[nf])/(2*N[nf]))]++;
        }
        
        // Distributions for polymer location
        for(nf=0;nf<NFil;nf++)
        {
            for(i=0;i<N[nf];i++)
            {
                binCurrent = (long)floor(((double)r[nf][i][2]+N[nf])/binSize[nf])
                polymerLocationCounts[nf][i][binCurrent]++;
            }
        }
        
        
        
        /* LOCAL CONCENTRATION */
        
        // initialize local concentration arrays
        
        if (MULTIPLE)
        {
            for(nf=0;nf<NFil;nf++)
            {
                for(iy=0; iy<iSiteTotal[nf];iy++)
                {
                    for(ib=0;ib<bSiteTotal[nf];ib++)
                    {
                        distiSiteToLigand[nf][iy][ib]  = 0;
                    }
                }
            }
            // calculate distance between each ligand and iSite center
            for(nf=0;nf<NFil;nf++)
            {
                for(iy=0;iy<iSiteTotal[nf];iy++)
                {
                    for(ib=0;ib<bSiteTotal[nf];ib++)
                    {
                        distiSiteToLigand[nf][iy][ib] = sqrt((iLigandCenter[nf][iy][0]-bLigandCenter[nf][ib][0])*
                                                             (iLigandCenter[nf][iy][0]-bLigandCenter[nf][ib][0])+
                                                             (iLigandCenter[nf][iy][1]-bLigandCenter[nf][ib][1])*
                                                             (iLigandCenter[nf][iy][1]-bLigandCenter[nf][ib][1])+
                                                             (iLigandCenter[nf][iy][2]-bLigandCenter[nf][ib][2])*
                                                             (iLigandCenter[nf][iy][2]-bLigandCenter[nf][ib][2]));
                    }
                }
            }
            
            // determine if distance is less than given cutoff AND iSite is not occluded
            for(nf=0;nf<NFil;nf++)
            {
                for(iy=0;iy<iSiteTotal[nf];iy++)
                {
                    // if iSite is not occluded
                    if(membraneAndSegmentOcclusion[nf][iy] == 0)
                    {
                        // test how far away each ligand is
                        for(ib=0;ib<bSiteTotal[nf];ib++)
                        {
                            if(distiSiteToLigand[nf][iy][ib] < localConcCutoff)
                            {
                                // if ligand is close enough (distance less than cutoff), ligand can bind
                                selfBind[nf][iy][ib]++;
                            }
                        }
                    }
                }
            }
            
        }
                                     
    }
	return;
	
}

