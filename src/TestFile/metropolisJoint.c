/*** Allard Group jun.allard@uci.edu                    ***/

void metropolisJoint();
void stationarity();
void appendBins();
void rotate(double *tIn, double *e1In, double *e2In, double *tOut, double *e1Out, double *e2Out, double phiHere, double thetaHere,double psiHere);

void metropolisJoint()
{

    
//    /********* INITIALIZE ISITES AND BSITES *******************/

        getSites();
    
    
        /********* STIFFEN SEGMENTS *******************/
    
    if (StiffenRange > -1) //stiffen only if StiffenRange is 0 or greater
    {
        initializeStiffSites();
    }
    
    /********* INITIALIZE CONFIGURATION *******************/

    //initalize base
	tBase[0]=0;
	tBase[1]=0;
	tBase[2]=1;
	rBase[0]=0;
	rBase[1]=0;
	rBase[2]=0;
	e1Base[0]=1;
	e1Base[1]=0;
	e1Base[2]=0;
	e2Base[0]=0;
	e2Base[1]=1;
	e2Base[2]=0;
	
	// initial configuration
	for(i=0;i<N;i++)
	{
		phi[i]   = 0; 
		theta[i] = 0;
		psi[i]   = 0;
		
		r[i][0]  = 0; 
		r[i][1]  = 0;
		r[i][2]  = i+1;
		
		t[i][0]  = 0;
		t[i][1]  = 0;
		t[i][2]  = 1;
        
        e1[i][0] = 1;  
        e1[i][1] = 0;
        e1[i][2] = 0;
        
        e2[i][0] = 0;
        e2[i][1] = 1;
        e2[i][2] = 0;
        
	}
    
	
	convergedTF=0;
	nt=0;
    E = INF;

    // stuff needed for automatic perturbation size adaptation
	for(iParam=0;iParam<2;iParam++)
	{
		dChi[iParam] = DCHIINIT;
		proposals[iParam] = 0;
		accepts[iParam] = 0;
	}
    
    // stuff needed for automatic stationarity (convergence) check
	ntNextStationarityCheck = 3*NTCHECK;
	for(iBin=0;iBin<NBINS;iBin++)
		rMCounts[iBin] = 0;
	
    // summary variables
    initializeSummary();

    /********* BEGIN LOOP THROUGH ITERATIONS! *******************/
    //for (m=0;m<100;m++)
	while(!convergedTF && nt < NTMAX) // Time loop!
    {

        // adapt step size
        if (!(nt % NTADAPT))
        {
            rate[0] = (double)accepts[0]/(double)proposals[0];
            rate[1] = (double)accepts[1]/(double)proposals[1];
            
            //printf("current accepts/proposals: %d/%d=%f, %d/%d=%f\n", accepts[0], proposals[0], rate[0], accepts[1], proposals[1], rate[1]);
            
            for(iParam=0;iParam<2;iParam++)
            {
                accepts[iParam] = 0;
                proposals[iParam] = 0;
                if (nt < NTCHECK)
                {
                    if (rate[iParam] > 0.5 || rate[iParam] < 0.3)
                    {
                        dChi[iParam] = dChi[iParam]*rate[iParam]/0.44;
                        if (dChi[iParam] > PI/4)
                            dChi[iParam] = PI/4;
                        if (dChi[iParam] < DCHIMIN)
                            dChi[iParam] = DCHIMIN;
                    }
                }
            }
            
        } // done adapting step

        /********* OUTLINE OF ALGORITHM *******************/
        // 1. We create a new proposal configuration and
        // then decide:
        // 2. whether the proposal configuration satisfies any constraints (e.g., floor?), and
        // 3. whether the proposal configuration is accepted by Metropolis test (always accept for a freely-jointed chain)
        // After the configuration update, we:
        // 4. collect any data of interest (e.g., ree; whether there is steric occlusion), and write it to file
        // 5. test for stationarity (convergence)
        // 6. increment the iteration and repeat
		
        /********* 1. Generate proposal configuration *******************/

		for(i=0;i<N;i++)
		{
			phiPropose[i]   = phi[i];
			thetaPropose[i] = theta[i];
			psiPropose[i]   = psi[i];
		}
		
		
		iPropose = floor(N*TWISTER); // Initialize. This is the joint we will adjust this time.
        
        //set joints to stiff based on which iSites are occupied and the stiffness range
        if (StiffenRange > -1) //stiffen only if StiffenRange is 0 or greater
        {
            while(Stiff[iPropose]==1) //Test if proposed joint is stiff.
            {
                iPropose = floor(N*TWISTER); //If stiff, propose new joint until propose one not stiff.
            }
        }
        
		if(iPropose==0)
		{
			dChiHere = dChi[0];
			proposals[0]++;
		}
		else 
		{
			dChiHere = dChi[1];
			proposals[1]++;
		}
		
        constraintProposalsTotal = 0;
        long rounds = 0;
		constraintSatisfiedTF = 0;
        //for(j=0;j<5;j++)
		while(!constraintSatisfiedTF && constraintProposalsTotal < CPMAX) // keep generating proposal configurations until we get one that satisfies the constraint
        {
            iChi = floor(3*TWISTER);
				
            //printf("iChi=%d, iPropose=%d\n", iChi, iPropose);
				
            // propose new configuration angles
            // We can use normal (Gaussian) random variables using the Box-Muller method (see eg wikipedia)
            // or uniform random variable
            switch (iChi)
            {
                case 0:
                    //phiPropose[iPropose]   = phiPropose[iPropose]   + dChiHere*sqrt(-2*log(TWISTER))*cos(2*PI*TWISTER);
                    phiPropose[iPropose]   = phiPropose[iPropose]   + dChiHere*(2*TWISTER-1);
                    break;
                case 1:
                    //thetaPropose[iPropose] = thetaPropose[iPropose] + dChiHere*sqrt(-2*log(TWISTER))*cos(2*PI*TWISTER);
                    thetaPropose[iPropose] = thetaPropose[iPropose] + dChiHere*(2*TWISTER-1);
                    break;
                case 2:
                    //psiPropose[iPropose]   = psiPropose[iPropose]   + dChiHere*sqrt(-2*log(TWISTER))*cos(2*PI*TWISTER);
                    psiPropose[iPropose]   = psiPropose[iPropose]   + dChiHere*(2*TWISTER-1);
                    break;
            }

            // -- translate to proposal configuration --
            // rotate base
            rotate(&tBase[0], &e1Base[0], &e2Base[0], &tPropose[0][0], &e1Propose[0][0], &e2Propose[0][0], phiPropose[0], thetaPropose[0], psiPropose[0]);
            
        
//            printf("This is e1Propose: \n 1: %f\n 2: %f\n 3:%f\n", e1Propose[0][0],e1Propose[0][1],e1Propose[0][2]);
            
            
            for(ix=0;ix<3;ix++)
                rPropose[0][ix] = rBase[ix] + tPropose[0][ix];
				
            // for i<iPropose, proposal configuration is the same as current configuration
            for(i=1;i<iPropose;i++)
            {
                rPropose[i][0] = r[i][0];
                rPropose[i][1] = r[i][1];
                rPropose[i][2] = r[i][2];
                
                tPropose[i][0] = t[i][0];
                tPropose[i][1] = t[i][1];
                tPropose[i][2] = t[i][2];
                
                e1Propose[i][0] = e1[i][0];
                e1Propose[i][1] = e1[i][1];
                e1Propose[i][2] = e1[i][2];
                
                e2Propose[i][0] = e2[i][0];
                e2Propose[i][1] = e2[i][1];
                e2Propose[i][2] = e2[i][2];
            }
            
            // rotate all segments above (and including) iPropose
            if (iPropose==0)
                iStart=1;
            else
                iStart=iPropose;
            for(i=iStart;i<N;i++)
            {
                rotate(&tPropose[i-1][0], &e1Propose[i-1][0], &e2Propose[i-1][0], &tPropose[i][0], &e1Propose[i][0], &e2Propose[i][0],
                    phiPropose[i], thetaPropose[i], psiPropose[i]);
                for (ix=0;ix<3;ix++)
                    rPropose[i][ix] = rPropose[i-1][ix] + tPropose[i][ix];
            }
            
            //how does this work?  only returns first component of t, e1, e2? where are the other components?
//            

				
            /********* 2. Test constraints *******************/
            constraintSatisfiedTF=1;
            if (MEMBRANE)
            {
                printf("Testing Membrane");
                for(i=0;i<N;i++)
                {
                    if (rPropose[i][2] < 0)
                    {
                        constraintSatisfiedTF=0;
                        i=N+1; // shortcut out of the loop
                    }
                } // done checking constraint
            } //finished first constraint
            
            //if (MULTIPLE && constraintSatisfiedTF)
            if (1)
            {
                
                //printf("Testing bound ligands.");
                
                for(ib=0;ib<bSiteTotal;ib++) //for each bound iSite, find the center of the attached ligand
                {
                    switch (ib % 4)
                    {
                            case 0:
                    
                    bSiteCurrent = bSite[ib];
                    bLigandCenter[ib][0] = rPropose[bSiteCurrent][0] + rLigand*e1Propose[bSiteCurrent][0];
                    bLigandCenter[ib][1] = rPropose[bSiteCurrent][1] + rLigand*e1Propose[bSiteCurrent][1];
                    bLigandCenter[ib][2] = rPropose[bSiteCurrent][2] + rLigand*e1Propose[bSiteCurrent][2];
                            break;
                            
                            case 1:
                            
                        bSiteCurrent = bSite[ib];
                        bLigandCenter[ib][0] = rPropose[bSiteCurrent][0] - rLigand*e1Propose[bSiteCurrent][0];
                        bLigandCenter[ib][1] = rPropose[bSiteCurrent][1] - rLigand*e1Propose[bSiteCurrent][1];
                        bLigandCenter[ib][2] = rPropose[bSiteCurrent][2] - rLigand*e1Propose[bSiteCurrent][2];
                            break;
                            
                            case 2:
                        
                        bSiteCurrent = bSite[ib];
                        bLigandCenter[ib][0] = rPropose[bSiteCurrent][0] + rLigand*e2Propose[bSiteCurrent][0];
                        bLigandCenter[ib][1] = rPropose[bSiteCurrent][1] + rLigand*e2Propose[bSiteCurrent][1];
                        bLigandCenter[ib][2] = rPropose[bSiteCurrent][2] + rLigand*e2Propose[bSiteCurrent][2];
                            break;
                            
                            case 3:
                            
                            bSiteCurrent = bSite[ib];
                            bLigandCenter[ib][0] = rPropose[bSiteCurrent][0] - rLigand*e2Propose[bSiteCurrent][0];
                            bLigandCenter[ib][1] = rPropose[bSiteCurrent][1] - rLigand*e2Propose[bSiteCurrent][1];
                            bLigandCenter[ib][2] = rPropose[bSiteCurrent][2] - rLigand*e2Propose[bSiteCurrent][2];
                            break;
                    }
                    
                }
                

                for (ib=0;ib<bSiteTotal;ib++) //for each bound ligand
                {
                    
                    bSiteCurrent = bSite[ib];
                    //printf("This is the current bound site: %ld \n", bSiteCurrent);
                    
                    for(i=0;i<N;i++)// for each joint
                    {
                        if ( (bLigandCenter[ib][0]-rPropose[i][0])*(bLigandCenter[ib][0]-rPropose[i][0]) +
                            (bLigandCenter[ib][1]-rPropose[i][1])*(bLigandCenter[ib][1]-rPropose[i][1]) +
                            (bLigandCenter[ib][2]-rPropose[i][2])*(bLigandCenter[ib][2]-rPropose[i][2]) <= rLigand*rLigand
                            && i != bSite[ib]) //if proposed joint is inside ligand sphere AND joint is not where ligand is attached
                        {
                            constraintSatisfiedTF=0; //constraint not satisfied
                            
                            //printf("Constraint not satisfied");
                            i=N; //shortcut out of loop
                        }
                    }
                    
                    
                    if (1) //if constraint is still satisfied, test ligand sphere with other ligands
                    {
                        for (ib2=(ib+1);ib2<bSiteTotal;ib2++) //for each next ligand
                        {
                            
                            if ((bLigandCenter[ib][0]-bLigandCenter[ib2][0])*(bLigandCenter[ib][0]-bLigandCenter[ib2][0])+(bLigandCenter[ib][1]-bLigandCenter[ib2][1])*(bLigandCenter[ib][1]-bLigandCenter[ib2][1])+(bLigandCenter[ib][2]-bLigandCenter[ib2][2])*(bLigandCenter[ib][2]-bLigandCenter[ib2][2])<= (2*rLigand)*(2*rLigand)) //if distance between centers is less than 2*rLigand, then ligands are intersecting
                                
                            {
                                constraintSatisfiedTF=0; //constraint not satisfied
                                ib2=bSiteTotal; //shortcut out of loop
                            }
                        }
                        
                    }
                    
                }

            } //finished second constraint
            
            constraintProposalsTotal++;
        } //finish constraint while loop
        
        if (constraintProposalsTotal >= CPMAX)
        {
            printf("Exceeded maximum proposals.\n");
            fflush(stdout);
            
            exit(0);
            
        }
        /********* 3. Metropolis test *******************/
        // We now have a propsoal configuration that passes the constraints.
        // Step 3 is to see if it passes our acceptance test (Metropolis test).
        
        // Compute energy
        ENew = -rPropose[N-1][2]*Force; // Energy in units of kBT. Force in units of kBT/Kuhn
		
        if (  TWISTER < exp(E-ENew) ) //always accepts if ENew<E, accepts with normal (?) probability if ENew>E
		{
            E = ENew;
            
            // Make configuration into the proposal configuration
            for(i=iPropose;i<N;i++)
			{
				phi[i]   = phiPropose[i];
				theta[i] = thetaPropose[i];
				psi[i]   = psiPropose[i];
				
				r[i][0] = rPropose[i][0];
				r[i][1] = rPropose[i][1];
				r[i][2] = rPropose[i][2];
                
                t[i][0] = tPropose[i][0];
                t[i][1] = tPropose[i][1];
                t[i][2] = tPropose[i][2];
                
                e1[i][0] = e1Propose[i][0];
                e1[i][1] = e1Propose[i][1];
                e1[i][2] = e1Propose[i][2];
                
                e2[i][0] = e2Propose[i][0];
                e2[i][1] = e2Propose[i][1];
                e2[i][2] = e2Propose[i][2];
             
			}
			if(iPropose==0)
				accepts[0] ++;
			else 		
				accepts[1] ++;
			
		}
        
        /********* 4. Data collection and output to file *******************/
        // check if blocking sphere
        if (1)
        {
            for(iy=0;iy<iSiteTotal;iy++)
            {
                iSiteCurrent = iSite[iy];
                iLigandCenter[iy][0] = r[iSiteCurrent][0] + rLigand*e1[iSiteCurrent][0];
                iLigandCenter[iy][1] = r[iSiteCurrent][1] + rLigand*e1[iSiteCurrent][1];
                iLigandCenter[iy][2] = r[iSiteCurrent][2] + rLigand*e1[iSiteCurrent][2];
            
                stericOcclusion[iy] = 0;
            }
//            
//            baseLigandCenter[0] = rBase[0] + rLigand*e1Base[0];
//            baseLigandCenter[1] = rBase[1] + rLigand*e1Base[1];
//            baseLigandCenter[2] = rBase[2] + rLigand*e1Base[2];
//            
//            stericOcclusionBase = 0;
            
            
            for(iy=0; iy<iSiteTotal;iy++)
            {
                for (ib=0;ib<bSiteTotal;ib++)
                {
                    if(iSite[iy]==bSite[ib]) //test if iSite is bound already
                    {
                        stericOcclusion[iy]++;
                        ib=bSiteTotal;
                    }
                }//didn't include base - assuming can't be bound to base
                
                
                if (stericOcclusion[iy]==0) //if not occluded yet, test membrane or base
                {
                if (MEMBRANE)
                {
                    // check if sphere violates membrane
                    if (iLigandCenter[iy][2]<rLigand)
                        stericOcclusion[iy]++;

                }
                else
                {
                    // check if sphere violates base
                    if ( (iLigandCenter[iy][0])*(iLigandCenter[iy][0]) +
                        (iLigandCenter[iy][1])*(iLigandCenter[iy][1]) +
                        (iLigandCenter[iy][2])*(iLigandCenter[iy][2]) <= rLigand*rLigand )
                        stericOcclusion[iy]++;
                    
                    //didn't include base - don't want the base to violate the base
                }
                }
                
                if(stericOcclusion[iy]==0) //if not occluded yet, test joints
                {
                for(i=0;i<N;i++)
                {
                    if ( (iLigandCenter[iy][0]-r[i][0])*(iLigandCenter[iy][0]-r[i][0]) +
                         (iLigandCenter[iy][1]-r[i][1])*(iLigandCenter[iy][1]-r[i][1]) +
                         (iLigandCenter[iy][2]-r[i][2])*(iLigandCenter[iy][2]-r[i][2]) <= rLigand*rLigand
                        && i != iSite[iy])
                    {
                        stericOcclusion[iy]++;
                        i=N; // shortcut out of the loop
                    }
                }
                }
                if (MULTIPLE && (stericOcclusion[iy]==0)) //if there are multiple ligands and not occluded yet, test other ligands
                {
                    for(ib=0;ib<bSiteTotal;ib++) //for each bound iSite, find the center of the attached ligand
                    {
                        bSiteCurrent = bSite[ib];

                        bLigandCenter[ib][0] = r[bSiteCurrent][0] + rLigand*e1[bSiteCurrent][0];
                        bLigandCenter[ib][1] = r[bSiteCurrent][1] + rLigand*e1[bSiteCurrent][1];
                        bLigandCenter[ib][2] = r[bSiteCurrent][2] + rLigand*e1[bSiteCurrent][2];
                        
                        //would it be better to reinitialize bLigandCenter before writing over it?
                    
                        if ((iLigandCenter[iy][0]-bLigandCenter[ib][0])*(iLigandCenter[iy][0]-bLigandCenter[ib][0])+(iLigandCenter[iy][1]-bLigandCenter[ib][1])*(iLigandCenter[iy][1]-bLigandCenter[ib][1])+(iLigandCenter[iy][2]-bLigandCenter[ib][2])*(iLigandCenter[iy][2]-bLigandCenter[ib][2])<=(2*rLigand)*(2*rLigand))
                        // if potential ligand intersects with bound ligand
                        
                        {
                            stericOcclusion[iy]++;
                            ib=bSiteTotal; //shortcut out of the loop
                        }
                    }
                }
            }
        
        

            if (!MEMBRANE) //check occlusion at base if there is no membrane
            {
                //initialize ligand center and stericOcclusion
                baseLigandCenter[0] = rBase[0] + rLigand*e1Base[0];
                baseLigandCenter[1] = rBase[1] + rLigand*e1Base[1];
                baseLigandCenter[2] = rBase[2] + rLigand*e1Base[2];
                
                stericOcclusionBase = 0;

                
                //check occlusion with joints
                for(i=0;i<N;i++)
                {
                    if ( (baseLigandCenter[0]-r[i][0])*(baseLigandCenter[0]-r[i][0]) +
                        (baseLigandCenter[1]-r[i][1])*(baseLigandCenter[1]-r[i][1]) +
                        (baseLigandCenter[2]-r[i][2])*(baseLigandCenter[2]-r[i][2]) <= rLigand*rLigand)
                    {
                        stericOcclusionBase+=N;
                        i=N; // shortcut out of the loop
                    }
                }
                
                
                
                //check occlusion with other ligands  - to test if they "deliver" their cargo
                //or do we want to test Occlusion of base with ligands? Could do both?
                if (MULTIPLE && stericOcclusionBase==0)
                {
                    //initialize
                    for (ib=0;ib<bSiteTotal;ib++)
                    {
                        bSiteCurrent = bSite[ib];
                        
                        bLigandCenter[ib][0] = r[bSiteCurrent][0] + rLigand*e1[bSiteCurrent][0];
                        bLigandCenter[ib][1] = r[bSiteCurrent][1] + rLigand*e1[bSiteCurrent][1];
                        bLigandCenter[ib][2] = r[bSiteCurrent][2] + rLigand*e1[bSiteCurrent][2];
                        
                        boundToBaseDeliver[ib]=0;
                    }
                    
                    
//                    switch (deliveryMethod)
//                    {
//                            case 0:
                                //for each bound iSite, test if bound ligand intersects with base ligand site
                                for (ib=0; ib<bSiteTotal;ib++)
                                {
                                    if ((baseLigandCenter[0]-bLigandCenter[ib][0])*(baseLigandCenter[0]-bLigandCenter[ib][0])+(baseLigandCenter[1]-bLigandCenter[ib][1])*(baseLigandCenter[1]-bLigandCenter[ib][1])+(baseLigandCenter[2]-bLigandCenter[ib][2])*(baseLigandCenter[2]-bLigandCenter[ib][2]) <= (2*rLigand)*(2*rLigand))
                                    {
                                        boundToBaseDeliver[ib]++;
                                        stericOcclusionBase++;
                                    }
                                }
                            //break;
                            
//                            case 1:
//
//                                //for each bound iSite, test if bound ligand is within "delivery" distance
//                                for (ib=0;ib<bSiteTotal;ib++)
//                                {
//                                    if ((bLigandCenter[ib][0])*(bLigandCenter[ib][0])+(bLigandCenter[ib][1])*(bLigandCenter[ib][1])+(bLigandCenter[ib][2])*(bLigandCenter[ib][2]) <= deliveryDistance)
//                                    {
//                                        boundToBaseDeliver[ib]++;
//                                    }
//                                    
//                                    //test if bound ligand intersects base site
//                                    if ((baseLigandCenter[0]-bLigandCenter[ib][0])*(baseLigandCenter[0]-bLigandCenter[ib][0])+(baseLigandCenter[1]-bLigandCenter[ib][1])*(baseLigandCenter[1]-bLigandCenter[ib][1])+(baseLigandCenter[2]-bLigandCenter[ib][2])*(baseLigandCenter[2]-bLigandCenter[ib][2]) <= (2*rLigand)*(2*rLigand))
//                                    {
//                                        stericOcclusionBase++;
//                                    }
//                                }
//                            break;
//                     }
                }
                
            
            }
        
        } //end checking occlusion

		if (1)
		{
            /********* 5. Test for stationarity *******************/
            // This test for stationarity automatically stops the iterations when the distribution is converged to within the tolerance KSCRITICAL.
            if (nt == 2*NTCHECK)
				appendBins();
			
			if (nt == ntNextStationarityCheck)
				stationarity();
		}
		
        // output to time series file
		dataRecording();
        
        /********* 6. Increment time *******************/
		nt++;
		
	} // done time loop
    
    // finalize summary statistics
    finalizeSummary();

}



/********************************************************************************************************/
// Test for convergence ("stationarity")
void stationarity()
{
	double cdf1, cdf2;
	ksStatistic = 0;
	cdf1=0; 
	cdf2=0;
	for(iBin=0;iBin<NBINS;iBin++)
	{
		
		cdf1 += (double)rMCountsPrevious[iBin]/(((double)nt-(double)NTCHECK)/2.0);
		cdf2 += (double)rMCounts[iBin]/(((double)nt-(double)NTCHECK)/2.0);
		
		if (fabs(cdf1-cdf2)>ksStatistic)
			ksStatistic = fabs(cdf1-cdf2);
		//printf("rMCounts[%d]: %d, %d \t\t\t cdf: %f, %f\n", iBin, rMCountsPrevious[iBin],rMCounts[iBin], cdf1, cdf2);
		
	}
	
	//printf("ksStatistic=%f\n", ksStatistic);
	
	if (ksStatistic<KSCRITICAL)
		convergedTF=1;
	
	appendBins();
	
	ntNextStationarityCheck = 2*(ntNextStationarityCheck-NTCHECK)+NTCHECK;
	return;
}

// Helper function needed by test for convergence ("stationarity")
void appendBins()
{
    for(iBin=0;iBin<NBINS;iBin++)
    {
        rMCountsPrevious[iBin] += rMCounts[iBin];
        rMCounts[iBin] = 0;
    }
}

/********************************************************************************************************/
// This function rotates the orthogonal vectors tIn, e1In and e2In by angles (phiHere, thetaHere, psiHere) in their own frame of reference
void rotate(double *tIn, double *e1In, double *e2In, double *tOut, double *e1Out, double *e2Out, double phiHere, double thetaHere,double psiHere)
{	
	// R Local
    RLocal[0][0] =   cos(thetaHere)*cos(psiHere);
    RLocal[0][1] =   cos(phiHere)*sin(psiHere)     + sin(phiHere)*sin(thetaHere)*cos(psiHere);
    RLocal[0][2] =   sin(phiHere)*sin(psiHere)     - cos(phiHere)*sin(thetaHere)*cos(psiHere);
    RLocal[1][0] =  -cos(thetaHere)*sin(psiHere);
    RLocal[1][1] =   cos(phiHere)*cos(psiHere)     - sin(phiHere)*sin(thetaHere)*sin(psiHere);
    RLocal[1][2] =   sin(phiHere)*cos(psiHere)     + cos(phiHere)*sin(thetaHere)*sin(psiHere);
    RLocal[2][0] =   sin(thetaHere);
    RLocal[2][1] =  -sin(phiHere)*cos(thetaHere);
    RLocal[2][2] =   cos(phiHere)*cos(thetaHere);
    
    
	// R Global
	RGlobal[0][0] = (*(e1In+0))* ((*(e1In+0))* RLocal[0][0] + (*(e2In+0))* RLocal[1][0] + RLocal[2][0]* (*(tIn+0)))
                  + (*(e2In+0))* ((*(e1In+0))* RLocal[0][1] + (*(e2In+0))* RLocal[1][1] + RLocal[2][1]* (*(tIn+0)))
                  + (*(tIn+0))* ((*(e1In+0))* RLocal[0][2] + (*(e2In+0))* RLocal[1][2] + RLocal[2][2]* (*(tIn+0)));
	RGlobal[0][1] = (*(e1In+1))* ((*(e1In+0))* RLocal[0][0] + (*(e2In+0))* RLocal[1][0] + RLocal[2][0]* (*(tIn+0)))
                  + (*(e2In+1))* ((*(e1In+0))* RLocal[0][1] + (*(e2In+0))* RLocal[1][1] + RLocal[2][1]* (*(tIn+0)))
                  + (*(tIn+1))* ((*(e1In+0))* RLocal[0][2] + (*(e2In+0))* RLocal[1][2] + RLocal[2][2]* (*(tIn+0)));
	RGlobal[0][2] = (*(e1In+2))* ((*(e1In+0))* RLocal[0][0] + (*(e2In+0))* RLocal[1][0] + RLocal[2][0]* (*(tIn+0)))
                  + (*(e2In+2))* ((*(e1In+0))* RLocal[0][1] + (*(e2In+0))* RLocal[1][1] + RLocal[2][1]* (*(tIn+0)))
                  + (*(tIn+2))* ((*(e1In+0))* RLocal[0][2] + (*(e2In+0))* RLocal[1][2] + RLocal[2][2]* (*(tIn+0)));
	RGlobal[1][0] = (*(e1In+0))* ((*(e1In+1))* RLocal[0][0] + (*(e2In+1))* RLocal[1][0] + RLocal[2][0]* (*(tIn+1)))
                  + (*(e2In+0))* ((*(e1In+1))* RLocal[0][1] + (*(e2In+1))* RLocal[1][1] + RLocal[2][1]* (*(tIn+1)))
                  + (*(tIn+0))* ((*(e1In+1))* RLocal[0][2] + (*(e2In+1))* RLocal[1][2] + RLocal[2][2]* (*(tIn+1)));
	RGlobal[1][1] = (*(e1In+1))* ((*(e1In+1))* RLocal[0][0] + (*(e2In+1))* RLocal[1][0] + RLocal[2][0]* (*(tIn+1)))
                  + (*(e2In+1))* ((*(e1In+1))* RLocal[0][1] + (*(e2In+1))* RLocal[1][1] + RLocal[2][1]* (*(tIn+1)))
                  + (*(tIn+1))* ((*(e1In+1))* RLocal[0][2] + (*(e2In+1))* RLocal[1][2] + RLocal[2][2]* (*(tIn+1)));
	RGlobal[1][2] = (*(e1In+2))* ((*(e1In+1))* RLocal[0][0] + (*(e2In+1))* RLocal[1][0] + RLocal[2][0]* (*(tIn+1)))
                  + (*(e2In+2))* ((*(e1In+1))* RLocal[0][1] + (*(e2In+1))* RLocal[1][1] + RLocal[2][1]* (*(tIn+1)))
                  + (*(tIn+2))* ((*(e1In+1))* RLocal[0][2] + (*(e2In+1))* RLocal[1][2] + RLocal[2][2]* (*(tIn+1)));
	RGlobal[2][0] = (*(e1In+0))* ((*(e1In+2))* RLocal[0][0] + (*(e2In+2))* RLocal[1][0] + RLocal[2][0]* (*(tIn+2)))
                  + (*(e2In+0))* ((*(e1In+2))* RLocal[0][1] + (*(e2In+2))* RLocal[1][1] + RLocal[2][1]* (*(tIn+2)))
                  + (*(tIn+0))* ((*(e1In+2))* RLocal[0][2] + (*(e2In+2))* RLocal[1][2] + RLocal[2][2]* (*(tIn+2)));
	RGlobal[2][1] = (*(e1In+1))* ((*(e1In+2))* RLocal[0][0] + (*(e2In+2))* RLocal[1][0] + RLocal[2][0]* (*(tIn+2)))
                  + (*(e2In+1))* ((*(e1In+2))* RLocal[0][1] + (*(e2In+2))* RLocal[1][1] + RLocal[2][1]* (*(tIn+2)))
                  + (*(tIn+1))* ((*(e1In+2))* RLocal[0][2] + (*(e2In+2))* RLocal[1][2] + RLocal[2][2]* (*(tIn+2)));
	RGlobal[2][2] = (*(e1In+2))* ((*(e1In+2))* RLocal[0][0] + (*(e2In+2))* RLocal[1][0] + RLocal[2][0]* (*(tIn+2)))
                  + (*(e2In+2))* ((*(e1In+2))* RLocal[0][1] + (*(e2In+2))* RLocal[1][1] + RLocal[2][1]* (*(tIn+2)))
                  + (*(tIn+2))* ((*(e1In+2))* RLocal[0][2] + (*(e2In+2))* RLocal[1][2] + RLocal[2][2]* (*(tIn+2)));
	
	// New unit vectors
	*(tOut+0) = RGlobal[0][0]*(*(tIn+0)) + RGlobal[0][1]*(*(tIn+1)) + RGlobal[0][2]*(*(tIn+2));
	*(tOut+1) = RGlobal[1][0]*(*(tIn+0)) + RGlobal[1][1]*(*(tIn+1)) + RGlobal[1][2]*(*(tIn+2));
	*(tOut+2) = RGlobal[2][0]*(*(tIn+0)) + RGlobal[2][1]*(*(tIn+1)) + RGlobal[2][2]*(*(tIn+2));	
	
	*(e1Out+0) = RGlobal[0][0]*(*(e1In+0)) + RGlobal[0][1]*(*(e1In+1)) + RGlobal[0][2]*(*(e1In+2));
	*(e1Out+1) = RGlobal[1][0]*(*(e1In+0)) + RGlobal[1][1]*(*(e1In+1)) + RGlobal[1][2]*(*(e1In+2));
	*(e1Out+2) = RGlobal[2][0]*(*(e1In+0)) + RGlobal[2][1]*(*(e1In+1)) + RGlobal[2][2]*(*(e1In+2));	
	
	*(e2Out+0) = RGlobal[0][0]*(*(e2In+0)) + RGlobal[0][1]*(*(e2In+1)) + RGlobal[0][2]*(*(e2In+2));
	*(e2Out+1) = RGlobal[1][0]*(*(e2In+0)) + RGlobal[1][1]*(*(e2In+1)) + RGlobal[1][2]*(*(e2In+2));
	*(e2Out+2) = RGlobal[2][0]*(*(e2In+0)) + RGlobal[2][1]*(*(e2In+1)) + RGlobal[2][2]*(*(e2In+2));
    

    // Re-normalize and re-orthogonalize unit vectors (using Modified Gram Schmidt algorithm)
    norm = sqrt((*(tOut+0))*(*(tOut+0)) + (*(tOut+1))*(*(tOut+1)) + (*(tOut+2))*(*(tOut+2)));
    *(tOut+0) = *(tOut+0)/norm;
    *(tOut+1) = *(tOut+1)/norm;
    *(tOut+2) = *(tOut+2)/norm;
        
    e1_dot_t = (*(tOut+0))*(*(e1Out+0))+(*(tOut+1))*(*(e1Out+1))+(*(tOut+2))*(*(e1Out+2));
    
    *(e1Out+0) = *(e1Out+0) - e1_dot_t*(*(tOut+0));
    *(e1Out+1) = *(e1Out+1) - e1_dot_t*(*(tOut+1));
    *(e1Out+2) = *(e1Out+2) - e1_dot_t*(*(tOut+2));
    
    norm = sqrt((*(e1Out+0))*(*(e1Out+0)) + (*(e1Out+1))*(*(e1Out+1)) + (*(e1Out+2))*(*(e1Out+2)));
    *(e1Out+0) = *(e1Out+0)/norm;
    *(e1Out+1) = *(e1Out+1)/norm;
    *(e1Out+2) = *(e1Out+2)/norm;

    e2_dot_t  = (*(tOut+0))*(*(e2Out+0))+(*(tOut+1))*(*(e2Out+1))+(*(tOut+2))*(*(e2Out+2));
    e2_dot_e1 =(*(e1Out+0))*(*(e2Out+0))+(*(e1Out+1))*(*(e2Out+1))+(*(e1Out+2))*(*(e2Out+2));
    
    *(e2Out+0) = *(e2Out+0) - e2_dot_t*(*(tOut+0)) - e2_dot_e1*(*(e1Out+0));
    *(e2Out+1) = *(e2Out+1) - e2_dot_t*(*(tOut+1)) - e2_dot_e1*(*(e1Out+1));
    *(e2Out+2) = *(e2Out+2) - e2_dot_t*(*(tOut+2)) - e2_dot_e1*(*(e1Out+2));
    
    norm = sqrt((*(e2Out+0))*(*(e2Out+0)) + (*(e2Out+1))*(*(e2Out+1)) + (*(e2Out+2))*(*(e2Out+2)));
    *(e2Out+0) = *(e2Out+0)/norm;
    *(e2Out+1) = *(e2Out+1)/norm;
    *(e2Out+2) = *(e2Out+2)/norm;
   
    if (0) // print stuff out for debugging
    {
        // PRINT RGlobal
        printf("At nt=%ld, i=%ld\n", nt, i);
        printf(" Rotation matrix:\n");
        printf(" %f,", RGlobal[0][0]);
        printf(" %f,", RGlobal[0][1]);
        printf(" %f,", RGlobal[0][2]);
        printf(" \n");
        printf(" %f,", RGlobal[1][0]);
        printf(" %f,", RGlobal[1][1]);
        printf(" %f,", RGlobal[1][1]);
        printf(" \n");
        printf(" %f,", RGlobal[2][0]);
        printf(" %f,", RGlobal[2][1]);
        printf(" %f,", RGlobal[2][2]);
        printf(" \n");
        
        // PRINT RTR (should be the identity matrix)
        printf(" RTR:\n");
        
        printf(" %f,", RGlobal[0][0]*RGlobal[0][0]+RGlobal[1][0]*RGlobal[1][0]+RGlobal[2][0]*RGlobal[2][0]);
        printf(" %f,", RGlobal[0][0]*RGlobal[0][1]+RGlobal[1][0]*RGlobal[1][1]+RGlobal[2][0]*RGlobal[2][1]);
        printf(" %f,", RGlobal[0][0]*RGlobal[0][2]+RGlobal[1][0]*RGlobal[1][2]+RGlobal[2][0]*RGlobal[2][2]);
        printf(" \n");
        printf(" %f,", RGlobal[0][1]*RGlobal[0][0]+RGlobal[1][1]*RGlobal[1][0]+RGlobal[2][1]*RGlobal[2][0]);
        printf(" %f,", RGlobal[0][1]*RGlobal[0][1]+RGlobal[1][1]*RGlobal[1][1]+RGlobal[2][1]*RGlobal[2][1]);
        printf(" %f,", RGlobal[0][1]*RGlobal[0][2]+RGlobal[1][1]*RGlobal[1][2]+RGlobal[2][1]*RGlobal[2][2]);
        printf(" \n");
        printf(" %f,", RGlobal[0][2]*RGlobal[0][0]+RGlobal[1][2]*RGlobal[1][0]+RGlobal[2][2]*RGlobal[2][0]);
        printf(" %f,", RGlobal[0][2]*RGlobal[0][1]+RGlobal[1][2]*RGlobal[1][1]+RGlobal[2][2]*RGlobal[2][1]);
        printf(" %f,", RGlobal[0][2]*RGlobal[0][2]+RGlobal[1][2]*RGlobal[1][2]+RGlobal[2][2]*RGlobal[2][2]);
        printf(" \n");
        
        // PRINT dot products of Cosserat vectors (should have zero dot products)
        printf("t dot e1 = %lf, \ne1 dot e2 = %lf\n",
               (*(tOut+0))*(*(e1Out+0))+(*(tOut+1))*(*(e1Out+1))+(*(tOut+2))*(*(e1Out+2)),
               (*(e2Out+0))*(*(e1Out+0))+(*(e2Out+1))*(*(e1Out+1))+(*(e2Out+2))*(*(e1Out+2)));
        
        printf(" \n");
        printf(" \n");

    }
    
	
	return;
	
}

