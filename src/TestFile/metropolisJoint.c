/*** Allard Group jun.allard@uci.edu                    ***/

void metropolisJoint();
void stationarity();
void appendBins();
void rotate(double *tIn, double *e1In, double *e2In, double *tOut, double *e1Out, double *e2Out, double phiHere, double thetaHere,double psiHere);

void metropolisJoint()
{
    
    /********* INITIALIZE ISITES *******************/
//    
//    
//    for (iy=0;iy<iSiteTot;iy++)
//    {
//        iSite[iy]=0;
//    }
    
    //switch (commandiSites)
    //{
            //case 0:
    
                switch (testRun)
                {
                    case 0:  // iSites initialized for human CD3Zeta-Chain
                        
                        iSiteTot = 7;
                        
                        
                        for (iy=0;iy<iSiteTot;iy++)
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
                        
                        iSiteTot = 3;
                        
                        for (iy=0;iy<iSiteTot;iy++)
                        {
                            iSite[iy]=0;
                        }
                        
                        iSite[0]=0;
                        iSite[1]=3;
                        iSite[2]=7;
                        break;
                        
                        
                    case 2: //test case 2 - stiffen none, but test all iSites, make Ratio half of Ratio for case 1
                        
                        iSiteTot = 7;
                        
                        for (iy=0;iy<iSiteTot;iy++)
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
                        
                        iSiteTot = 2;
                        
                        for(iy=0;iy<iSiteTot;iy++)
                        {
                            iSite[iy]=0;
                        }
                        
                        iSite[0]=2;
                        iSite[1]=4;
                }
            
            //break;
            
            
            
//            case 1:
//            
//                for (iy=0;iy<iSiteTot;iy++)
//                {
//                    iSite[iy]=0;
//                }
//            
//            printf("Total iSites: %ld", iSiteTot);
//            
//            //for debugging
//            for(iy=0;iy<iSiteTot;iy++)
//        {
//            printf("iSite[%ld] = %ld", iy, iSite[iy]);
//        }
//            
//                //char input[] = iSiteLocations;
//                printf("I want to split this into tokens: %s", input);
//                char* strArray[NMAX];
//                char *token = strtok(input, " ");
//                
//                //for(int j = 0; j<NMAX;j++)
//                //{
//                //    strArray[j] = new char[4];
//                //}
//                
//                while(token != NULL)
//                {
//                    strcpy(strArray[st],token);
//                    printf("This is the next token: %s\n",token); //for debugging
//                    token = strtok(NULL, " ");
//                    st++;
//                }
//            
//                //for debugging
//            
//                if (iSiteTot!=st)
//                {
//                    printf("Warning! iSite Total is %ld but Number of iSites in String is %ld ", iSiteTot, st);
//                }
//            
//                //reassign strings as doubles
//                for(iy=0;iy<st;iy++)
//                {
//                    iSite[iy]=atof(strArray[iy]);
//                }
//            
//                //for debugging
//                for(iy=0;iy<iSiteTot;iy++)
//                {
//                    printf("iSite[%ld] = %ld", iy, iSite[iy]);
//                }
//            
//            break;
//    }

    /********* INITIALIZE BOUND ISITES *******************/
    
    if (MULTIPLE)
    {
    
    //switch () //add more cases later
    //{
        //case 0: // arbitrary subset are occupied
            
            boundTotal = 1; //total number of iSites bound
            
            iSiteBound[0]=2; //currently identifying by location
            //what is the best way to do this - identify by location or identify by iSite number?
            //pro for location - can specify locations other than iSites to be bound - but then might want to change name
    //}
    }
    
    /********* STIFFEN SEGMENTS *******************/
    
    
    if (StiffenRange > -1) //stiffen only if StiffenRange is 0 or greater
    {
    
    
        //initializes phosiSites to 0 (none phosphorylated)
        for(ty=0;ty<iSiteTot;ty++)
        {
            phosiSites[ty]=0;
        }
    

    
        //initializes phosiSites to 0 (none phosphorylated)
        for(ty=0;ty<iSiteTot;ty++)
        {
            phosiSites[ty]=0;
        }
    
    
        printf("This is a string: %s\n", phosphorylatediSites);
    
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
            
                sscanf(phosphorylatediSites,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &phosiSites[0],&phosiSites[1],&phosiSites[2],&phosiSites[3], &phosiSites[4],&phosiSites[5],&phosiSites[6],&phosiSites[7],&phosiSites[8],&phosiSites[9],&phosiSites[10],&phosiSites[11],&phosiSites[12],&phosiSites[13]);
                break;
        }
    
    
        for (iy=0;iy<iSiteTot;iy++)
        {
            printf("phosiSites[ %ld ] =  %f\n",iy, phosiSites[iy]);
        }
    
        //initializes stiffened rods to 0 (none stiff)
        for(i=0;i<N;i++)
        {
            Stiff[i] =0;
        }
    

    
        //Stiffen segments
        for(ty=0;ty<iSiteTot;ty++)
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
		r[i][2]  = 0;
		
		t[i][0]  = 0;
		t[i][1]  = 0;
		t[i][2]  = 1;
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
        printf("This is the joint we are rotating: %ld", iPropose);
        
        
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
		
        constraintProposals = 0;
        constraintProposalsTotal = 0;
        long rounds = 0;
		constraintSatisfiedTF = 0;
		while(!constraintSatisfiedTF) // keep generating proposal configurations until we get one that satisfies the constraint
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
            rotate(&tBase[0], &e1Base[0], &e2Base[0], &tPropose[0][0], &e1Propose[0][0], &e2Propose[0][0],
                    phiPropose[0], thetaPropose[0], psiPropose[0]);
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
            
            
            
            printf("This is rPropose: \n %f\n 2: %f\n 3:%f\n", rPropose[42][0],rPropose[42][1],rPropose[42][2]);
            printf("This is e1Propose: \n %f\n 2: %f\n 3:%f\n", e1Propose[42][0],e1Propose[42][1],e1Propose[42][2]);
				
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
            if (0)
            {
                
                //printf("Testing bound ligands.");
                
                for(ib=0;ib<boundTotal;ib++) //for each bound iSite, find the center of the attached ligand
                {
                    currentBoundSite = iSiteBound[ib];
                    rLigandCenterBound[ib][0] = rPropose[currentBoundSite][0] + rLigand*e1Propose[currentBoundSite][0];
                    rLigandCenterBound[ib][1] = rPropose[currentBoundSite][1] + rLigand*e1Propose[currentBoundSite][1];
                    rLigandCenterBound[ib][2] = rPropose[currentBoundSite][2] + rLigand*e1Propose[currentBoundSite][2];

                }
                
                
               // printf("This is rPropose: \n %f\n 2: %f\n 3:%f\n", rPropose[42][0],rPropose[42][1],rPropose[42][2]);
               // printf("This is e1Propose: \n %f\n 2: %f\n 3:%f\n", e1Propose[42][0],e1Propose[42][1],e1Propose[42][2]);

                printf("This is the center:\n 1: %f\n 2: %f\n 3:%f\n", rLigandCenterBound[0][0],rLigandCenterBound[0][1],rLigandCenterBound[0][2]);
                for (ib=0;ib<boundTotal;ib++) //for each bound ligand
                {
                    
                    currentBoundSite = iSiteBound[ib];
                    //printf("This is the current bound site: %ld", currentBoundSite);
                    
                    for(i=0;i<N;i++)// for each joint
                    {
                        if ( (rLigandCenterBound[ib][0]-rPropose[i][0])*(rLigandCenterBound[ib][0]-rPropose[i][0]) +
                            (rLigandCenterBound[ib][1]-rPropose[i][1])*(rLigandCenterBound[ib][1]-rPropose[i][1]) +
                            (rLigandCenterBound[ib][2]-rPropose[i][2])*(rLigandCenterBound[ib][2]-rPropose[i][2]) <= rLigand*rLigand
                            && i != iSiteBound[ib]) //if proposed joint is inside ligand sphere AND joint is not where ligand is attached
                        {
                            constraintSatisfiedTF=0; //constraint not satisfied
                            
                            //printf("Constraint not satisfied");
                            i=N; //shortcut out of loop
                        }
                    }
                    
                    
                    if (0) //if constraint is still satisfied, test ligand sphere with other ligands
                    {
                        //printf("Testing everything");
                        for (ib2=(ib+1);ib2<boundTotal;ib2++) //for each next ligand
                        {
                            printf("Ooops!  There shouldn't be anything for me to test!");
                            
                            if ((rLigandCenterBound[ib][0]-rLigandCenterBound[ib2][0])*(rLigandCenterBound[ib][0]-rLigandCenterBound[ib2][0])+(rLigandCenterBound[ib][1]-rLigandCenterBound[ib2][1])*(rLigandCenterBound[ib][1]-rLigandCenterBound[ib2][1])+(rLigandCenterBound[ib][2]-rLigandCenterBound[ib2][2])*(rLigandCenterBound[ib][2]-rLigandCenterBound[ib2][2])<= (2*rLigand)*(2*rLigand)) //if distance between centers is less than 2*rLigand, then ligands are intersecting
                                
                            {
                                constraintSatisfiedTF=0; //constraint not satisfied
                                ib2=boundTotal; //shortcut out of loop
                            }
                        }
                        
                    }
                    
                }

            } //finished second constraint
            
//            constraintProposals++;
            constraintProposalsTotal++;
//            
//            if ( ( (double) constraintProposals/100000 ) > 3)
//            {
//                printf("WHILE LOOP: %ld ", rounds);
//                rounds++;
//                constraintProposals=0;
//            }
            
        } // finished while-loop to impose constraint
        
        printf("YAY!  I'm done with the WHILE LOOP!");
        
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
            for(iy=0;iy<iSiteTot;iy++)
            {
                iSiteCurrent = iSite[iy];
                rLigandCenter[iy][0] = r[iSiteCurrent][0] + rLigand*e1[iSiteCurrent][0];
                rLigandCenter[iy][1] = r[iSiteCurrent][1] + rLigand*e1[iSiteCurrent][1];
                rLigandCenter[iy][2] = r[iSiteCurrent][2] + rLigand*e1[iSiteCurrent][2];
            
                stericOcclusion[iy] = 0;
            }
//            
//            rLigandCenterBase[0] = rBase[0] + rLigand*e1Base[0];
//            rLigandCenterBase[1] = rBase[1] + rLigand*e1Base[1];
//            rLigandCenterBase[2] = rBase[2] + rLigand*e1Base[2];
//            
//            stericOcclusionBase = 0;
            
            
            for(iy=0; iy<iSiteTot;iy++)
            {
                for (ib=0;ib<boundTotal;ib++)
                {
                    if(iSite[iy]==iSiteBound[ib]) //test if iSite is bound already
                    {
                        stericOcclusion[iy]++;
                        ib=boundTotal;
                    }
                }//didn't include base - assuming can't be bound to base
                
                
                if (stericOcclusion[iy]==0) //if not occluded yet, test membrane or base
                {
                if (MEMBRANE)
                {
                    // check if sphere violates membrane
                    if (rLigandCenter[iy][2]<rLigand)
                        stericOcclusion[iy]++;

                }
                else
                {
                    // check if sphere violates base
                    if ( (rLigandCenter[iy][0])*(rLigandCenter[iy][0]) +
                        (rLigandCenter[iy][1])*(rLigandCenter[iy][1]) +
                        (rLigandCenter[iy][2])*(rLigandCenter[iy][2]) <= rLigand*rLigand )
                        stericOcclusion[iy]++;
                    
                    //didn't include base - don't want the base to violate the base
                }
                }
                
                if(stericOcclusion[iy]==0) //if not occluded yet, test joints
                {
                for(i=0;i<N;i++)
                {
                    if ( (rLigandCenter[iy][0]-r[i][0])*(rLigandCenter[iy][0]-r[i][0]) +
                         (rLigandCenter[iy][1]-r[i][1])*(rLigandCenter[iy][1]-r[i][1]) +
                         (rLigandCenter[iy][2]-r[i][2])*(rLigandCenter[iy][2]-r[i][2]) <= rLigand*rLigand
                        && i != iSite[iy])
                    {
                        stericOcclusion[iy]++;
                        i=N; // shortcut out of the loop
                    }
                }
                }
                if (MULTIPLE && (stericOcclusion[iy]==0)) //if there are multiple ligands and not occluded yet, test other ligands
                {
                    for(ib=0;ib<boundTotal;ib++) //for each bound iSite, find the center of the attached ligand
                    {
                        currentBoundSite = iSiteBound[ib];

                        rLigandCenterBound[ib][0] = r[currentBoundSite][0] + rLigand*e1[currentBoundSite][0];
                        rLigandCenterBound[ib][1] = r[currentBoundSite][1] + rLigand*e1[currentBoundSite][1];
                        rLigandCenterBound[ib][2] = r[currentBoundSite][2] + rLigand*e1[currentBoundSite][2];
                        
                        //would it be better to reinitialize rLigandCenterBound before writing over it?
                    
                        if ((rLigandCenter[iy][0]-rLigandCenterBound[ib][0])*(rLigandCenter[iy][0]-rLigandCenterBound[ib][0])+(rLigandCenter[iy][1]-rLigandCenterBound[ib][1])*(rLigandCenter[iy][1]-rLigandCenterBound[ib][1])+(rLigandCenter[iy][2]-rLigandCenterBound[ib][2])*(rLigandCenter[iy][2]-rLigandCenterBound[ib][2])<=(2*rLigand)*(2*rLigand))
                        // if potential ligand intersects with bound ligand
                        
                        {
                            stericOcclusion[iy]++;
                            ib=boundTotal; //shortcut out of the loop
                        }
                    }
                }
            }
        

            if (!MEMBRANE) //check occlusion at base if there is no membrane
            {
                //initialize ligand center and stericOcclusion
                rLigandCenterBase[0] = rBase[0] + rLigand*e1Base[0];
                rLigandCenterBase[1] = rBase[1] + rLigand*e1Base[1];
                rLigandCenterBase[2] = rBase[2] + rLigand*e1Base[2];
                
                stericOcclusionBase = 0;

                
                //check occlusion with joints
                for(i=0;i<N;i++)
                {
                    if ( (rLigandCenterBase[0]-r[i][0])*(rLigandCenterBase[0]-r[i][0]) +
                        (rLigandCenterBase[1]-r[i][1])*(rLigandCenterBase[1]-r[i][1]) +
                        (rLigandCenterBase[2]-r[i][2])*(rLigandCenterBase[2]-r[i][2]) <= rLigand*rLigand)
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
                    for (ib=0;ib<boundTotal;ib++)
                    {
                        currentBoundSite = iSiteBound[ib];
                        
                        rLigandCenterBound[ib][0] = r[currentBoundSite][0] + rLigand*e1[currentBoundSite][0];
                        rLigandCenterBound[ib][1] = r[currentBoundSite][1] + rLigand*e1[currentBoundSite][1];
                        rLigandCenterBound[ib][2] = r[currentBoundSite][2] + rLigand*e1[currentBoundSite][2];
                        
                        boundToBaseDeliver[ib]=0;
                    }
                    
                    
//                    switch (deliveryMethod)
//                    {
//                            case 0:
                                //for each bound iSite, test if bound ligand intersects with base ligand site
                                for (ib=0; ib<boundTotal;ib++)
                                {
                                    if ((rLigandCenterBase[0]-rLigandCenterBound[ib][0])*(rLigandCenterBase[0]-rLigandCenterBound[ib][0])+(rLigandCenterBase[1]-rLigandCenterBound[ib][1])*(rLigandCenterBase[1]-rLigandCenterBound[ib][1])+(rLigandCenterBase[2]-rLigandCenterBound[ib][2])*(rLigandCenterBase[2]-rLigandCenterBound[ib][2]) <= (2*rLigand)*(2*rLigand))
                                    {
                                        boundToBaseDeliver[ib]++;
                                        stericOcclusionBase++;
                                    }
                                }
                            //break;
                            
//                            case 1:
//
//                                //for each bound iSite, test if bound ligand is within "delivery" distance
//                                for (ib=0;ib<boundTotal;ib++)
//                                {
//                                    if ((rLigandCenterBound[ib][0])*(rLigandCenterBound[ib][0])+(rLigandCenterBound[ib][1])*(rLigandCenterBound[ib][1])+(rLigandCenterBound[ib][2])*(rLigandCenterBound[ib][2]) <= deliveryDistance)
//                                    {
//                                        boundToBaseDeliver[ib]++;
//                                    }
//                                    
//                                    //test if bound ligand intersects base site
//                                    if ((rLigandCenterBase[0]-rLigandCenterBound[ib][0])*(rLigandCenterBase[0]-rLigandCenterBound[ib][0])+(rLigandCenterBase[1]-rLigandCenterBound[ib][1])*(rLigandCenterBase[1]-rLigandCenterBound[ib][1])+(rLigandCenterBase[2]-rLigandCenterBound[ib][2])*(rLigandCenterBase[2]-rLigandCenterBound[ib][2]) <= (2*rLigand)*(2*rLigand))
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

