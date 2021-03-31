#include "MD_run.h"

//define in MD_run.h
//class MD_run


void MD_run::runnvt(int step1){
//------------------------------------------------------------------------------
// DESCRIPTION:
// Molecular dynamics calculation with fixed Volume and slight fluctuated temperature
// Using thermostat : 1, Nose-Hoover Chains; 2, Langevin; 3, Anderson
// Normal Nose-Hoover thermostat method is retained for test.
//------------------------------------------------------------------------------

	TITLE("MD_run","runnvt");
	timer::tick("MD_run","runnvt",'C');
	step=step1+step_rst;
	//the real MD step
	if (!MY_RANK)out<<"step: " << step <<endl;
	
	//change target temperature
	if (step!=1)ReadNewTemp( step );
	
	cout << " ---------------------------------------------------------" << endl;
	
	cout<< " Target temperature  : "<< temperature/K_BOLTZMAN_AU<< " (K)"<<endl;
	
	if(step==1||step%fixTemperature==1)
    {
		for(int k=0;k<numIon;k++)
        {
			if(ionmbl[k].x==0)vel[k].x=0;
			if(ionmbl[k].y==0)vel[k].y=0;
			if(ionmbl[k].z==0)vel[k].z=0;
		}
		scalevel();
	}

	// get the kinetic energy
	twiceKE = GetAtomKE();
	twiceKE = twiceKE * 2;
	
	//Set up forces 
	callforce();

	//print total stress + stress_MD
	if(STRESS)
	{
		cout<<"output Pressure for check!"<<endl;
		double press;
		for(int i=0;i<3;i++)
			press += stress_lcao(i,i)/3;
		press += twiceKE/3/ucell.omega; //output virtual press = 2/3 *Ek/V + sum(sigma[i][i])/3
		double unit_transform = RYDBERG_SI / pow(BOHR_RADIUS_SI,3) * eps8 ;
		cout<<"Virtual Pressure is "<<press*unit_transform<<" Kbar "<<endl;
	}

	for(int k=0;k<numIon;k++){
		if(ionmbl[k].x==0)force[k].x=0;
		if(ionmbl[k].y==0)force[k].y=0;
		if(ionmbl[k].z==0)force[k].z=0;
	}
	maxForce = MAXVALF();
	mdenergy=en.etot/2;
	if(NVT_control == 1) hamiltonian = NHChamiltonian(twiceKE/2, mdenergy);
	else hamiltonian = Conserved(twiceKE/2, mdenergy);
	//----------------------------------------------
	// big loop
	//-----------------------------------------------
	cout<<" "<<std::left<<setw(12)<<"MD_STEP"<<std::left<<setw(12)<< "SystemE"<<std::left<<setw(12)<< "Conserved"<<std::left<<setw(12)<< "DeltaE"<<std::left<<setw(12)<< "Temperature"<<endl;
	
	cout<<" "<<std::left<<setw(12)<<step<<std::left<<setw(12)<< mdenergy<<std::left<<setw(12)<< hamiltonian<<std::left<<setw(12)<< mdenergy-oldEtot<<std::left<<setw(12)<<twiceKE/(3*double(numIon-nfrozen))/K_BOLTZMAN_AU<<endl;
	
	oldEtot=mdenergy;

	if (!MY_RANK)
    {
	      out<<" --------------------------------------------------"<<endl;
	      out<<" Molecular Dynamics (NVT) STEP "<< step<<endl;
	      out<< " --------------------------------------------------"<<endl;
	      ofs_running<<" --------------------------------------------------"<<endl;
	      ofs_running<<" Molecular Dynamics (NVT) STEP "<< step<<endl;
	      ofs_running<< " --------------------------------------------------"<<endl;
	}
	
	
	// Calculate the Mean-Square-Displacement.
	if(step==1&&rstMD==0)
    { 
        //Note: side scheme position before
        //Now turn to middle scheme.
	
	    // New vel obtained (Verlet-like step)
	    for(int  ii=0;ii<numIon;ii++){ 
			mass = allmass[ii];
			vel[ii] = vel[ii] + force[ii]/mass*dt/2.0;
		}
	}
	else 
    {
		for(int  ii=0;ii<numIon;ii++)
        { 
			mass = allmass[ii];
			vel[ii] = vel[ii] + force[ii]/mass*dt/2.0;
		}
		
		twiceKE=GetAtomKE();
		twiceKE = 2 * twiceKE;
		if(NVT_control==1) hamiltonian = NHChamiltonian(twiceKE/2, mdenergy);
		else hamiltonian = Conserved(twiceKE/2, mdenergy);
        //Note: side scheme position before
        //Now turn to middle scheme. 
		for(int  ii=0;ii<numIon;ii++){ 
			mass = allmass[ii];
			vel[ii] = vel[ii] + force[ii]/mass*dt/2.0;
		}
	}

	// Update the Non-Wrapped cartesion coordinates
    if(mdtype==2) scalevel();//choose velocity scaling method

	for(int  ii=0;ii<numIon;ii++)
    { 
        cartNoWrap[ii] = vel[ii]*dt/2.0 +cartNoWrap[ii];
	}
	
	//by zifei
	Integrator(NVT_control);
    twiceKE=GetAtomKE();
    twiceKE = 2 * twiceKE;
	
	for(int  ii=0;ii<numIon;ii++){ 
		cartNoWrap[ii] = vel[ii]*dt/2.0 +cartNoWrap[ii];
	}

	// Calculate the maximal velocities.
	// The step in fractional coordinates
	maxStep = 0;
	Vector3<double> fracStep;
	for(int  ii=0;ii<numIon;ii++)
    { 
		if((pow(vel[ii].x,2)+pow(vel[ii].y,2)+pow(vel[ii].z,2))>maxStep)
        {
		    maxStep = pow(vel[ii].x,2)+pow(vel[ii].y,2)+pow(vel[ii].z,2);
        }
		Mathzone::Cartesian_to_Direct(vel[ii].x*dt/ucell.lat0,vel[ii].y*dt/ucell.lat0,vel[ii].z*dt/ucell.lat0,
					ionlatvec.e11,ionlatvec.e12,ionlatvec.e13,
					ionlatvec.e21,ionlatvec.e22,ionlatvec.e23,
					ionlatvec.e31,ionlatvec.e32,ionlatvec.e33,
					fracStep.x,fracStep.y,fracStep.z);

		taudirac[ii] = taudirac[ii] + fracStep;
	}

	//save the atom position change to DFT module
	moveatoms(step);
	string t("md_pos_");
	t=intTurnTostring(step,t);
	connection2();
	printpos(t,step);
	maxStep = sqrt(maxStep)*dt;

    if (!MY_RANK){
         out<<setw(15)<<"maxForce=     "<<setw(15)<<"maxstep=      "<<setw(15)<<"step=     "<<endl;
         out<<setw(15)<<maxForce<<setw(15)<<maxStep<<setw(15)<<step<<endl;
         out<<step<<" "<< mdenergy<<" "<< hamiltonian<<" "<< mdenergy-oldEtot<<" "<< twiceKE/(3*double(numIon-nfrozen))/K_BOLTZMAN_AU<<endl;
    }
    oldEtot=mdenergy;
    mstout(step);

//2015-09-25, xiaohui
#ifdef __MPI
    if(step%100==0)
    {
        MPI_Bcast(vel,numIon*3,MPI_DOUBLE,0,MPI_COMM_WORLD);
        MPI_Bcast(cartNoWrap,numIon*3,MPI_DOUBLE,0,MPI_COMM_WORLD);
        MPI_Bcast(taudirac,numIon*3,MPI_DOUBLE,0,MPI_COMM_WORLD);
    //    MPI_Bcast(&xLogS,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    //    MPI_Bcast(&vLogS,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    }
#endif
    timer::tick("MD_run","runnvt",'C');
    return;
}

void MD_run::runNVE(int step1){
//-------------------------------------------------------------------------------
// Adiabatic ensemble 
// Molecular dynamics calculation with Verlet algorithm
//-------------------------------------------------------------------------------

    TITLE("MD_run","runNVE");
    timer::tick("MD_run","runNVE",'C');
    step=step1+step_rst;
    if(step==1)
    {
        cout<<" (NVE) Start."<<endl;
    }
    else 
    {
        cout<<" (NVE) step: "<<step<<endl;
    }
  
  cout<<" (NVE) Fronzen Atom Number : "<<nfrozen<<endl;

  // Set up extra output to ion optimizer / MD header
  cout<<"Time interval       : "<< dt*fundamentalTime*1e15<< " (fs)"<<endl;
  cout<<"Target temperature  : "<< temperature/K_BOLTZMAN_AU<< " (K)"<<endl;

    if(step==1){
        for(int k=0;k<numIon;k++)
        {
            if(ionmbl[k].x==0)vel[k].x=0;
            if(ionmbl[k].y==0)vel[k].y=0;
            if(ionmbl[k].z==0)vel[k].z=0;
        }
        scalevel();
    }
    moveatoms(step);
    twiceKE=GetAtomKE();
    twiceKE = twiceKE * 2;

    tempNow = twiceKE/(3*double(numIon-nfrozen))/K_BOLTZMAN_AU;
    cout<<" start temperature = "<< tempNow/K_BOLTZMAN_AU<< " (k)"<<endl;

    // Set up forces 
    callforce();
	if(STRESS)
	{
		cout<<"output Pressure for check!"<<endl;
		double press;
		for(int i=0;i<3;i++)
			press += stress_lcao(i,i)/3;
		press += twiceKE/3/ucell.omega; //output virtual press = 2/3*Ek/V + sum(sigma[i][i])/3
		double unit_transform = RYDBERG_SI / pow(BOHR_RADIUS_SI,3) * eps8 ;
		cout<<"Virtual Pressure is "<<press*unit_transform<<" Kbar "<<endl;
	}
  
    for(int k=0;k<numIon;k++)
    {
        if(ionmbl[k].x==0)force[k].x=0;
        if(ionmbl[k].y==0)force[k].y=0;
        if(ionmbl[k].z==0)force[k].z=0;
    }

    //cout<<"begin maxForce"<<endl;
    maxForce = MAXVALF();
    cout<<"maxForce: "<<sqrt(maxForce)<<endl; 
    mdenergy=en.etot/2;
    conservedE = Conserved(twiceKE/2, mdenergy);

    cout<< "NVE_STEP "<<" "<<"SystemEnergy"<<" "<< "Conserved"<<" "<< "DeltaE"<<" "<< "Temperature"<<endl;

    // Output the message to the screen.
    cout<<step<<" "<< mdenergy<<" "<< conservedE<<" "
        << mdenergy-oldEtot<<" "<< twiceKE/(3*double(numIon-nfrozen))/K_BOLTZMAN_AU<<endl;

    cout<<"------------------------------------------------------------------------------"<<endl;
    cout<< "MD(NVE) STEP "<< step<<endl;
    cout<<"------------------------------------------------------------------------------"<<endl;
    if (!MY_RANK){
         out<<"------------------------------------------------------------------------------"<<endl;
         out<< "MD(NVE) STEP "<< step<<endl;
         out<<"------------------------------------------------------------------------------"<<endl;
         ofs_running<<"------------------------------------------------------------------------------"<<endl;
         ofs_running<< "MD(NVE) STEP "<< step<<endl;
         ofs_running<<"------------------------------------------------------------------------------"<<endl;
    }
  
    // (1) 1st step of Verlet-Velocity
    // New velocity are obtained
	for(int k=0;k<numIon;k++)
    {       
        mass = allmass[k];
        vel[k] = vel[k] + force[k]/mass*dt/2.0;
    } 
    // (2) 2nd step of Verlet-Velocity 
    // Update the Non-Wrapped cartesion coordinate
    twiceKE = GetAtomKE();
    twiceKE = 2 * twiceKE;
    if(step!=1||rstMD==1)for(int k=0;k<numIon;k++)
    {       
        mass = allmass[k];
        vel[k] = vel[k] + force[k]/mass*dt/2.0;
    } 
    for(int i=0;i<numIon;i++)
    {
        cartNoWrap[i] = vel[i]*dt + cartNoWrap[i];	
	}
    // Calculate the maximal velocities.
    // The step in fractional coordinates
    maxStep = 0;
    for(int  i = 0;i< numIon;i++)
    {
        if((pow(vel[i].x,2.0)+pow(vel[i].y,2.0)+pow(vel[i].z,2.0))>maxStep)
        {
            maxStep = pow(vel[i].x,2.0)+pow(vel[i].y,2.0)+pow(vel[i].z,2.0);
        }
        Mathzone::Cartesian_to_Direct(vel[i].x*dt/ucell.lat0,vel[i].y*dt/ucell.lat0,vel[i].z*dt/ucell.lat0,
                                        ionlatvec.e11,ionlatvec.e12,ionlatvec.e13,
                                        ionlatvec.e21,ionlatvec.e22,ionlatvec.e23,
                                        ionlatvec.e31,ionlatvec.e32,ionlatvec.e33,
                                        fracStep.x,fracStep.y,fracStep.z);
        taudirac[i] = taudirac[i] + fracStep;
    }
    moveatoms(step);
    string t("md_pos_");
    t=intTurnTostring(step,t);
	connection2();
    printpos(t,step);
    maxStep = sqrt(maxStep)*dt;


    // calculate the conserved quantity during MD 
    hamiltonian = Conserved(twiceKE/2, mdenergy);

    cout<< setprecision (9)<<hamiltonian<<" "<< setprecision (9)<<twiceKE/2<<endl;

    
    // Output the message to the screen.
    if (!MY_RANK)
    { 
         out<<step<<" "<< mdenergy<<" "<< hamiltonian<<" "<< mdenergy-oldEtot<<" "<< twiceKE/(3*double(numIon-nfrozen))/K_BOLTZMAN_AU<<endl;
         out<<step<<" "<< mdenergy<<" "<< conservedE<<" "<< mdenergy-oldEtot<<" "<< twiceKE/(3*double(numIon-nfrozen))/K_BOLTZMAN_AU<<endl;
         ofs_running<<step<<" "<< mdenergy<<" "<< hamiltonian<<" "<< mdenergy-oldEtot<<" "<< twiceKE/(3*double(numIon-nfrozen))/K_BOLTZMAN_AU<<endl;
         ofs_running<<step<<" "<< mdenergy<<" "<< conservedE<<" "<< mdenergy-oldEtot<<" "<< twiceKE/(3*double(numIon-nfrozen))/K_BOLTZMAN_AU<<endl;
    }
    oldEtot=mdenergy;
    mstout(step);

    cout<<"(NVE): this step finished."<<endl;
    if (!MY_RANK)
    {
        out<<"(NVE): this step finished."<<endl;
        ofs_running<<"(NVE): this step finished."<<endl;
    }
//2015-09-25, xiaohui
#ifdef __MPI
    if(step%100==0)
    {
        MPI_Bcast(vel,numIon*3,MPI_DOUBLE,0,MPI_COMM_WORLD);
        MPI_Bcast(cartNoWrap,numIon*3,MPI_DOUBLE,0,MPI_COMM_WORLD);
        MPI_Bcast(taudirac,numIon*3,MPI_DOUBLE,0,MPI_COMM_WORLD);
 //       MPI_Bcast(&xLogS,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
 //       MPI_Bcast(&vLogS,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    }
#endif

    timer::tick("MD_run","runNVE",'C');
    return;
}

bool MD_run::runFIRE(int step1){
//-------------------------------------------------------------------------------
// REFERENCES:
//   
    TITLE("MD_run","runFIRE");
    timer::tick("MD_run","runFIRE",'C');
    if(step==1)
    {
        cout<<" FIRE Start."<<endl;
    }
    else 
    {
        cout<<" FIRE step: "<<step<<endl;
    }
  
    cout<<" FIRE Fronzen Atom Number : "<<nfrozen<<endl;
    

    // Set up extra output to ion optimizer / MD header
    cout<<"Time interval       : "<< dt*fundamentalTime*1e15<< " (fs)"<<endl;

    if(step==1)
    {
        for(int k=0;k<numIon;k++)
        {
            if(ionmbl[k].x==0)vel[k].x=0;
            if(ionmbl[k].y==0)vel[k].y=0;
            if(ionmbl[k].z==0)vel[k].z=0;
        }
    //      if(nfrozen==0)RemoveMovementOfCenterOfMass();
        scalevel();
    }
    moveatoms(step);
    twiceKE=GetAtomKE();
    twiceKE = twiceKE * 2;

    tempNow = twiceKE/(3*double(numIon-nfrozen))/K_BOLTZMAN_AU;
    cout<<" start temperature = "<< tempNow/K_BOLTZMAN_AU<< " (k)"<<endl;





    // Set up forces 
    callforce();
    
    for(int k=0;k<numIon;k++)
    {
        if(ionmbl[k].x==0)force[k].x=0;
        if(ionmbl[k].y==0)force[k].y=0;
        if(ionmbl[k].z==0)force[k].z=0;
    }

    //cout<<"begin maxForce"<<endl;
    maxForce = MAXVALF();
    cout<<"maxForce: "<<sqrt(maxForce)<<endl; 
    mdenergy = en.etot/2;
    conservedE = Conserved(twiceKE/2, mdenergy);

    cout<<"------------------------------------------------------------------------------"<<endl;
    cout<< "FIRE STEP "<< step<<endl;
    cout<<"------------------------------------------------------------------------------"<<endl;
    if (!MY_RANK){
         ofs_running<<"------------------------------------------------------------------------------"<<endl;
         ofs_running<< "FIRE STEP "<< step<<endl;
         ofs_running<<"------------------------------------------------------------------------------"<<endl;
    }
  
    // (1) 1st step of Verlet-Velocity
    // New velocity are obtained
	for(int k=0;k<numIon;k++)
    {       
        mass = allmass[k];
        vel[k] = vel[k] + force[k]/mass*dt/2.0;
    } 

	//addedy by zifei
	//initialise FIRE parameters
	if(step1 == 1 ){
		dt_max=2.5*dt;
		alpha_start=0.10;
		alpha = alpha_start;

		finc=1.1;
		fdec=0.5;
		f_alpha=0.99;
		N_min=4;
		negative_count=0;

		ofs_running<< "dt_max = " <<dt_max <<endl;
		ofs_running<< "alpha_start begin = " <<alpha_start <<endl;
		ofs_running<< "alpha begin = " <<alpha <<endl;
		ofs_running<< "finc begin = " <<finc <<endl;
		ofs_running<< "fdec begin = " <<fdec <<endl;
		ofs_running<< "N_min = " <<N_min <<endl;
		ofs_running<< "Negative count = " <<negative_count <<endl;
	}

    // (2) 2nd step of Verlet-Velocity 
    // Update the Non-Wrapped cartesion coordinate
    twiceKE = GetAtomKE();
    twiceKE = 2 * twiceKE;
    if(step!=1)for(int k=0;k<numIon;k++)
    {       
        mass = allmass[k];
        vel[k] = vel[k] + force[k]/mass*dt/2.0;
    } 

	double largest_grad_FIRE = 0.0;
	for(int i=0;i<numIon;i++)
    {
        if(largest_grad_FIRE < abs(force[i].x))
        {
            largest_grad_FIRE= abs(force[i].x);
        } 
        if(largest_grad_FIRE < abs(force[i].y))
        {
            largest_grad_FIRE= abs(force[i].y);
        } 
        if(largest_grad_FIRE < abs(force[i].z))
        {
            largest_grad_FIRE= abs(force[i].z);
        } 
	}
	
	ofs_running << " LARGEST GRAD (eV/A)  : " << largest_grad_FIRE * Ry_to_eV / 0.529177 << endl;

	if(largest_grad_FIRE*Ry_to_eV/0.529177 < 0.01)
    {
	//"convergency reach"
		//cout <<"CONVERGENCY REACH of FIRE in the "<<step <<" steps " <<endl;
		timer::tick("MD_run","runFIRE",'C');
		return true;
	}

	check_FIRE();
    for(int i=0;i<numIon;i++)
    {
        cartNoWrap[i] = vel[i]*dt + cartNoWrap[i];	
	}
    // Calculate the maximal velocities.
    // The step in fractional coordinates
    maxStep = 0;
    for(int  i = 0;i< numIon;i++)
    {
        if((pow(vel[i].x,2.0)+pow(vel[i].y,2.0)+pow(vel[i].z,2.0))>maxStep)
        {
            maxStep = pow(vel[i].x,2.0)+pow(vel[i].y,2.0)+pow(vel[i].z,2.0);
        }
        Mathzone::Cartesian_to_Direct(vel[i].x*dt/ucell.lat0,vel[i].y*dt/ucell.lat0,vel[i].z*dt/ucell.lat0,
                                    ionlatvec.e11,ionlatvec.e12,ionlatvec.e13,
                                    ionlatvec.e21,ionlatvec.e22,ionlatvec.e23,
                                    ionlatvec.e31,ionlatvec.e32,ionlatvec.e33,
                                    fracStep.x,fracStep.y,fracStep.z);
        taudirac[i] = taudirac[i] + fracStep;
    }

    moveatoms(step);
    string t("md_pos_");
    t=intTurnTostring(step,t);
	connection2();
    printpos(t,step);
    maxStep = sqrt(maxStep)*dt;

    hamiltonian = Conserved(twiceKE/2, mdenergy);



    
    // Output the message to the screen.
    oldEtot=mdenergy;

    cout<<"(NVE): this step finished."<<endl;
    if (!MY_RANK){
        ofs_running<<"(NVE): this step finished."<<endl;
    }
//2015-09-25, xiaohui
#ifdef __MPI
    if(step%100==0){
        MPI_Bcast(vel,numIon*3,MPI_DOUBLE,0,MPI_COMM_WORLD);
        MPI_Bcast(cartNoWrap,numIon*3,MPI_DOUBLE,0,MPI_COMM_WORLD);
        MPI_Bcast(taudirac,numIon*3,MPI_DOUBLE,0,MPI_COMM_WORLD);
    }
#endif

    timer::tick("MD_run","runFIRE",'C');
    return false;
}
