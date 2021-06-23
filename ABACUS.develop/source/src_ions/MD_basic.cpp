#include "src_ions/MD_basic.h"
#include "src_pw/global.h"

//define in MD_basic.h
//class MD_basic

const double MD_basic::fundamentalTime = 2.418884326505e-17;

MD_basic::MD_basic(MD_parameters& MD_para_in, UnitCell_pseudo &unit_in):
    mdp(MD_para_in),
    ucell(unit_in)
{	
	mdp.dt=mdp.dt/fundamentalTime/1e15;
	temperature_=mdp.tfirst/3.1577464e5;


	//other parameter default
	outputstressperiod_ = 1 ;// The period to output stress
//	ucell.nat=ucell.nat;
//	ntype=ucell.ntype;
	step_rst_=0;
    step_=0;
//	ucell.latvec=ucell.latvec;

    vel=new Vector3<double>[ucell.nat];
	cart_change=new Vector3<double>[ucell.nat];
	tauDirectChange=new Vector3<double>[ucell.nat];
	allmass=new double[ucell.nat];
	ionmbl=new Vector3<int>[ucell.nat];
	force=new Vector3<double>[ucell.nat];

	energy_=en.etot/2;

	//MD starting setup
	if(!mdp.rstMD){
		nfrozen_ = mdf.getMassMbl(ucell, allmass, ionmbl);
		//mdf.gettauDirectChange(ucell, tauDirectChange);
		//A fresh new MD: Do not restart MD
		mdf.InitVelocity(ucell.nat, temperature_, fundamentalTime, allmass, vel) ;
		// Initialize thermostat, and barostat
		
	}
	else{
		nfrozen_ = mdf.getMassMbl(ucell, allmass, ionmbl);
		//mdf.gettauDirectChange(ucell, tauDirectChange);
		mdf.InitVelocity(ucell.nat, temperature_, fundamentalTime, allmass, vel) ;
		if(!mdf.RestartMD(ucell.nat, vel, step_rst_)){
			cout<<"error in restart MD!"<<endl;
			return;
		}
	}
	if(mdp.NVT_tau<1e-10)
	{
		if(mdp.NVT_control == 1) mdp.NVT_tau = 20 * mdp.dt; //NHC
		else if(mdp.NVT_control == 2) mdp.NVT_tau = 20 * mdp.dt; //LGV
		else if(mdp.NVT_control == 3) mdp.NVT_tau = 20 * mdp.dt; //ADS
		else WARNING_QUIT("md:InitMD", "please choose available reservoir!!!");
	}
	else
	{
		mdp.NVT_tau /= fundamentalTime*1e15;
	}
	//tau = 1.0/(NVT_tau*fundamentalTime*1e15);
	// Q=KbT*tau**2
    mdp.Qmass=mdp.Qmass/6.02/9.109*1e5;
	if(mdp.Qmass<1e-10)//Qmass=dt*fundamentalTime*1e15/6.02/9.109*1e5;
	mdp.Qmass = pow(mdp.NVT_tau,2)*temperature_;///beta;
	//NHC setup
	if(mdp.mdtype==1 && mdp.NVT_control == 1)
	{
		mdt.init_NHC(
        mdp.MNHC, 
        mdp.Qmass, 
        mdp.NVT_tau, 
        mdp.dt,
        mdp.NVT_control, 
        ofs_running, 
        ucell.nat,
        temperature_,
        vel,
        allmass
        );
        if(mdp.rstMD) mdt.NHC_restart();
	}
}

MD_basic::~MD_basic()
{
    delete []force;
    delete []ionmbl;
    delete []cart_change;
    delete []vel;
    delete []tauDirectChange;
    delete []allmass;
}

void MD_basic::runNVT(int step1){
//------------------------------------------------------------------------------
// DESCRIPTION:
// Molecular dynamics calculation with fixed Volume and slight fluctuated temperature
// Using thermostat : 1, Nose-Hoover Chains; 2, Langevin; 3, Anderson
// Normal Nose-Hoover thermostat method is retained for test.
//------------------------------------------------------------------------------

	TITLE("MD_basic","runnvt");
	timer::tick("MD_basic","runnvt",'C');
	step_=step1+step_rst_;
	//the real MD step
	
	//change target temperature
//	if (step_!=1)mdf.ReadNewTemp( step_ );
	
	cout << " ---------------------------------------------------------" << endl;
	
	cout<< " Target temperature  : "<< temperature_/K_BOLTZMAN_AU<< " (K)"<<endl;
	
	if(step_==1||step_%mdp.fixTemperature==1)
    {
		for(int k=0;k<ucell.nat;k++)
        {
			if(ionmbl[k].x==0)vel[k].x=0;
			if(ionmbl[k].y==0)vel[k].y=0;
			if(ionmbl[k].z==0)vel[k].z=0;
		}
		mdf.scalevel(ucell.nat, nfrozen_, temperature_, vel, allmass);
	}

	// get the kinetic energy
	double twiceKE = mdf.GetAtomKE(ucell.nat, vel, allmass);
	twiceKE = twiceKE * 2;
	
	//Set up forces and stress
    if(INPUT.basis_type == "lcao")
    {
	    mdf.callInteraction_LCAO(ucell.nat, force, stress);
    }
    else if(INPUT.basis_type == "pw")
    {
        mdf.callInteraction_PW(ucell.nat, force, stress);
    }
	//print total stress + stress_MD
	if(STRESS)
	{
		outStressMD(stress, twiceKE);
	}

	for(int k=0;k<ucell.nat;k++){
		if(ionmbl[k].x==0)force[k].x=0;
		if(ionmbl[k].y==0)force[k].y=0;
		if(ionmbl[k].z==0)force[k].z=0;
	}
	double maxForce = mdf.MAXVALF(ucell.nat, force);
	energy_=en.etot/2;
    double hamiltonian;
	//----------------------------------------------
	// big loop
	//-----------------------------------------------
	cout<<" "<<std::left<<setw(12)<<"MD_STEP"<<std::left<<setw(12)<< "SystemE"<<std::left<<setw(12)<< "Conserved"<<std::left<<setw(12)<< "DeltaE"<<std::left<<setw(12)<< "Temperature"<<endl;
	
	cout<<" "<<std::left<<setw(12)<<step_<<std::left<<setw(12)<< energy_<<std::left<<setw(12)<< hamiltonian<<std::left<<setw(12)<< energy_-oldEtot_<<std::left<<setw(12)<<twiceKE/(3*double(ucell.nat-nfrozen_))/K_BOLTZMAN_AU<<endl;
	
	oldEtot_=energy_;

	if (!MY_RANK)
    {
        ofs_running<<" --------------------------------------------------"<<endl;
        ofs_running<<" Molecular Dynamics (NVT) STEP "<< step_<<endl;
        ofs_running<< " --------------------------------------------------"<<endl;
	}
	
	
	// Calculate the Mean-Square-Displacement.
	if(step_==1&&mdp.rstMD==0)
    { 
        //Note: side scheme position before
        //Now turn to middle scheme.
	
	    // New vel obtained (Verlet-like step)
        if(mdp.NVT_control == 1)
        {
            hamiltonian = mdt.NHChamiltonian(
                twiceKE/2, 
                energy_,
                temperature_,
                nfrozen_);
        }
        else hamiltonian = mdf.Conserved(twiceKE/2, energy_, ucell.nat-nfrozen_);
	    this->update_half_velocity();
	}
	else 
    {
		this->update_half_velocity();
		
		twiceKE=mdf.GetAtomKE(ucell.nat, vel, allmass);
		twiceKE = 2 * twiceKE;
		if(mdp.NVT_control==1)
        {
            hamiltonian = mdt.NHChamiltonian(
                twiceKE/2, 
                energy_,
                temperature_,
                nfrozen_);
        } 
		else hamiltonian = mdf.Conserved(twiceKE/2, energy_, ucell.nat-nfrozen_);
        //Note: side scheme position before
        //Now turn to middle scheme. 
		this->update_half_velocity();
	}

	// Update the Non-Wrapped cartesion coordinates
    if(mdp.mdtype==2) mdf.scalevel(ucell.nat, nfrozen_, temperature_, vel, allmass);//choose velocity scaling method

    update_half_direct(1);
	
	mdt.Integrator(mdp.NVT_control, temperature_, vel, allmass);//thermostat interact with velocity
    twiceKE=mdf.GetAtomKE(ucell.nat, vel, allmass);
    twiceKE = 2 * twiceKE;
	
	update_half_direct(0);

	// Calculate the maximal velocities.
	// The step in fractional coordinates
	double maxStep = 0;
    for(int  ii=0;ii<ucell.nat;ii++)
    { 
		if((pow(vel[ii].x,2)+pow(vel[ii].y,2)+pow(vel[ii].z,2))>maxStep)
        {
		    maxStep = pow(vel[ii].x,2)+pow(vel[ii].y,2)+pow(vel[ii].z,2);
        }
    }
	getTaudUpdate();

	//save the atom position change to DFT module
	save_output_position();
	maxStep = sqrt(maxStep)*mdp.dt;

    if (!MY_RANK){
        ofs_running<<setw(15)<<"maxForce=     "<<setw(15)<<"maxstep=      "<<setw(15)<<"step=     "<<endl;
        ofs_running<<setw(15)<<maxForce<<setw(15)<<maxStep<<setw(15)<<step_<<endl;
        ofs_running<<step_<<" "<< energy_<<" "<< hamiltonian<<" "<< energy_-oldEtot_<<" "<< twiceKE/(3*double(ucell.nat-nfrozen_))/K_BOLTZMAN_AU<<endl;
    }
    oldEtot_=energy_;
    //output basic restart info
    mdf.mdRestartOut(step_, mdp.recordFreq, ucell.nat, vel);
    //output NHC thermo restart info
    if(mdp.NVT_control==1) mdt.NHC_info_out(step_, mdp.recordFreq, ucell.nat);

//2015-09-25, xiaohui
#ifdef __MPI
    if(step_%100==0)
    {
        MPI_Bcast(vel,ucell.nat*3,MPI_DOUBLE,0,MPI_COMM_WORLD);
    //    MPI_Bcast(&xLogS,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    //    MPI_Bcast(&vLogS,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    }
#endif
    timer::tick("MD_basic","runnvt",'C');
    return;
}

void MD_basic::runNVE(int step1){
//-------------------------------------------------------------------------------
// Adiabatic ensemble 
// Molecular dynamics calculation with Verlet algorithm
//-------------------------------------------------------------------------------

    TITLE("MD_basic","runNVE");
    timer::tick("MD_basic","runNVE",'C');
    step_=step1+step_rst_;
    if(step_==1)
    {
        cout<<" (NVE) Start."<<endl;
    }
    else 
    {
        cout<<" (NVE) step: "<<step_<<endl;
    }

  // Set up extra output to ion optimizer / MD header
  cout<<"Time interval       : "<< mdp.dt*fundamentalTime*1e15<< " (fs)"<<endl;
  cout<<"Target temperature  : "<< temperature_/K_BOLTZMAN_AU<< " (K)"<<endl;

    if(step_==1){
        for(int k=0;k<ucell.nat;k++)
        {
            if(ionmbl[k].x==0)vel[k].x=0;
            if(ionmbl[k].y==0)vel[k].y=0;
            if(ionmbl[k].z==0)vel[k].z=0;
        }
        mdf.scalevel(ucell.nat, nfrozen_, temperature_, vel, allmass);
    }
    double twiceKE=mdf.GetAtomKE(ucell.nat, vel, allmass);
    twiceKE = twiceKE * 2;

    double tempNow = twiceKE/(3*double(ucell.nat-nfrozen_))/K_BOLTZMAN_AU;
    cout<<" start temperature = "<< tempNow/K_BOLTZMAN_AU<< " (k)"<<endl;

    // Set up forces
    if (INPUT.basis_type == "lcao")
    {
        mdf.callInteraction_LCAO(ucell.nat, force, stress);
    }
    else if (INPUT.basis_type == "pw")
    {
        mdf.callInteraction_PW(ucell.nat, force, stress);
    }
    if(STRESS)
	{
		outStressMD(stress, twiceKE);
	}
  
    for(int k=0;k<ucell.nat;k++)
    {
        if(ionmbl[k].x==0)force[k].x=0;
        if(ionmbl[k].y==0)force[k].y=0;
        if(ionmbl[k].z==0)force[k].z=0;
    }

    //cout<<"begin maxForce"<<endl;
    double maxForce = mdf.MAXVALF(ucell.nat, force);
    cout<<"maxForce: "<<sqrt(maxForce)<<endl; 
    energy_=en.etot/2;
    double conservedE = mdf.Conserved(twiceKE/2, energy_, ucell.nat-nfrozen_);

    cout<< "NVE_STEP "<<" "<<"SystemEnergy"<<" "<< "Conserved"<<" "<< "DeltaE"<<" "<< "Temperature"<<endl;

    // Output the message to the screen.
    cout<<step_<<" "<< energy_<<" "<< conservedE<<" "
        << energy_-oldEtot_<<" "<< twiceKE/(3*double(ucell.nat-nfrozen_))/K_BOLTZMAN_AU<<endl;

    cout<<"------------------------------------------------------------------------------"<<endl;
    cout<< "MD(NVE) STEP "<< step_<<endl;
    cout<<"------------------------------------------------------------------------------"<<endl;
    if (!MY_RANK){
        ofs_running<<"------------------------------------------------------------------------------"<<endl;
        ofs_running<< "MD(NVE) STEP "<< step_<<endl;
        ofs_running<<"------------------------------------------------------------------------------"<<endl;
    }
  
    // (1) 1st step of Verlet-Velocity
    // New velocity are obtained
	this->update_half_velocity();
    // (2) 2nd step of Verlet-Velocity 
    // Update the Non-Wrapped cartesion coordinate
    twiceKE = mdf.GetAtomKE(ucell.nat, vel, allmass);
    twiceKE = 2 * twiceKE;
    if(step_!=1||mdp.rstMD==1)this->update_half_velocity();
    update_half_direct(1);
    update_half_direct(0);
    // Calculate the maximal velocities.
    // The step in fractional coordinates
    double maxStep = 0;
    for(int  i = 0;i< ucell.nat;i++)
    {
        if((pow(vel[i].x,2.0)+pow(vel[i].y,2.0)+pow(vel[i].z,2.0))>maxStep)
        {
            maxStep = pow(vel[i].x,2.0)+pow(vel[i].y,2.0)+pow(vel[i].z,2.0);
        }
    }
    getTaudUpdate();
    save_output_position();
    maxStep = sqrt(maxStep)*mdp.dt;


    // calculate the conserved quantity during MD 
    double hamiltonian = mdf.Conserved(twiceKE/2, energy_, ucell.nat-nfrozen_);

    cout<< setprecision (9)<<hamiltonian<<" "<< setprecision (9)<<twiceKE/2<<endl;

    
    // Output the message to the screen.
    if (!MY_RANK)
    { 
        ofs_running<<step_<<" "<< energy_<<" "<< hamiltonian<<" "<< energy_-oldEtot_<<" "<< twiceKE/(3*double(ucell.nat-nfrozen_))/K_BOLTZMAN_AU<<endl;
        ofs_running<<step_<<" "<< energy_<<" "<< conservedE<<" "<< energy_-oldEtot_<<" "<< twiceKE/(3*double(ucell.nat-nfrozen_))/K_BOLTZMAN_AU<<endl;
    }
    oldEtot_=energy_;
    mdf.mdRestartOut(step_, mdp.recordFreq, ucell.nat, vel);

    cout<<"(NVE): this step finished."<<endl;
    if (!MY_RANK)
    {
        ofs_running<<"(NVE): this step finished."<<endl;
    }
//2015-09-25, xiaohui
#ifdef __MPI
    if(step_%100==0)
    {
        MPI_Bcast(vel,ucell.nat*3,MPI_DOUBLE,0,MPI_COMM_WORLD);
 //       MPI_Bcast(&xLogS,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
 //       MPI_Bcast(&vLogS,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    }
#endif

    timer::tick("MD_basic","runNVE",'C');
    return;
}

bool MD_basic::runFIRE(int step1){
//-------------------------------------------------------------------------------
// REFERENCES:
//   
    TITLE("MD_basic","runFIRE");
    timer::tick("MD_basic","runFIRE",'C');
    step_ = step1;
    if(step_==1)
    {
        cout<<" FIRE Start."<<endl;
    }
    else 
    {
        cout<<" FIRE step: "<<step_<<endl;
    }
  
    cout<<" FIRE Fronzen Atom Number : "<<nfrozen_<<endl;
    

    // Set up extra output to ion optimizer / MD header
    cout<<"Time interval       : "<< mdp.dt*fundamentalTime*1e15<< " (fs)"<<endl;

    if(step_==1)
    {
        for(int k=0;k<ucell.nat;k++)
        {
            if(ionmbl[k].x==0)vel[k].x=0;
            if(ionmbl[k].y==0)vel[k].y=0;
            if(ionmbl[k].z==0)vel[k].z=0;
        }
    //      if(nfrozen_==0)RemoveMovementOfCenterOfMass();
        mdf.scalevel(ucell.nat, nfrozen_, temperature_, vel, allmass);
    }
    double twiceKE=mdf.GetAtomKE(ucell.nat, vel, allmass);
    twiceKE = twiceKE * 2;

    double tempNow = twiceKE/(3*double(ucell.nat-nfrozen_))/K_BOLTZMAN_AU;
    cout<<" start temperature = "<< tempNow/K_BOLTZMAN_AU<< " (k)"<<endl;





    // Set up forces
    if (INPUT.basis_type == "lcao")
    {
        mdf.callInteraction_LCAO(ucell.nat, force, stress);
    }
    else if (INPUT.basis_type == "pw")
    {
        mdf.callInteraction_PW(ucell.nat, force, stress);
    }

    for(int k=0;k<ucell.nat;k++)
    {
        if(ionmbl[k].x==0)force[k].x=0;
        if(ionmbl[k].y==0)force[k].y=0;
        if(ionmbl[k].z==0)force[k].z=0;
    }

    //cout<<"begin maxForce"<<endl;
    double maxForce = mdf.MAXVALF(ucell.nat, force);
    cout<<"maxForce: "<<sqrt(maxForce)<<endl; 
    energy_ = en.etot/2;
    double conservedE = mdf.Conserved(twiceKE/2, energy_, ucell.nat-nfrozen_);

    cout<<"------------------------------------------------------------------------------"<<endl;
    cout<< "FIRE STEP "<< step_<<endl;
    cout<<"------------------------------------------------------------------------------"<<endl;
    if (!MY_RANK){
         ofs_running<<"------------------------------------------------------------------------------"<<endl;
         ofs_running<< "FIRE STEP "<< step_<<endl;
         ofs_running<<"------------------------------------------------------------------------------"<<endl;
    }
  
    // (1) 1st step of Verlet-Velocity
    // New velocity are obtained
	this->update_half_velocity();

    // (2) 2nd step of Verlet-Velocity 
    // Update the Non-Wrapped cartesion coordinate
    twiceKE = mdf.GetAtomKE(ucell.nat, vel, allmass);
    twiceKE = 2 * twiceKE;
    if(step_!=1)this->update_half_velocity();

	double largest_grad_FIRE = 0.0;
	for(int i=0;i<ucell.nat;i++)
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
		timer::tick("MD_basic","runFIRE",'C');
		return true;
	}

	fire.check_FIRE(ucell.nat, force, mdp.dt, vel);
    
    update_half_direct(1);
    update_half_direct(0);
    // Calculate the maximal velocities.
    // The step in fractional coordinates
    double maxStep = 0;
    for(int  i = 0;i< ucell.nat;i++)
    {
        if((pow(vel[i].x,2.0)+pow(vel[i].y,2.0)+pow(vel[i].z,2.0))>maxStep)
        {
            maxStep = pow(vel[i].x,2.0)+pow(vel[i].y,2.0)+pow(vel[i].z,2.0);
        }
    }
    getTaudUpdate();

    save_output_position();
    maxStep = sqrt(maxStep)*mdp.dt;

    double hamiltonian = mdf.Conserved(twiceKE/2, energy_, ucell.nat-nfrozen_);



    
    // Output the message to the screen.
    oldEtot_=energy_;

    cout<<"(NVE): this step finished."<<endl;
    if (!MY_RANK){
        ofs_running<<"(NVE): this step finished."<<endl;
    }
//2015-09-25, xiaohui
#ifdef __MPI
    if(step1%100==0){
        MPI_Bcast(vel,ucell.nat*3,MPI_DOUBLE,0,MPI_COMM_WORLD);
    }
#endif

    timer::tick("MD_basic","runFIRE",'C');
    return false;
}

//update velocities of ions for half MD step
void MD_basic::update_half_velocity()
{
    for(int  ii=0;ii<ucell.nat;ii++){ 
        vel[ii] = vel[ii] + force[ii]/allmass[ii]*mdp.dt/2.0;
    }
}

//update cartesian coordinate changes for half MD step, would turn to direct coordinate changes later
void MD_basic::update_half_direct(const bool is_restart)
{
    if(is_restart)
    {//first half
        for(int  ii=0;ii<ucell.nat;ii++)
        { 
            cart_change[ii] = vel[ii]*mdp.dt/2.0/ucell.lat0;
        }
    }
    else
    {//second half
        for(int  ii=0;ii<ucell.nat;ii++)
        { 
            cart_change[ii] += vel[ii]*mdp.dt/2.0/ucell.lat0;
        }
    }
}

//output Structure information for each MD step
void MD_basic::save_output_position()
{
    string posOutName("md_pos_");
	posOutName=mdf.intTurnTostring(step_,posOutName);
	ucell.update_pos_taud(tauDirectChange);
	mdf.printpos(posOutName, step_, mdp.recordFreq, ucell);
}

int MD_basic::getRealStep()
{
    return this->step_;
}

//output pressure of total MD system, P = tr(stress) + P_kin
void MD_basic::outStressMD(const matrix& stress, const double& twiceKE)
{
    ofs_running<<"output Pressure for check!"<<endl;
    double press;
    for(int i=0;i<3;i++)
        press += stress(i,i)/3;
    press += twiceKE/3/ucell.omega; //output virtual press = 2/3 *Ek/V + sum(sigma[i][i])/3
    double unit_transform = RYDBERG_SI / pow(BOHR_RADIUS_SI,3) * 1.0e-8 ;
    ofs_running<<"Virtual Pressure is "<<press*unit_transform<<" Kbar "<<endl;
}

//turn cartesian coordinate changes to direct changes
void MD_basic::getTaudUpdate()
{
    Vector3<double> fracStep;
	for(int  ii=0;ii<ucell.nat;ii++)
    { 
		Mathzone::Cartesian_to_Direct(cart_change[ii].x,cart_change[ii].y,cart_change[ii].z,
					ucell.latvec.e11,ucell.latvec.e12,ucell.latvec.e13,
					ucell.latvec.e21,ucell.latvec.e22,ucell.latvec.e23,
					ucell.latvec.e31,ucell.latvec.e32,ucell.latvec.e33,
					fracStep.x,fracStep.y,fracStep.z);

		tauDirectChange[ii] = fracStep;
	}
}
