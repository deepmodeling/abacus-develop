#include "MD_func.h"



void MD_func::initMD(){
	//parameter from input file
	mdtype=INPUT.md_mdtype;
	Qmass=INPUT.md_qmass/6.02/9.109*1e5;
	dt=INPUT.md_dt/fundamentalTime/1e15;
	temperature=INPUT.md_tfirst/3.1577464e5;
	rstMD=INPUT.md_rstmd;
	dumpmdfreq=INPUT.md_dumpmdfred;
	fixTemperature=INPUT.md_fixtemperature;
	ediff=INPUT.md_ediff;
	ediffg=INPUT.md_ediffg;
	mdoutputpath=INPUT.md_mdoutpath;
	NVT_control = INPUT.md_nvtcontrol;
   	NVT_tau = INPUT.md_nvttau;
	MNHC = INPUT.md_mnhc;


	//other parameter default
	outputstressperiod = 1 ;// The period to output stress
	numIon=ucell.nat;
	ntype=ucell.ntype;
	step_rst=0;
	ionlatvec=ucell.latvec;

	mdenergy=en.etot/2;

	//allocate for MD
	this->md_allocate();

	//MD starting setup
	if(rstMD==0){
		connection0();
		connection1();
		//A fresh new MD: Do not restart MD
		InitVelocity() ;
		// Initialize thermostat, and barostat
		
	}
	else if ( rstMD>0 ){
		connection0();
		connection1();
		InitVelocity() ;
		if(!RestartMD()){
			cout<<"error in restart MD!"<<endl;
			return;
		}
		connection2();
	}
	if(NVT_tau<1e-10)
	{
		if(NVT_control == 1) NVT_tau = 20 * dt; //NHC
		else if(NVT_control == 2) NVT_tau = 20 * dt; //LGV
		else if(NVT_control == 3) NVT_tau = 20 * dt; //ADS
		else WARNING_QUIT("md:InitMD", "please choose available reservoir!!!");
	}
	else
	{
		NVT_tau *= fundamentalTime*1e15;
	}
	//tau = 1.0/(NVT_tau*fundamentalTime*1e15);
	// Q=KbT*tau**2
	if(Qmass<1e-10)//Qmass=dt*fundamentalTime*1e15/6.02/9.109*1e5;
	Qmass = pow(1.0/NVT_tau,2)*temperature;///beta;
	//NHC setup
	if(mdtype==1 && NVT_control == 1)
	{
		unsigned long init[4]={0x123, 0x234, 0x345, 0x456}, length=4;
		init_by_array(init, length);
		ofs_running<<" ...............Nose-Hoover Chain parameter initialization...............  " << endl;
		ofs_running<<" Temperature =    "<< temperature << endl;
		ofs_running<<" Temperature2 =    "<< temperature/K_BOLTZMAN_AU << endl;
		ofs_running<<" NHC frequency =    "<< 1.0/NVT_tau << endl;
		ofs_running<<" NHC chain =    "<< MNHC << endl;
		ofs_running<<" Qmass  =    "<< Qmass << endl;
		ofs_running<<" ...............................................................  " << endl;
		w[0]=0.7845136105;
		w[6]=0.7845136105;
		w[1]=0.2355732134;
		w[5]=0.2355732134;
		w[2]=-1.177679984;
		w[4]=-1.177679984;
		w[3]=1-w[0]-w[1]-w[2]-w[4]-w[5]-w[6];
		//
		for(int i=0;i<nsy;i++){
			delta[i]=w[i]*dt;
		}
		
		
		for(int j=0;j<MNHC;j++){
			for(int k=0;k<numIon;k++){
				NHCeta[j*numIon+k].x=genrand_real2()-0.5; 
				NHCeta[j*numIon+k].y=genrand_real2()-0.5; 
				NHCeta[j*numIon+k].z=genrand_real2()-0.5;
				NHCpeta[j*numIon+k].x=gaussrand()*sqrt(Qmass * temperature);
				NHCpeta[j*numIon+k].y=gaussrand()*sqrt(Qmass * temperature);
				NHCpeta[j*numIon+k].z=gaussrand()*sqrt(Qmass * temperature);
			}
		}
		
		for(int j=MNHC-1;j>=1;j--){
			for(int k=0;k<numIon;k++){
				G[j*numIon+k].x=pow(NHCpeta[(j-1)*numIon+k].x,2)/Qmass - 1.0 * temperature;
				G[j*numIon+k].y=pow(NHCpeta[(j-1)*numIon+k].y,2)/Qmass - 1.0 * temperature;
				G[j*numIon+k].z=pow(NHCpeta[(j-1)*numIon+k].z,2)/Qmass - 1.0 * temperature;
			}
		}
		
		for(int k=0;k<numIon;k++){
			G[0*numIon+k].x=pow(vel[k].x,2)*allmass[k]-1.0 * temperature;
			G[0*numIon+k].y=pow(vel[k].y,2)*allmass[k]-1.0 * temperature;
			G[0*numIon+k].z=pow(vel[k].z,2)*allmass[k]-1.0 * temperature;
		}
		ofs_running << "finish NHC thermostat define" << endl; //zifei
	}

	return;
}
	
/*void md::runMD(int step1){
       if(rstMD==1)MDNVT.runnvt();
       else if(mdtype==3)MDNVE.runNVE(step1);
       return;
}*///we need this form

bool  MD_func::RestartMD(){
	int i;
	double *vell=new double[numIon*3];
	double *cart=new double[numIon*3];
	double *dira=new double[numIon*3];
	double *nhcg=new double[numIon*3*MNHC];
	double *nhcp=new double[numIon*3*MNHC];
	double *nhce=new double[numIon*3*MNHC];
	if (!MY_RANK){
		stringstream ssc;
		ssc << global_out_dir << "Restart_md.dat";
		ifstream file(ssc.str().c_str());

		if(!file){
			cout<<"please ensure whether 'Restart_md.dat' exists"<<endl;
			exit(0);
		}
		if (file.good())
		{
			file.ignore(11, '\n');
		}

//----------------------------------------------------------
// main parameters
//----------------------------------------------------------
		file.ignore(7, '\n');
		//                file>>xLogS;
		file.get();	
		file.ignore(7, '\n');
		//                file>>vLogS;
		file.get();
		file.ignore(23, '\n');
		for(i = 0;i<numIon*3;i++){
			file>>vell[i];
		}
		/*file.get();	file.ignore(17, '\n');
		for(i=0;i<numIon*3;i++){
			file>>cart[i];
		}
		file.get();	file.ignore(13, '\n');
		for(i=0;i<numIon*3;i++){
			file>>dira[i];
		}*/
		file.get();     
		file.ignore(6, '\n');
		file>>step_rst;
		file.close();
	}
//2015-09-05, xiaohui
#ifdef __MPI
	MPI_Bcast(&step_rst,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(vell,numIon*3,MPI_DOUBLE,0,MPI_COMM_WORLD);
#endif
	for(i=0;i<numIon;i++){
		vel[i].x=vell[i*3];
		vel[i].y=vell[i*3+1];
		vel[i].z=vell[i*3+2];
	}
	// DEBUG
	//int debug = 0;
	//if (debug) {
	//  cout<<" dump cell + ion infomation file to -> echo_md_restart.dat, for debug purpose. "<<endl;
	//	fstream debugfile;
	//	debugfile.open("echo_md_restart.dat");
	//  debugfile<<"  "<<endl;
	//debugfile<<"cellReal[0] = "<<a.latvec.e11<<" "<<a.latvec.e12<<" "<<a.latvec.e13<<endl;
	//debugfile<<"cellReal[1] = "<< a.latvec.e21<<" "<<a.latvec.e22<<" "<<a.latvec.e23<<endl;
	//debugfile<<"cellReal[2] = "<< a.latvec.e31<<" "<<a.latvec.e32<<" "<<a.latvec.e33<<endl;
	//debugfile<<"ION FRACTIONAL COORDINATES:"<<endl;
	//for(ii=0;ii<numIon;ii++){
	//  debugfile<<taudirac[ii].x<<" "<<taudirac[ii].y<<" "<<taudirac[ii].z<<endl;
	//}
	//debugfile<<"ION VELOCITIES:"<<endl;
	//for( ii=0;ii<ionIon;ii++){
	// debugfile<<vel[ii].x<<" "<<vel[ii].y<<" "<<vel[ii].z<<endl;
	//	}
	//  debugfile.close();
	//}
	// END DEBUG
	if(vell!=NULL)delete []vell;
	if(cart!=NULL)delete []cart;
	if(dira!=NULL)delete []dira;
	if(nhcg!=NULL)delete []nhcg;
	if(nhcp!=NULL)delete []nhcp;
	if(nhce!=NULL)delete []nhce;
	return true;
}

void MD_func::mstout(int step){
//this function used for outputting the information of restart file
	bool pass;
	pass = 0;

	if (dumpmdfreq==1||step==1||( dumpmdfreq > 1&& step%dumpmdfreq==0 ) )
		pass =1;
	if (!pass) return;

	if(!MY_RANK){
		stringstream ssc;
		ssc << global_out_dir << "Restart_md.dat";
		ofstream file(ssc.str().c_str());
		file<<"MD_RESTART"<<endl;
		file<<"ION_VELOCITIES_(a.u.): "<<endl;
		for(int i=0;i<numIon;i++){
			file<<setprecision (12)<<vel[i].x<<" "<<setprecision (12)<<vel[i].y<<" "<<setprecision (12)<<vel[i].z<<endl;
		}
/*	file<<"CARTESIAN_COORD: "<<endl;
	for(i=0;i<numIon;i++){
		file<<cartNoWrap[i].x<<" "<<cartNoWrap[i].y<<" "<<cartNoWrap[i].z<<endl;
	}
	file<<"DIRAC_COORD: "<<endl;
	for(i=0;i<numIon;i++){
		file<<taudirac[i].x<<" "<<taudirac[i].y<<" "<<taudirac[i].z<<endl;
	}*/
		file<<"step: "<<step<<endl;                
		file.close();
	}
	return;
}

double MD_func::GetAtomKE(){
//---------------------------------------------------------------------------
// DESCRIPTION:
//   This function calculates the classical kinetic energy of a group of atoms.
//----------------------------------------------------------------------------

	double mass;
	double ke = 0.0;

	// kinetic energy = \sum_{i} 1/2*m_{i}*(vx^2+vy^2+vz^2)
	for(int ion=0;ion<numIon;ion++){
		mass = allmass[ion];
		ke +=   0.5 * mass * (vel[ion].x*vel[ion].x+vel[ion].y*vel[ion].y+vel[ion].z*vel[ion].z);
	}
	//cout<<"in GetAtomKE KE="<< ke<<endl;
	return ke;
}

void MD_func::InitVelocity()
{
	if(!MY_RANK){ //xiaohui add 2015-09-25
		srand(time(0));
		int iy=rand()%21,im=rand()%11,id=rand()%21,ih=rand()%31,in=rand()%61,is=rand()%41;
		int jy=rand()%21,jm=rand()%11,jd=rand()%21,jh=rand()%31,jn=rand()%61,js=rand()%41;
		//    int iy=13,im=10,id=1,ih=21,in=54,is=23;
		//  int jy=14,jm=6,jd=19,jh=4,jn=29,js=32;
		long int i,j;
		long int i0,i1,i2,j0,j1,j2;
		long int a=16807,m=2147483647,q=127773,r=2836;
		long int M;  
		double x,y,u,v;
		i=iy+70*(im+12*(id+31*(ih+23*(in+59*is))));				//initial the seeds
		j=jy+70*(jm+12*(jd+31*(jh+23*(jn+59*js))));	
		i0=a*(i%q)-r*((i-i%q)/q);
		j0=a*(j%q)-r*((j-j%q)/q);
		if(i0>=0)	i1=i0;  	else i1=i0+m;	
		i0=a*(j%q)-r*((j-j%q)/q);
		if(j0>=0)	j1=j0;  	else j1=j0+m;	
		for(M=0;M<3*numIon;){	
			i0=a*(i1%q)-r*((i1-i1%q)/q);                        
			if(i0>=0)	i2=i0;  
			else i2=i0+m;									
			u=((double)i1)/((double)m);		       //first ramdon number        
			i1=i2;
			j0=a*(j1%q)-r*((j1-j1%q)/q);
			if(j0>=0)	j2=j0;  
			else j2=j0+m;										
			v=((double)j1)/((double)m);		        //second ramdon number
			j1=j2;
			x=tan(PI*(u-0.5));
			y=1.6*v/(PI*(x*x+1.0));
			if(y<=(1.0/sqrt(2.0*PI)*exp(-0.5*x*x))){
				if(M<numIon) vel[M].x=x;
				else if(M<2*numIon) vel[M-numIon].y=x;
				else vel[M-2*numIon].z=x;
					// cout<<"this vel is:"<<x<<endl;
					M++;
			}
		}
      
		int ion;
		for(ion=0;ion<numIon;ion++){
			vel[ion]/=sqrt(allmass[ion]*9.01938e-31);				 
		}
		for(ion=0;ion<numIon;ion++){
			vel[ion]*=sqrt(temperature/K_BOLTZMAN_AU * K_BOLTZMAN_SI);
			vel[ion]*= fundamentalTime/BOHR_RADIUS_SI;
		//rescale the velocities to a.u.
		}
	} //xiaohui add 2 lines 2015-09-25, bcast vel
//2015-09-25, xiaohui
#ifdef __MPI
	MPI_Bcast(vel,numIon*3,MPI_DOUBLE,0,MPI_COMM_WORLD);
#endif
	return;
}


void MD_func::OutputMDHeader()
{
//---------------------------------------------------
// Output the Header for MD simulation
//---------------------------------------------------

	// Set up extra output to ion optimizer / MD header
	if (!MY_RANK)
		out<<"Q ="<<Qmass<<"a.u.\n"<<" dt = "<<dt*fundamentalTime*1e15<<"fs\n"<< "T ="<<temperature/K_BOLTZMAN_AU<<endl;
		ofs_running<<"Q ="<<Qmass<<"a.u.\n"<<" dt = "<<dt*fundamentalTime*1e15<<"fs\n"<< "T ="<<temperature/K_BOLTZMAN_AU<<endl;
		return;
} 

void MD_func::OutputMDGeom(int iter)
{
//---------------------------------------------------------------------------
// DESCRIPTION:
//  Output the ion's position and cell's structure for monitoring and MD
//  restart.
//---------------------------------------------------------------------------
  
	bool pass;
	int i,it;
	ofstream uf;
	string iteration_number;
	string eleName;
	   
	//! >> FUNCTION << !
	pass = 0;
	// if(rankGlobal!=0) return;
	if(dumpmdfreq==1||iter==1||( dumpmdfreq > 1&&(iter%dumpmdfreq)==0 ) ) 
		pass = 1;
	if(!pass)return;
	// cout<<"test the function OutputMDGeom"<<endl;

	//the final file name is the md_geometry plus the iteraction number
	iteration_number=intTurnTostring(iter,iteration_number);
	string  message=mdoutputpath+"_ion"+iteration_number+".txt";
	//  cout<<"file name is:"<<message<<endl;
	// open the file: 'ion.1.dat'
	uf.open(message.c_str());
	//  cout<<"already open file"<<endl;
	uf<< "BLOCK LATTICE_CART"<<endl;
	uf<<"    "<<ionlatvec.e11*ucell.lat0<<"  "<<ionlatvec.e12*ucell.lat0<<"  "<<ionlatvec.e13*ucell.lat0<<endl;
	  
	uf<<"    "<<ionlatvec.e21*ucell.lat0<<"  "<<ionlatvec.e22*ucell.lat0<<"  "<<ionlatvec.e23*ucell.lat0<<endl;

	uf<<"    "<<ionlatvec.e31*ucell.lat0<<"  "<<ionlatvec.e32*ucell.lat0<<"  "<<ionlatvec.e33*ucell.lat0<<endl;

	uf<< "END BLOCK LATTICE_CART"<<endl;
	//  cout<<"already out put lattice cart"<<endl;
  
	//output the atom positions.
	uf<<"BLOCK POSITIONS_CART"<<endl;
	int ion=0;
	for(it=0;it<ntype;it++){
		for(i = 0;i< na[it];i++){
			eleName=ucell.atoms[it].label; 
			uf<<"  "<<eleName 
				<<"  "<<cartNoWrap[ion].x  // unit of cartNorap is 'Bohr' 
				<<"  "<<cartNoWrap[ion].y   
				<<"  "<<cartNoWrap[ion++].z<<endl;
		}
	}
	uf<< "END BLOCK POSITIONS_CART"<<endl;
	uf<<"already output positions cart"<<endl;
	double distance;
	distance=pow(cartNoWrap[1].x-cartNoWrap[0].x,2.0)+pow(cartNoWrap[1].y-cartNoWrap[0].y,2.0)+pow(cartNoWrap[1].z-cartNoWrap[0].z,2);
	distance=sqrt(distance);
	//  cout<<"distance is "<<distance<<" ";
	uf.close();
	return;
}

void MD_func::ReadNewTemp(int step )
{
//-----------------------------------------------------------------
//  If fixTemperature == 0, then this subroutine will be skipped
//  otherwise, we are going to read a new temperature from the disk.
//  You only need to create a file named ChangeTemp.dat, and then 
//  change the temperature in it. 
//-----------------------------------------------------------------

	double intemp;
	int sstep=0;

	if ( fixTemperature > 0 )
	{

	// Change the temperature every 'fixTemperature' steps.
		if( (step!=1)&&(step%fixTemperature == 1) )
		{
			
			// Read in new temperature from file.
			ifstream file;
			file.open("ChangeTemp.dat");
			if (!file){
				cout<<"ERROR IN OPENING ChangeTemp.dat, CODE STOP!"<<endl;
				exit(0);
			}
			while((sstep+fixTemperature)<step){
				file>> intemp;
				sstep+=fixTemperature;
			}
			file.close();

			// Renew information.
			intemp =  intemp * K_BOLTZMAN_AU;
			if ( fabs(intemp-temperature) >1e-6 ) {
				cout <<"(ReadNewTemp): Read in new temp:"<< intemp/K_BOLTZMAN_AU 
					<<" previous temp:"<< temperature/K_BOLTZMAN_AU<<endl;
				temperature = intemp;
			}
			else{
				cout<<"(ReadNewTemp): new temp:"<< intemp/K_BOLTZMAN_AU
					<<" previous temp:"<<temperature/K_BOLTZMAN_AU
					<< ". No change of temp."<<endl;
			}
		}
	}

	return;
}

//int to string and add to output path
string MD_func::intTurnTostring(long int iter, string path)
{
	long int i[10],k=0;
	if(iter>9999999999) return "error!";
	for(int j=9; j>-1; j--)
	{
		if(iter==0) continue;
		if(iter>pow(10,j)-1)
		{
			i[k] = iter%10;
			iter /= 10;
			k++;
		}
	}
	for(int j=k-1;j>-1;j--){
		path+=(i[j]+48);
	}
	return path;
}

void MD_func::connection0(){
//some prepared information
//mass and degree of freedom
	int ion=0;
	nfrozen=0;
	for(int it=0;it<ntype;it++){
		for(int i=0;i<na[it];i++){
			allmass[ion]=ucell.atoms[it].mass/6.02/9.109*1e5;
			ionmbl[ion]=ucell.atoms[it].mbl[i];
			if (ionmbl[ion].x==0||ionmbl[ion].y==0||ionmbl[ion].z==0) nfrozen++;

			ion++;
		}
	}
	return;
}
	
void MD_func::connection1()
{
//call each atom's position
	int i,it,ion;
	ion=0;
	for(it=0;it<ntype;it++){
		for(i=0;i<na[it];i++){
			cartNoWrap[ion]=ucell.atoms[it].tau[i]*ucell.lat0;
			taudirac[ion]=ucell.atoms[it].taud[i];
			ion++;
		}
	}
	return;
}

void MD_func::connection2()
{
//to refresh the atoms' positions
	int i,it,ion;
	ion=0;
	for(it=0;it<ntype;it++){
		for(i=0;i<na[it];i++){
			ucell.atoms[it].tau[i]=cartNoWrap[ion]/ucell.lat0;
			ucell.atoms[it].taud[i]=taudirac[ion];
			ion++;
		}
	}
	return;
}

void MD_func::callforce()
{
//to call the force of each atom
	int ion;
	Force_LCAO FL;
	FL.allocate (); 
	FL.start_force();
	for(ion=0;ion<numIon;ion++){
		force[ion].x=FL.fcs(ion,0)/2.0;
		force[ion].y=FL.fcs(ion,1)/2.0;
		force[ion].z=FL.fcs(ion,2)/2.0;
		// cout<<force[ion].x/0.0195<<" "<<force[ion].y/0.0195<<" "<<force[ion].z/0.0195<<endl;
	}
#ifdef __MPI //2015-10-01, xiaohui
	atom_arrange::delete_vector( SEARCH_RADIUS );
#endif //2015-10-01, xiaohui
	//to call the stress
	if(STRESS)
	{
	//	matrix stress_lcao;//this is the stress matrix same as src_pw/ion.cpp
		stress_lcao.create(3,3);
		FL.cal_stress(stress_lcao);
	}
	return;
}

void MD_func::moveatoms(int step){
//intent to ensure that no atom is out of the cell

	for(int ia=0;ia<numIon;ia++)
	{
//changed by a constrain box
                 
		if(taudirac[ia].x<0) taudirac[ia].x += 1.0;
		if(taudirac[ia].y<0) taudirac[ia].y += 1.0;
		if(taudirac[ia].z<0) taudirac[ia].z += 1.0;
		if(taudirac[ia].x>=1.0) taudirac[ia].x -= 1.0;
		if(taudirac[ia].y>=1.0) taudirac[ia].y -= 1.0;
		if(taudirac[ia].z>=1.0) taudirac[ia].z -= 1.0;

		if(taudirac[ia].x<0 || taudirac[ia].y<0
			|| taudirac[ia].z<0 ||
			taudirac[ia].x>=1.0 ||
			taudirac[ia].y>=1.0 ||
			taudirac[ia].z>=1.0)
		{
			cout << " ion: "<< ia  << endl;
			cout << "d=" << taudirac[ia].x << " " <<
				taudirac[ia].y << " " << taudirac[ia].z << endl;
			WARNING_QUIT("md::moveatoms","the movement of atom is larger than the length of cell.");
		}

		cartNoWrap[ia] = taudirac[ia] * ionlatvec*ucell.lat0;
	}
	return;
}

void MD_func::printpos(string file,int iter)
{
//intend to output the positions of atoms to ordered file
	bool pass;
	pass = 0;

	if (dumpmdfreq==1||iter==1||( dumpmdfreq > 1&& iter%dumpmdfreq==0 ) )
		pass =1;
	if (!pass) return;

	string file1=file+".xyz";
	string file2=file+".cif";

	//xiaohui add 'OUT_LEVEL', 2015-09-16
	if(OUT_LEVEL == "i"||OUT_LEVEL == "ie") ucell.print_tau();
	if(OUT_LEVEL == "i"||OUT_LEVEL == "ie") ucell.print_cell_xyz(file1);
	ucell.print_cell_cif(file2);
	stringstream ss;

	ss << global_out_dir << "STRU_MD";

	//zhengdy modify 2015-05-06, outputfile "STRU_Restart"
	ucell.print_stru_file(ss.str(),2);

	return;
}

//rescale velocities to target temperature.
void MD_func::scalevel()
{
	double ke=GetAtomKE();
	if(ke>1e-9)
		for(int i=0;i<numIon;i++){
			vel[i]*=sqrt(3*(numIon-nfrozen)*temperature/ke/2);
		}
	return;
}