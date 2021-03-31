//Description:this module is the base of MD simulation
#include<iostream>
#include<string>
#include<fstream>
#include<stdlib.h>
#include<time.h>
#include<cmath>
#include"md.h"
#include"../input.h"
#include"../src_lcao/LOOP_ions.h"
#include "../src_global/sltk_atom_arrange.h" //2015-10-01, xiaohui

md::md(int n)
{
	xLogS=0; // thermostat's position
	vLogS=0 ;// thermostat's speed
	tauthermo = 0;  // The thermostat fluctuation period, for Al fcc at 300K ~ 1000fs
	taubaro   = 0 ; // The barostat   fluctuation period, for Al fcc at 300K ~ 500 fs
	Qmass      = -1;  // Inertia of extended system variable
	dt         = -1; // Time increment (hbar/E_hartree)
	temperature= -1;  //Temperature (in Hartree, 1 Hartree ~ 3E5 K)
	doMSD = 1;   // compute the mean square displacement.
	doMSDAtom = 1;        // compute the mean square displacement for each atom.
	msdstartTime = 1;   //beyond this time, msd is going to be computed
	diffuCoeff    = 0;
	nResn=3;   // Coarseness of Trotter factorization of Nose-Hoover thermostat
	nYosh=3;  // Order of thermostat scheme corresponding with "intCoeff"
	outputstressperiod = 1 ;// The period to output stress
	dumpmdfreq = 1  ;       // The period to dump MD information for monitoring and restarting MD
	fixTemperature = 1;       // The period to read in new temperature during MD run.
	rstMD = 0  ;// 1 : restart MD, vel.restart and ion.restart files will be read
	// 0 : not restart from ion and vel files
	// -n: restart from vel.n and ion.n in directory
	// "md_output_path"

	startStep=0 ;// the first step number of MD
	velRescale = 0 ;// whether to rescale velocities at the first step
	mdoutputpath="mdDataOut";
	step_rst=0;
	numIon=ucell.nat;
	ntype=ucell.ntype;
};

void md::md_allocate()
{
	na=new int[ntype];
	if (!MY_RANK)
	{
		stringstream ssc;
		ssc << global_out_dir << "MD_OUT";
		out.open(ssc.str().c_str());
	}
	for(int i=0;i<ntype;i++)
	{
		na[i]=ucell.atoms[i].na;
	} 
	ionlatvec=ucell.latvec;
	msd=new double[ntype];
	vel=new Vector3<double>[numIon];
	cartNoWrap=new Vector3<double>[numIon];
	taudirac=new Vector3<double>[numIon];
	allmass=new double[numIon];
	ionmbl=new Vector3<int>[numIon];
	force=new Vector3<double>[numIon];
	return;	
}

const double md::fundamentalTime = 2.418884326505e-17;

void md::initMD(void)
{

	mdtype=INPUT.md_mdtype;

	mdenergy=en.etot/2;

	Qmass=INPUT.md_qmass/6.02/9.109*1e5;
	dt=INPUT.md_dt/fundamentalTime/1e15;
	temperature=INPUT.md_tfirst/3.1577464e5;
	rstMD=INPUT.md_rstmd;

	nResn=INPUT.md_nresn;
	nYosh=INPUT.md_nyosh;
	dumpmdfreq=INPUT.md_dumpmdfred;
	fixTemperature=INPUT.md_fixtemperature;
	ediff=INPUT.md_ediff;
	ediffg=INPUT.md_ediffg;
	doMSD=INPUT.md_domsd;
	doMSDAtom=INPUT.md_domsdatom;
	mdoutputpath=INPUT.md_mdoutpath;
	msdstartTime=INPUT.md_msdstartTime;

	if(rstMD==0)
	{
		connection0();
		connection1();
		//A fresh new MD: Do not restart MD
		InitVelocity() ;
		// Initialize thermostat, and barostat
		xLogS = 0.0 ;      // position of thermostat
		vLogS = 0.0;        // velocity of thermostat
	}
	else if( rstMD>0 )
	{
		connection0();
		connection1();
		InitVelocity() ;
		if(!RestartMD()){
			cout<<"error in restart MD!"<<endl;
			return;
		}
		connection2();
	}
	if(Qmass==0)
	{
		Qmass=dt*fundamentalTime*1e15/6.02/9.109*1e5;
	}

	return;
}
	
/*void md::runMD(int step1){
       if(rstMD==1)MDNVT.runnvt();
       else if(mdtype==3)MDNVE.runNVE(step1);
       return;
}*///we need this form

bool md::RestartMD(void)
{
	int i;
	double *vell=new double[numIon*3];
	double *cart=new double[numIon*3];
	double *dira=new double[numIon*3];

	if (!MY_RANK)
	{
		stringstream ssc;
		ssc << global_out_dir << "Restart_md.dat";
		ifstream file(ssc.str().c_str());

		if(!file)
		{
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
		file>>xLogS;
		file.get();	file.ignore(7, '\n');
		file>>vLogS;
		file.get();	file.ignore(23, '\n');
		for(i = 0;i<numIon*3;i++)
		{
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
	MPI_Bcast(&xLogS,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&vLogS,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&step_rst,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(vell,numIon*3,MPI_DOUBLE,0,MPI_COMM_WORLD);
	//      MPI_Bcast(cart,numIon*3,MPI_DOUBLE,0,MPI_COMM_WORLD);
	//      MPI_Bcast(dira,numIon*3,MPI_DOUBLE,0,MPI_COMM_WORLD);
#endif

	for(i=0;i<numIon;i++)
	{
		vel[i].x=vell[i*3];
		vel[i].y=vell[i*3+1];
		vel[i].z=vell[i*3+2];
		//     cartNoWrap[i].x=cart[i*3];
		//     cartNoWrap[i].y=cart[i*3+1];
		//     cartNoWrap[i].z=cart[i*3+2];
		//     taudirac[i].x=dira[i*3];
		//     taudirac[i].y=dira[i*3+1];
		//     taudirac[i].z=dira[i*3+2];
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

	return true;
}


void md::mstout(int step)
{
	//this function used for outputting the information of restart file
	bool pass;
	pass = 0;


	if (dumpmdfreq==1||step==1||( dumpmdfreq > 1&& step%dumpmdfreq==0 ) )
		pass =1;
	if (!pass) return;

	int i;
	if(!MY_RANK)
	{
		stringstream ssc;
		ssc << global_out_dir << "Restart_md.dat";
		ofstream file(ssc.str().c_str());
		file<<"MD_RESTART"<<endl;
		file<<"xLogS: "<<setprecision (12)<<xLogS<<endl;
		file<<"vLogS: "<<setprecision (12)<<vLogS<<endl;
		file<<"ION_VELOCITIES_(a.u.): "<<endl;
		for( i=0;i<numIon;i++)
		{
			file<<setprecision (12)<<vel[i].x
			<<" "<<setprecision (12)
			<<vel[i].y<<" "<<setprecision (12)<<vel[i].z<<endl;
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

double md::GetAtomKE()
{
//---------------------------------------------------------------------------
// DESCRIPTION:
//   This function calculates the classical kinetic energy of a group of atoms.
//----------------------------------------------------------------------------
//!>> EXTERNAL VARIABLES <<!
//vel[3][size(cell%ionTable,1)]  ! Velocities of atoms
// ke ! kinetic energy we want.

//!>> INTERIOR VARIABLES <<!
	double mass;
	int it,ii,ion;

	//!>> INITIALIZATION <<!
	//!>> FUNCTION BODY <<!
	double ke = 0.0;

	//  WRITE(outputUnit,*) "velocities"
	// kinetic energy = \sum_{i} 1/2*m_{i}*(vx^2+vy^2+vz^2)
	for(ion=0;ion<numIon;ion++){
		mass = allmass[ion];
		ke +=   0.5 * mass * (vel[ion].x*vel[ion].x+vel[ion].y*vel[ion].y+vel[ion].z*vel[ion].z);
		//    WRITE(outputUnit,*) vel(1,ion), vel(2,ion), vel(3,ion)
	}
	//cout<<"in GetAtomKE KE="<< ke<<endl;
	return ke;
}


void md::MakeIntCoeff()
{
//---------------------------------------------------------------------------
// DESCRIPTION:
//  Initialize the Yoshida/Suzuki coefficient
//----------------------------------------------------------------------------
	// nYosh ! approximation order
	int i;
	switch (nYosh)
	{
		case 1:
			intCoeff[0] = 1.0;
			break;
		case 3:
			for(i=0;i<3;i++)intCoeff[i] = 1.0/(2.0-pow(2,1/3.0));
			intCoeff[1] = 1.0 - 2.0*intCoeff[1];
			break;
		case 5:
			for(i=0;i<5;i++)intCoeff[i] = 1.0/(4.0-pow(4,1/3.0));
			intCoeff[2] = 1.0 - 4.0*intCoeff[2];
			break;
		default:
			cout<<"ERROR: NYOSH MUST BE 1, 3, or 5.  STOPPING."<<endl;
			exit(0);
	}

	return;
}


//  double md::InitRandomSeed(){
//---------------------------------------------------------------------------
// DESCRIPTION:
//  Initialize the random seed used for random number generator.
//----------------------------------------------------------------------------

  //}//don't know whether it is necessary,so I neglect it.
void md::InitVelocity(void)
{
	if(!MY_RANK)
	{ //xiaohui add 2015-09-25
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
		for(M=0;M<3*numIon;)
		{	
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
			if(y<=(1.0/sqrt(2.0*PI)*exp(-0.5*x*x)))
			{
				if(M<numIon) vel[M].x=x;
				else if(M<2*numIon) vel[M-numIon].y=x;
				else vel[M-2*numIon].z=x;
				// cout<<"this vel is:"<<x<<endl;
				M++;
			}
		}

		int ion;
		for(ion=0;ion<numIon;ion++)
		{
			vel[ion]/=sqrt(allmass[ion]*9.01938e-31);
		}

		for(ion=0;ion<numIon;ion++)
		{
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


/*void md::InitVelocity(){
//----------------------------------------------------------------------------
// DESCRIPTION:
//   This function initializes the velocities of atom before 
//   doing MD simulation.
//----------------------------------------------------------------------------
  
  //! >> EXTERNAL VARIABLES <<!
  //temperature ! temperature to be simulated
  //double v[3][numIon] // velocity of each atom, in a.u.

  //! >> LOCAL VARIABLES << !!
  int i, ii, it,ion,nfrozen;
  double  vel1, vel2, mass,  KE  ;              //KE: to store kinetic energy

  TITLE("md","InitVelocity");

  // note from mohan 2013-03
  // get the random wave functions from processor 1!!
  // otherwise the velocity on each node will be different!
  // then the geometry will be different at each MD step
  // the results are totally wrong.
      srand(   (unsigned)time(   NULL   )   );
    // Assign random numbers to velocities
	  for(i = 1; i<numIon*3/2 + 2;i++){
      
      vel1 = sqrt(-2*log((rand()%10000)/10000.0))*cos(2*PI*(rand()%10000)/10000.0);
      vel2 = sqrt(-2*log((rand()%10000)/10000.0))*sin(2*PI*(rand()%10000)/10000.0);
      if ((2*i-2)/3 < numIon){
		  if(2*i-2%3==0)vel[(2*i-2)/3].x = vel1;
		  else if((2*i-2)%3==1)vel[(2*i-2)/3].y = vel1;
		  else if((2*i-2)%3==2)vel[(2*i-2)/3].z = vel1;
	  }
	  if ((2*i-1)/3 < numIon){ 
		  if(2*i-1%3==0)vel[(2*i-1)/3].x = vel2;
		  else if((2*i-1)%3==1)vel[(2*i-1)/3].y = vel2;
		  else if((2*i-1)%3==2)vel[(2*i-1)/3].z = vel2;
	  }

    // Scale velocity according to mass of atoms
	  for(ion = 0;ion<numIon;ion++){
              mass = allmass[ion];
	          vel[ion] = vel[ion]/mass;			  
	  }
	  
  
    RemoveMovementOfCenterOfMass();
  
    // Do not move atoms which are frozen
    nfrozen = 0;
	for(ion=0;ion<numIon;ion++){
         if(ionmbl[ion].x==0)vel[ion].x=0;
         if(ionmbl[ion].y==0)vel[ion].y=0;
         if(ionmbl[ion].z==0)vel[ion].z=0;
		 nfrozen++;
	}

  // Rescale all velocities to start to exact temperature desired
    KE=GetAtomKE();
    cout<<"Kinetic energy from random distribution of velocities = "<< KE;
    for(ion=0;ion<numIon;ion++){
       vel[ion] = vel[ion] * sqrt(temperature*1.5 * double( numIon-nfrozen) / KE);
	}
  //mohan add 2013-04-17
    KE= GetAtomKE();
    out<<"In order to start from the input temperature, kinetic energy should be = "<<KE<<endl;

  // mohan fix bug 2013-03
  // distribute the random velocity from first processor to other processors

//  WRITE(outputUnit,*) " init velocity:"
// DO i=1, numIon
//   WRITE(outputUnit,*) v(1,i), v(2,i), v(3,i)
// ENDDO
    
  return ;
  }
}
*/


void md::RemoveMovementOfCenterOfMass(void)
{
//---------------------------------------------------------------------------
// DESCRIPTION:
//    Compute the velocity for the center of mass, and remove it
//----------------------------------------------------------------------------
  //double vel[3][numIon] // velocties of each atom.

  //!! >> LOCAL VARIABLES << !!
  int ii,it, id,j,ion;
  Vector3<double> avgvel;
  double mass,totMass;

  //!! >> FUNCTION << !!
  
  avgvel.x = 0;
  avgvel.y=0;
  avgvel.z=0;
  totMass = 0;

  // sum up the total momentum of all the atoms
  // in the system. (Px,Py,Pz) =\sum m_{i}*(vx,vy,vz)_{i}
  for(ion=0;ion<numIon;ion++){
        if(ionmbl[ion].x==0||ionmbl[ion].y==0||ionmbl[ion].z==0)
        {
        continue;
        }

        mass = allmass[ion];
        avgvel += mass * vel[ion];
        totMass = totMass + mass;
  }
	//cout<<"avgvel: "<<avgvel<<endl;
    //    cout<<"numIon:"<<numIon<<endl;
    //    cout<<"last mass:"<<mass<<endl;
    //    cout<<"last vel:"<<vel[ion].x<<" "<<vel[ion].y<<endl;
//	cout<<"totMass: "<<totMass<<endl;
    avgvel = avgvel / totMass;
//  cout<<setw(15)<<" Average velocity : "<<avgvel.x<<" "<< avgvel.y<<" "<< avgvel.z<<endl;
  for(ion=0;ion<numIon;ion++){
    if(ionmbl[ion].x==0||ionmbl[ion].y==0||ionmbl[ion].z==0) 
    {
        continue;
    }
    vel[ion] = vel[ion] - avgvel;
  }
  return ;
}


  void md::OutputMDHeader(){
//---------------------------------------------------
// Output the Header for MD simulation
//---------------------------------------------------

  
  
  //!!>>> FUNCTION BODY <<<!!

  // Set up extra output to ion optimizer / MD header
  if (!MY_RANK)
  out<<"Q ="<<Qmass<<"a.u.\n"<<" dt = "<<dt*fundamentalTime*1e15<<"fs\n"<< "T ="<<temperature/K_BOLTZMAN_AU<<endl;
  ofs_running<<"Q ="<<Qmass<<"a.u.\n"<<" dt = "<<dt*fundamentalTime*1e15<<"fs\n"<< "T ="<<temperature/K_BOLTZMAN_AU<<endl;
  return;
 } 



  void md::MonitorMeanSquareDisplacement(int step){
//--------------------------------------------------------------------------
// DESCRIPTION:
//  Calculate the Mean Square Displcacement of particles
//  Compute the mean square displacement of atoms according to formula:
//    MSD = <|r(t)-r(0)|^2>
//  <...> denotes the averaging over all the atoms
//
//  For solids, MSD saturates to a finte value, for liquids, 
//  MSD grows linearly with time.
// 
//  Also compute the diffusion coefficient D:  
//   D = MSD/(6*t)
//----------------------------------------------------------------------------

  //! >> LOCAL VARIABLES << !
        Matrix3 box;
        int it,ii,jj,ion,i; 
        double timePassed;
        double *tmp;
        tmp =new double[ntype];
        Vector3<double> diff; 
  // INITIALIZATION
   //       box.e11 =ionlatvec.e11;
  //	  box.e12 =ionlatvec.e12;
//	  box.e13=ionlatvec.e13;
//	  box.e21 =ionlatvec.e21;
//	  box.e22 =ionlatvec.e22;
//	  box.e23 =ionlatvec.e23;
//	  box.e31 =ionlatvec.e31;
//	  box.e32 =ionlatvec.e32;
//	  box.e33 =ionlatvec.e33;
        box=ionlatvec;



  //! >> FUNCTION << !
        if(msdstartTime>step-step_rst)return;
        if (doMSD ==1) {
    
               if (step-step_rst == msdstartTime) {
                    msdcoords0=new Vector3<double>[numIon];
                    msdcoordst=new Vector3<double>[numIon];
                    if(rstMD==1){
                         double *msd0=new double[numIon*3];
                         double *msdt=new double[numIon*3];
                         if (!MY_RANK){
                              stringstream ssc;
                              ssc << global_out_dir << "Restart_msd.dat";
                              ifstream file(ssc.str().c_str());
                              if(!file){
                                   cout<<"please ensure whether 'Restart_msd.dat' exists"<<endl;
                                   rstMD=0;
                              }
                              if(rstMD==1){
                                   for(i=0;i<numIon*3;i++){
                                        file>>msd0[i];
                                   }
                                   for(i=0;i<numIon*3;i++){
                                        file>>msdt[i];
                                   }
                              }
                              file.close();
                         }
//2015-09-05, xiaohui
#ifdef __MPI
                         MPI_Bcast(msd0,numIon*3,MPI_DOUBLE,0,MPI_COMM_WORLD);
                         MPI_Bcast(msdt,numIon*3,MPI_DOUBLE,0,MPI_COMM_WORLD);
#endif
                         for(i=0;i<numIon;i++){
                              msdcoords0[i].x=msd0[i*3];
                              msdcoords0[i].y=msd0[i*3+1];
                              msdcoords0[i].z=msd0[i*3+2];
                              msdcoordst[i].x=msdt[i*3];
                              msdcoordst[i].y=msdt[i*3+1];
                              msdcoordst[i].z=msdt[i*3+2];
                         } 
                   }
                   if(rstMD==0){for(ion=0;ion<numIon;ion++){
                         msdcoords0[ion] = cartNoWrap[ion];
                         msdcoordst[ion] = cartNoWrap[ion];
		   }
                   }
	      }

	      if (step-step_rst > msdstartTime){
		   for(ion=0;ion<numIon;ion++){
                         msdcoordst[ion] += vel[ion]*dt;
		   }
	      }
      //-----------------------------------
      // compute msd and diffuCoeff 
              for(ii=0;ii<ntype;ii++){ 
                   msd[ii] = 0;
                   tmp[ii]=0;
              }

              ion=0;
              for( ii=0;ii<ntype;ii++){
                   for(jj=0;jj<na[ii];jj++){// loop over ions
                       diff = msdcoords0[ion] - msdcoordst[ion];
                       // the mean square displacement has unit of Bohr^2,
                       // so we should not use 'SQRT' here.
                       // mohan note 2013-05-01
                       //tmp = SQRT( SUM(diff**2) )
                       tmp[ii] += diff*diff;
                       ion++;
                   }
                   msd[ii] = msd[ii] + tmp[ii];
              }
              if (!MY_RANK){
                   for(ii=0;ii<ntype;ii++){ 
                        out<<"msd_single_atom:"<< ii<<" "<< tmp[ii]<<endl;
			//xiaohui remove, 2015-09-16
                        //ofs_running<<"msd_single_atom:"<< ii<<" "<< tmp[ii]<<endl;

                   }
                   if(rstMD==0)
                        timePassed = dt * ( step - step_rst - msdstartTime );
                   else timePassed = dt * ( step - msdstartTime );
                   for(ii=0;ii<ntype;ii++){
                        msd[ii] = msd[ii] / double(na[ii]);
                        diffuCoeff = msd[ii] / (6.0*timePassed);

			//xiaohui move 2 lines, add 2 another lines 2015-09-25, setprecision
                        //out<<ii<< "msd:STEP="<<step<< " msd= "<<msd[ii]<<" (Bohr^2), D= "<<diffuCoeff<< " a.u."<<endl;
                        //ofs_running<<" "<<ii<< "msd:STEP="<<step<< " msd= "<<msd[ii]<<" (Bohr^2), D= "<<diffuCoeff<< " a.u."<<endl;
                        out<<ii<< "msd:STEP="<<step<< " msd= "<<setprecision (12)<<msd[ii]<<" (Bohr^2), D= "<<setprecision (12)<<diffuCoeff<< " a.u."<<endl;
                        ofs_running<<ii<< "msd:STEP="<<step<< " msd= "<<setprecision (12)<<msd[ii]<<" (Bohr^2), D= "<<setprecision (12)<<diffuCoeff<< " a.u."<<endl;
                   }
              }

// NOTE !
// to convert the diffuCoeff to A^2/PIcosecond, 
// D * bohr * bohr / ( fundamentalTime / 1.0e-12 )
// where fundamentalTime is 2.418884326505E-17_DP
//!DEBUG          
//          DO ii = 1,numIon
//            WRITE(message,'(a,I6,a,ES12.4)') & 
//              ' ion:',ii,' ion_sd:', & 
//              SQRT(SUM((msd_coordst(:,ii)-msd_coords0(:,ii))**2))
//            CALL WRTOUT(6,Message)
//            WRITE(message,'(a,I6,a,3ES12.4)') & 
//              ' ion:',ii,' ion_fracCoord:', cell%ionTable(ii)%coord(:)
//              CALL WRTOUT(6,Message)
//          ENDDO ! loop over ions
//!~DEBUG

	 // if (step > msd_startTime)
               if(step==NSTEP){
                     if(msdcoords0!=NULL)delete []msdcoords0;
                     if(msdcoordst!=NULL)delete []msdcoordst;
               }
	  } // doMSD
     return;
  }


  void md::OutputIonVelocityAndStats(int iter){
//---------------------------------------------------------------------------
// DESCRIPTION:
//  Output the 
//   (1) ion's velocity 
//   (2) thermostat's xLogS, and vLogs (for NVT, NPT)
//  (3) barostat's vbox (box's velocity)  (NPT only)
// for restart purpose
//---------------------------------------------------------------------------
// REVISION LOG:
//    Dec/29/2008: Created by Chen Huang
//----------------------------------------------------------------------------
  
  // iter, mdtype //number of ions, iteration number, NVT/NPT

//2015/4/30:change to output msd data for restart.
  int i,j;

  bool pass;
  pass = 0;


  if (dumpmdfreq==1||iter==1||( dumpmdfreq > 1&& iter%dumpmdfreq==0 ) )
    pass =1;
  if (!pass) return;
  
  if (msdstartTime>iter)return;
  if (!MY_RANK){
        stringstream ssc;
        ssc << global_out_dir << "Restart_msd.dat";
        ofstream vfile(ssc.str().c_str());

/*  vfile<<"MD STEP= "<<iter<<endl;
  //----- thermo-stat's velocity and positions
  vfile<<" xLogS: "<< xLogS<< " <= thermostat position :(a.u.)"<<endl;
  vfile<<" vLogS: "<< vLogS<< " <= thermostat velocity :(a.u.)"<<endl;
	  switch(mdtype){
  case 1: // NVT
    vfile<<" 0.0 0.0 0.0 <= barostat velocity 1"<<endl;
    vfile<< " 0.0 0.0 0.0 <= barostat velocity 2"<<endl;
    vfile<< " 0.0 0.0 0.0 <= barostat velocity 3"<<endl;    // add addiontional outputs if needed
    break;  
  case 2: // NPT
    vfile<<vBoxG.e11<<" "<<vBoxG.e12<<" "<<vBoxG.e13<< " <= barostat velocity : vBoxG[0][]"<<endl;
    vfile<<vBoxG.e21<<" "<<vBoxG.e22<<" "<<vBoxG.e23<< " <= barostat velocity : vBoxG[1][]"<<endl; 
    vfile<<vBoxG.e31<<" "<<vBoxG.e32<<" "<<vBoxG.e33<< " <= barostat velocity : vBoxG[2][]"<<endl; 
    break; 
  case 3: // NVE
    vfile<< " 0.0 0.0 0.0 <= barostat velocity 1"<<endl;
    vfile<< " 0.0 0.0 0.0 <= barostat velocity 2"<<endl;
    vfile<< " 0.0 0.0 0.0 <= barostat velocity 3"<<endl;
	break;
  default: break;
  }
  // output velocities
  vfile<<"ION VELOCITIES (a.u.):"<<endl;
	  for( i = 0;i<numIon;i++){
          vfile<<"  "<< vel[i].x<<"  "<<vel[i].y<<"  "<<vel[i].z;
	  }
*/
   //vfile<<"msd0: " <<endl;
        for(i=0;i<numIon;i++){
                vfile<<setprecision (12)<<msdcoords0[i].x<<" "<<setprecision (12)<<msdcoords0[i].y<<" "<<setprecision (12)<<msdcoords0[i].z<<endl;
        }
     //   vfile<<"msdt: " <<endl;
        for(i=0;i<numIon;i++){
                vfile<<setprecision (12)<<msdcoordst[i].x<<" "<<setprecision (12)<<msdcoordst[i].y<<" "<<setprecision (12)<<msdcoordst[i].z<<endl;
        }

  vfile.close();
  }
   return;
}



void md::OutputMDGeom(int iter){
//---------------------------------------------------------------------------
// DESCRIPTION:
//  Output the ion's position and cell's structure for monitoring and MD
//  restart.
//---------------------------------------------------------------------------
// REVISION LOG:
//    Dec/29/2008: Created by Chen Huang
//----------------------------------------------------------------------------
  
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




void md::ReadNewTemp(int step ){
//-----------------------------------------------------------------
//  If fixTemperature == 0, then this subroutine will be skipped
//  otherwise, we are going to read a new temperature from the disk.
//  You only need to create a file named temperature.dat, and then 
//  change the temperature in it. 
//-----------------------------------------------------------------
 

    double intemp;
    int sstep=0;

    if ( fixTemperature > 0 ){

    // Change the temperature every 'fixTemperature' steps.
	  if( (step!=1)&&(step%fixTemperature == 1) ){
    
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
string md::intTurnTostring(int iter,string path){
//int number turn to string term
	int i[6],k;
	if(iter>99999){
        k=6;
        i[5]=iter%10;
		iter=iter/10;
		i[4]=iter%10;
		iter=iter/10;
		i[3]=iter%10;
		iter=iter/10;
		i[2]=iter%10;
		iter=iter/10;
		i[1]=iter%10;
		iter=iter/10;
		i[0]=iter%10;
	}
	else if(iter>9999){
		k=5;
		i[4]=iter%10;
		iter=iter/10;
		i[3]=iter%10;
		iter=iter/10;
		i[2]=iter%10;
		iter=iter/10;
		i[1]=iter%10;
		iter=iter/10;
		i[0]=iter%10;
	}
	else if(iter>999){
		k=4;
		i[3]=iter%10;
		iter=iter/10;
		i[2]=iter%10;
		iter=iter/10;
		i[1]=iter%10;
		iter=iter/10;
		i[0]=iter%10;
	}
	else if(iter>99){
		k=3;
		i[2]=iter%10;
		iter=iter/10;
		i[1]=iter%10;
		iter=iter/10;
		i[0]=iter%10;
	}
	else if(iter>9){
		k=2;
		i[1]=iter%10;
		iter=iter/10;
		i[0]=iter%10;
	}
	else if(iter>0){
		k=1;
		i[0]=iter;
	}
	for(int j=0;j<k;j++){
		path+=(i[j]+48);
	}
	return path;
}
void md::connection0(){
//some prepared imformation
    int i,it,ion;
	ion=0;
        nfrozen=0;
	for(it=0;it<ntype;it++){
		for(i=0;i<na[it];i++){
                        allmass[ion]=ucell.atoms[it].mass/6.02/9.109*1e5;
	         	ionmbl[ion]=ucell.atoms[it].mbl[i];
                        if (ionmbl[ion].x==0||ionmbl[ion].y==0||ionmbl[ion].z==0) nfrozen++;

			ion++;
		}
	}
	return;
}
	
void md::connection1(){
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
void md::connection2(){
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
void md::callforce(){
//to call the force of each atom
	int ion;
	Force_Stress_LCAO FSL;
	FSL.allocate (); 
	matrix forceMd;
	FSL.getForceStress(FORCE, STRESS, TEST_FORCE, TEST_STRESS, forceMd, stress_lcao);
	for(ion=0;ion<numIon;ion++){
			force[ion].x=forceMd(ion,0)/2.0;
			force[ion].y=forceMd(ion,1)/2.0;
			force[ion].z=forceMd(ion,2)/2.0;
                       // cout<<force[ion].x/0.0195<<" "<<force[ion].y/0.0195<<" "<<force[ion].z/0.0195<<endl;
	}
#ifdef __MPI //2015-10-01, xiaohui
        atom_arrange::delete_vector( SEARCH_RADIUS );
#endif //2015-10-01, xiaohui
	return;
}
void md::moveatoms(int step){
// intent to move atoms to the center of the cell
//intent to ensure that no atom is out of the cell
      int ia;
      int ii,it, id,j,ion;
      Vector3<double> avgdc;
      double mass,totMass;

  //!! >> FUNCTION << !!
 
// 2015/9/18 no move of center of mass
/* 
        avgdc.x = 0;
        avgdc.y=0;
        avgdc.z=0;
        totMass = 0;

  // sum up the total momentum of all the atoms
  // in the system. (Px,Py,Pz) =\sum m_{i}*(vx,vy,vz)_{i}
        for(ion=0;ion<numIon;ion++){
                         mass = allmass[ion];
                         avgdc += mass * taudirac[ion];
                         totMass = totMass + mass;
        }
	//cout<<"avgdc: "<<avgdc<<endl;
        avgdc = avgdc / totMass;
        avgdc.x += -0.5;
        avgdc.y += -0.5;
        avgdc.z += -0.5;
       // cout<<" Average dirac coords : "<<avgdc.x<<" "<< avgdc.y<<" "<< avgdc.z<<endl;
        for(ion=0;ion<numIon;ion++){ 
                        taudirac[ion] = taudirac[ion] - avgdc;
        }
      }*/
      for(ia=0;ia<numIon;ia++){
//changed by a constrain box
                 
       if(taudirac[ia].x<0) taudirac[ia].x += 1.0;
                        if(taudirac[ia].y<0) taudirac[ia].y += 1.0;
                        if(taudirac[ia].z<0) taudirac[ia].z += 1.0;
                        if(taudirac[ia].x>=1.0) taudirac[ia].x -= 1.0;
                        if(taudirac[ia].y>=1.0) taudirac[ia].y -= 1.0;
                        if(taudirac[ia].z>=1.0) taudirac[ia].z -= 1.0;
/*
                        if(taudirac[ia].x<0){
                                    taudirac[ia].x *= -1.0;
                                    vel[ia].x *=-1.0;
                        }
                        if(taudirac[ia].y<0) {
                                    taudirac[ia].y *= -1.0;
                                    vel[ia].y *=-1.0;
                        }
                        if(taudirac[ia].z<0) {
                                    taudirac[ia].z *= -1.0;
                                    vel[ia].z *=-1.0;
                        }
                        if(taudirac[ia].x>=1.0){
                                    taudirac[ia].x = 2.0-taudirac[ia].x;
                                    vel[ia].x *=-1.0;
                        }
                        if(taudirac[ia].y>=1.0){
                                    taudirac[ia].y = 2.0-taudirac[ia].y;
                                    vel[ia].y *=-1.0;
                        }
                        if(taudirac[ia].z>=1.0){
                                    taudirac[ia].z = 2.0-taudirac[ia].z;
                                    vel[ia].z *=-1.0;
                        }
*/
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
void md::printpos(string file,int iter){
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
void md::md_release(){
          
          if(force!=NULL)delete []force;
          if(ionmbl!=NULL)delete []ionmbl;
          if(cartNoWrap!=NULL)delete []cartNoWrap;
          if(vel!=NULL)delete []vel;
          if(taudirac!=NULL)delete []taudirac;
          if(allmass!=NULL)delete []allmass;
          
          return;
}
void md::scalevel(){
   int i;
   double ke=GetAtomKE();
   if(ke>1e-9)
      for(i=0;i<numIon;i++){
          vel[i]*=sqrt(3*(numIon-nfrozen)*temperature/ke/2);
      }
   return;
}
/*void md::PDF(int step){ 
   int i,j,k;
   double dis,c;
   c=(ionlatvec.e11+ionlatvec.e12+ionlatvec.e13)*ucell.lat0;
   if(step==1)for(i=0;i<100;i++){
            rdf[i]=0;
       }
   for(i=0;i<numIon;i++){
       for(j=i;j<numIon;j++){
          for(k=0;k<100;k++){
              dis=sqrt((cartNoWrap[i]-cartNoWrap[j])*(cartNoWrap[i]-cartNoWrap[j]));
              if(dis<((double)(k+1)*c/(double)100)&&dis>((double)k*c/(double)100))rdf[k]=rdf[k]+1;
          }
       }
   }
   return;
}
void md::printRDF(int step){
   int i;
   ofstream file;
//   if(step==1){
//        file.open("RDFnow.dat");
 //       file.close();
  // }
   string it="RDF";
   if(fixTemperature!=0)
        if(step/fixTemperature!=0) intTurnTostring(step/fixTemperature,it);
   it+=".dat";
   if(step%dumpmdfreq==1){
        file.open(it.c_str());
        file<<step<<endl;
        for(i=0;i<60;i++){
           file<<"number "<<i<<" RDF: "<<(double)rdf[i]<<endl;
        }
        file.close();
   }
   return;
}
void md::RDF(){
     int i,j,k;
     double dis;
     for(i=0;i<60;i++){
            rdf[i]=0;
     }
     for(i=0;i<numIon;i++){
        if(cartNoWrap[i].x-11.342>0&&cartNoWrap[i].x+11.342<0)
            for(j=0;j<numIon;j++){
                for(k=0;k<60;k++){
                     dis=sqrt((cartNoWrap[i]-cartNoWrap[j])*(cartNoWrap[i]-cartNoWrap[j]));
                     if(dis<((double)k/5.29+0.2835)&&dis>((double)k/5.29+0.0945))rdf[k]=rdf[k]+1;
                }
            }
     }
     return;
}*/
