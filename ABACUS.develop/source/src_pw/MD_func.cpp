#include "src_pw/MD_func.h"
#include "src_pw/tools.h"
#include "src_lcao/FORCE_STRESS.h"
#include "src_pw/forces.h"
#include "src_pw/stress_pw.h"
#include "src_global/sltk_atom_arrange.h"

bool MD_func::RestartMD(const int& numIon, Vector3<double>* vel, int& step_rst)
{
	int error(0);
	double *vell=new double[numIon*3];
	if (!MY_RANK)
	{
		stringstream ssc;
		ssc << global_readin_dir << "Restart_md.dat";
		ifstream file(ssc.str().c_str());

		if(!file)
		{
			cout<<"please ensure whether 'Restart_md.dat' exists!"<<endl;
			error = 1;
		}
		if (!error)
		{
			file.ignore(11, '\n');

        //----------------------------------------------------------
        // main parameters
        //----------------------------------------------------------
			file.ignore(13, '\n');//read and assert number of atoms
			int num;
			file>>num;
			if(num != numIon)
			{
				cout<<"please ensure whether 'Restart_md.dat' right!"<<endl;
				error = 1;
			}
		}
		if (!error)
		{
			file.get();
			file.ignore(23, '\n');//read velocities
			for(int i = 0;i<numIon*3;i++)
			{
				file>>vell[i];
			}
			file.get();
			file.ignore(6, '\n');//read start step of MD
			file>>step_rst;
			file.close();
		}
	}
#ifdef __MPI
	MPI_Bcast(&error,1,MPI_INT,0,MPI_COMM_WORLD);
#endif
	if(error)
	{
		delete[] vell;
		exit(0);
	}
#ifdef __MPI
	MPI_Bcast(&step_rst,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(vell,numIon*3,MPI_DOUBLE,0,MPI_COMM_WORLD);
#endif
	for(int i=0;i<numIon;i++)
	{
		vel[i].x=vell[i*3];
		vel[i].y=vell[i*3+1];
		vel[i].z=vell[i*3+2];
	}

	delete []vell;
	return true;
}

void MD_func::mdRestartOut(const int& step, const int& recordFreq, const int& numIon, Vector3<double>* vel)
{
//this function used for outputting the information of restart file
	bool pass;
	pass = 0;

	if (recordFreq==1||step==1||( recordFreq > 1&& step%recordFreq==0 ) )
		pass =1;
	if (!pass) return;

	if(!MY_RANK){
		stringstream ssc;
		ssc << global_out_dir << "Restart_md.dat";
		ofstream file(ssc.str().c_str());
		file<<"MD_RESTART"<<endl;
		file<<"ATOM_NUMBERS: "<<numIon<<endl;
		file<<"ION_VELOCITIES_(a.u.): "<<endl;
		for(int i=0;i<numIon;i++){
			file<<setprecision (12)<<vel[i].x<<" "<<setprecision (12)<<vel[i].y<<" "<<setprecision (12)<<vel[i].z<<endl;
		}
		file<<"step: "<<step<<endl;                
		file.close();
	}
	return;
}

double MD_func::GetAtomKE(const int& numIon, const Vector3<double>const* vel, const double const* allmass){
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

void MD_func::InitVelocity(
	const int& numIon, 
	const double& temperature, 
	const double& fundamentalTime, 
	const double const* allmass,
	Vector3<double>* vel)
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

/*void MD_func::ReadNewTemp(int step )
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
}*/

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

int MD_func::getMassMbl(const UnitCell_pseudo &unit_in, double* allmass, Vector3<int>* ionmbl)
{
//some prepared information
//mass and degree of freedom
	int ion=0;
	int nfrozen=0;
	for(int it=0;it<unit_in.ntype;it++){
		for(int i=0;i<unit_in.atoms[it].na;i++){
			allmass[ion]=unit_in.atoms[it].mass/6.02/9.109*1e5;
			ionmbl[ion]=unit_in.atoms[it].mbl[i];
			if (ionmbl[ion].x==0||ionmbl[ion].y==0||ionmbl[ion].z==0) nfrozen++;

			ion++;
		}
	}
	return nfrozen;
}

void MD_func::callInteraction_LCAO(const int& numIon, Vector3<double>* force, matrix& stress_lcao)
{
//to call the force of each atom
	matrix fcs;//temp force matrix
	Force_Stress_LCAO FSL;
	FSL.allocate (); 
	FSL.getForceStress(FORCE, STRESS, TEST_FORCE, TEST_STRESS, fcs, stress_lcao);
	for(int ion=0;ion<numIon;ion++){
		force[ion].x =fcs(ion, 0)/2.0;
		force[ion].y =fcs(ion, 1)/2.0;
		force[ion].z =fcs(ion, 2)/2.0;
	}
#ifdef __MPI //2015-10-01, xiaohui
	atom_arrange::delete_vector( SEARCH_RADIUS );
#endif //2015-10-01, xiaohui

	return;
}
void MD_func::callInteraction_PW(const int& numIon, Vector3<double>* force, matrix& stress_pw)
{
//to call the force of each atom
	matrix fcs;//temp force matrix
	Forces ff;
	ff.init(fcs);
	for(int ion=0;ion<numIon;ion++){
		force[ion].x =fcs(ion, 0)/2.0;
		force[ion].y =fcs(ion, 1)/2.0;
		force[ion].z =fcs(ion, 2)/2.0;
	}
	if(STRESS)
	{
		Stress_PW ss;
		ss.cal_stress(stress_pw);
	}
	return;
}

void MD_func::printpos(const string& file, const int& iter, const int& recordFreq, const UnitCell_pseudo& unit_in)
{
//intend to output the positions of atoms to ordered file
	bool pass;
	pass = 0;

	if (recordFreq==1||iter==1||( recordFreq > 1&& iter%recordFreq==0 ) )
		pass =1;
	if (!pass) return;

	string file1=file+".xyz";
	string file2=file+".cif";

	//xiaohui add 'OUT_LEVEL', 2015-09-16
	if(OUT_LEVEL == "i"||OUT_LEVEL == "ie") unit_in.print_tau();
	if(OUT_LEVEL == "i"||OUT_LEVEL == "ie") unit_in.print_cell_xyz(file1);
	unit_in.print_cell_cif(file2);
	stringstream ss;

	ss << global_out_dir << "STRU_MD";

	//zhengdy modify 2015-05-06, outputfile "STRU_Restart"
	unit_in.print_stru_file(ss.str(),2);

	return;
}

//rescale velocities to target temperature.
void MD_func::scalevel(
	const int& numIon,
	const int& nfrozen,
	const double& temperature,
	Vector3<double>* vel,
	const double* allmass
)
{
	double ke=GetAtomKE(numIon, vel, allmass);
	if(ke>1e-9)
	{
		for(int i=0;i<numIon;i++)
		{
			vel[i]*=sqrt(3*(numIon-nfrozen)*temperature/ke/2);
		}
	}
	return;
}

double MD_func::Conserved(const double KE, const double PE, const int number){
//---------------------------------------------------------------------------
//   This function calculates the conserved quantity for the NVE system. 
//----------------------------------------------------------------------------

   	// KE   // Kinetic energy of particles (due to their velocity)
   	// PE   // Potential energy (DFT total energy)
	// number //number of atoms which have full freedoms

  	double Conserved; // The conserved quantity

   	Conserved = KE + PE ;
   
	if (!MY_RANK)
	{                 
		ofs_running<< "NVE Conservation     : "<< Conserved<<" (Hartree)"<<endl;
		ofs_running<< "NVE Temperature      : "<< 2*KE/(3*number)/K_BOLTZMAN_AU<<" (K)"<<endl;
		ofs_running<< "NVE Kinetic energy   : "<< KE<<" (Hartree)"<<endl;
		ofs_running<< "NVE Potential energy : "<< PE<<" (Hartree)"<<endl;
	}
   	return Conserved;
}

double MD_func::MAXVALF(const int numIon, const Vector3<double>const* force){
	//cout<<"enter in MAXVALF"<<endl;
	double max=0;
	for(int i=0;i<numIon;i++){
		double force0 = pow(force[i].x,2)+pow(force[i].y,2)+pow(force[i].z,2);
		if(max<force0) max = force0;
	}
	return max;
}