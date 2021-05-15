#include "efield.h"
#include "tools.h"
#include "global.h"

Efield::Efield(){}
Efield::~Efield(){}

int Efield::edir; // direction of the field
double Efield::emaxpos; // position of the maximum of the field (0<emaxpos<1)
double Efield::eopreg; // amplitude of the inverse region (0<eopreg<1)
double Efield::eamp; // field amplitude (in a.u.) (1 a.u. = 51.44 10^10 V/m)
double Efield::etotefield; // energy correction due to the field
double Efield::bvec[3];
double Efield::bmod;

void Efield::add_efield(const double*const rho, double* v_in)
{
	TITLE("Efield","add_efield");

	if(!EFIELD) return;

	double avec[3];
	if(edir == 1)
	{
		this->bvec[0] = ucell.G.e11; avec[0] = ucell.latvec.e11;
		this->bvec[1] = ucell.G.e12; avec[1] = ucell.latvec.e12;
		this->bvec[2] = ucell.G.e13; avec[2] = ucell.latvec.e13;
	}
	else if(edir == 2)
	{
		this->bvec[0] = ucell.G.e21; avec[0] = ucell.latvec.e21;
		this->bvec[1] = ucell.G.e22; avec[1] = ucell.latvec.e22;
		this->bvec[2] = ucell.G.e23; avec[2] = ucell.latvec.e23;
	}
	else if(edir == 3)
	{
		this->bvec[0] = ucell.G.e31; avec[0] = ucell.latvec.e31;
		this->bvec[1] = ucell.G.e32; avec[1] = ucell.latvec.e32;
		this->bvec[2] = ucell.G.e33; avec[2] = ucell.latvec.e33;
	}
	else
	{
		throw range_error("Efield::add_efield, edir is < 1 or > 3. "+TO_STRING(__FILE__)+" line "+TO_STRING(__LINE__));
		//WARNING_QUIT("Efield::add_efield","edir is < 1 or > 3.");
	}
//	cout << " bvec=" << bvec[0] << " " << bvec[1] << " " << bvec[2] << endl;

	this->bmod = sqrt(pow(bvec[0],2) + pow(bvec[1],2) + pow(bvec[2],2));

	const double debye = 2.54176;

	double tot_dipole = 0.0;
	double e_dipole = 0.0;
	double ion_dipole = 0.0;

	// e2 is 2.0, means the square of the electron charge
	// calculation of dipole

	matrix fdip(ucell.nat, 3);

	if(DIPOLE)
	{
		// id dipole correction is active
		this->compute_el_dip(emaxpos, eopreg, edir, rho, e_dipole);
		this->compute_ion_dip(emaxpos, eopreg, edir, ion_dipole, this->bmod, bvec);
		tot_dipole = ion_dipole - e_dipole;
		Parallel_Reduce::reduce_double_all(tot_dipole);
		// calculate the correction to total energy
		etotefield = -e2*(eamp-tot_dipole/2.0)*tot_dipole*ucell.omega/FOUR_PI;
		// calculate the force
		if(FORCE)
		{
//			cout << "\n dipole force: " << endl;
			int iat = 0;
			for(int it=0; it<ucell.ntype; ++it)
			{
				for(int ia=0; ia<ucell.atoms[it].na; ++ia)
				{
					for(int jj=0; jj<3; ++jj)
					{
						fdip(iat, jj) = e2 *(eamp - tot_dipole)*ucell.atoms[it].zv*bvec[jj]/bmod;
					}
//					cout << setw(15) << fdip(iat, 0) 
//					<< setw(15) << fdip(iat, 1)
//					<< setw(15) << fdip(iat, 2) << endl;
					++iat;
				}
			}
		}
	}
	else
	{
		// dipole correction is not active
		this->compute_ion_dip(emaxpos, eopreg, edir, ion_dipole, bmod, bvec);
		this->etotefield=-e2*eamp*ion_dipole*ucell.omega/FOUR_PI;
		// calculate the force
		if(FORCE)
		{
			assert(bmod>0);
//			cout << "\n dipole force: " << endl;
			int iat = 0;
			for(int it=0; it<ucell.ntype; ++it)
			{
				for(int ia=0; ia<ucell.atoms[it].na; ++ia)
				{
					for(int jj=0; jj<3; ++jj)
					{
						fdip(iat, jj) = e2 *eamp * ucell.atoms[it].zv * bvec[jj]/bmod;
					}
//					cout << setw(15) << fdip(iat, 0) 
//					<< setw(15) << fdip(iat, 1)
//					<< setw(15) << fdip(iat, 2) << endl;
					++iat;
				}
			}
		}
	}


	// calculation of potentials.
	const double length = (1.0-eopreg)*(ucell.lat0*
		std::sqrt(avec[0]*avec[0]+avec[1]*avec[1]+avec[2]*avec[2]));

	const double vamp = e2*(eamp-tot_dipole)*length;
	//const double fac = ucell.omega/FOUR_PI;

//	cout << " Elec. dipole = " << e_dipole << " (Ry au), also = " 
//	<< e_dipole * debye << " (debye)" << endl; 
	
	ofs_running << " Ion dipole = " << ion_dipole << " (Ry au), also = " 
	<< ion_dipole * debye << " (debye)" << endl;

//	cout << " Dipole = " << (tot_dipole* fac) << " (Ry au), also = " 
//	<< (tot_dipole * fac) * debye << " (debye)" << endl;


	if( abs(eamp) > 0.0) 
	{
		OUT(ofs_running,"Amplitute of Efield (Hartree)",eamp);
		OUT(ofs_running,"Potential amplitute is (Ry)",vamp);
		OUT(ofs_running,"Total length is (Bohr)",length);
	}
	
	//-----------------------------------------------------------
	// add potential
	// V\left(ijk\right) = e^{2} \left( eamp - dip \right) z_{v}
	// Saw\left( \frac{k}{nr3} \right) \frac{alat}{bmod}
	//-----------------------------------------------------------

	int i,j,k,index;	//index0;
	double value;

	if(MY_RANK==0)
	{
		stringstream ss;
		ss << global_out_dir << "EFIELD.dat";
		ofstream ofs(ss.str().c_str());
		
		int npoi;
		if(edir==1) npoi = pw.ncx;
		else if(edir == 2) npoi = pw.ncy;
		else if(edir == 3) npoi = pw.ncz;
		else throw range_error("Efield::add_efield, edir is < 1 or > 3. "+TO_STRING(__FILE__)+" line "+TO_STRING(__LINE__));

		for(int ip=0; ip<npoi; ++ip)
		{
			const double sawarg = (double)ip/(double)npoi;
			value = e2 * (eamp - tot_dipole) *
				saw(emaxpos, eopreg, sawarg) * (ucell.lat0/bmod);
			ofs << ip << " " << value << endl;
		}
		ofs.close();
	}
	
	// compared this to gint_gamma_vl.cpp	
	const int yz = pw.ncy*pw.nczp;
	for(int ir=0; ir<pw.nrxx; ++ir)
	{
		index = ir;
		i     = index / yz; // get the z, z is the fastest
		index = index - yz * i;// get (x,y)
		j     = index / pw.nczp;// get y
		k 	  = index - pw.nczp*j + pw.nczp_start;// get x

		double sawarg;
		if (edir == 1) sawarg = (double)i/(double)pw.ncx;
		else if (edir == 2) sawarg = (double)j/(double)pw.ncy;
		else if (edir == 3) sawarg = (double)k/(double)pw.ncz;
		else throw range_error("Efield::add_efield, edir is < 1 or > 3. "+TO_STRING(__FILE__)+" line "+TO_STRING(__LINE__));

     	value = e2 * (eamp - tot_dipole) * 
			saw(emaxpos, eopreg, sawarg) * (ucell.lat0/bmod);

		/*
		cout << " eamp=" << eamp << endl;
		cout << " tot_dipole=" << tot_dipole << endl;
		cout << " emaxpos=" << emaxpos << endl;
		cout << " eopreg=" << eopreg << endl;
		cout << " sawarg=" << sawarg << endl;
		cout << " bmod=" << bmod << endl;
		*/

     	v_in[ir] += value;
	}

	return;
}

void Efield::compute_el_dip(const double &emaxpos, const double &eopreg,
	const int &edir, const double*const rho, double &e_dipole)const
{
	TITLE("Efield","compute_el_dip");


	return;
}

void Efield::compute_ion_dip(const double &emaxpos, const double &eopreg,
	const int &edir, double &ion_dipole, const double &bmod_in, const double *bvec_in)const
{
	TITLE("Efield","compute_ion_dip");

	// Calculate IONIC dipole
	// 
	// P_{ion} = \sum^{nat}_{s} z_{v} Saw\left( \vec{t_{s}}\cdot\vec{b_{edir}}}
	//                             \right) \frac{alat}{bmod} \frac{4\pi}{\Omega}
	// 

	double tvectb = 0.0;
	ion_dipole = 0.0;
	for(int it=0; it<ucell.ntype; ++it)
	{
		Atom* atom = &ucell.atoms[it];
		for(int ia=0; ia< atom->na; ++ia)
		{
			tvectb = atom->tau[ia].x*bvec_in[0] 
				+ atom->tau[ia].y * bvec_in[1] 
				+ atom->tau[ia].z * bvec_in[2];

			ion_dipole += atom->zv * this->saw(emaxpos,eopreg, tvectb )
			 * (ucell.lat0/bmod_in) * (FOUR_PI/ucell.omega);

//			cout << "--------------------------------------" << endl;
//			cout << " tvectb=" << tvectb << endl;
//			cout << " tau=" << atom->tau[ia].x << " " << atom->tau[ia].y << " " << atom->tau[ia].z << endl;
//			cout << " bvec=" << bvec[0] << " " << bvec[1] << " " << bvec[2] << endl;
//			cout << " saw=" << this->saw(emaxpos,eopreg, tvectb ) << endl;
//			cout << " ion_dipole=" << ion_dipole << endl;
		}	
	}

	return;
}

double Efield::saw(const double &emaxpos, const double &eopreg, const double &x)const
{
	double sawout = 0.0;
	const double z = x - emaxpos;
	const double y = z - std::floor(z);

	if (y<=eopreg)
	{
		sawout = (0.5 - y/eopreg) * (1-eopreg);
	}
	else
	{
		sawout = (-0.5 + (y-eopreg)/(1-eopreg)) * (1-eopreg);
	}
	return sawout;
}

void Efield::compute_force(matrix &fdip)
{
	assert(bmod>0);
	int iat = 0;
	for(int it=0; it<ucell.ntype; ++it)
	{
		for(int ia=0; ia<ucell.atoms[it].na; ++ia)
		{
			for(int jj=0; jj<3; ++jj)
			{
				fdip(iat, jj) = e2 * Efield::eamp * ucell.atoms[it].zv * Efield::bvec[jj]/Efield::bmod;
			}
			++iat;
		}
	}
}
