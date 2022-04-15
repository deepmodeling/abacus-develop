#include "efield.h"
#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "global.h"
#include "../src_parallel/parallel_reduce.h"

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
	ModuleBase::TITLE("Efield","add_efield");

	if(!GlobalV::EFIELD) return;

	double avec[3];
	if(edir == 1)
	{
		this->bvec[0] = GlobalC::ucell.G.e11; avec[0] = GlobalC::ucell.latvec.e11;
		this->bvec[1] = GlobalC::ucell.G.e12; avec[1] = GlobalC::ucell.latvec.e12;
		this->bvec[2] = GlobalC::ucell.G.e13; avec[2] = GlobalC::ucell.latvec.e13;
	}
	else if(edir == 2)
	{
		this->bvec[0] = GlobalC::ucell.G.e21; avec[0] = GlobalC::ucell.latvec.e21;
		this->bvec[1] = GlobalC::ucell.G.e22; avec[1] = GlobalC::ucell.latvec.e22;
		this->bvec[2] = GlobalC::ucell.G.e23; avec[2] = GlobalC::ucell.latvec.e23;
	}
	else if(edir == 3)
	{
		this->bvec[0] = GlobalC::ucell.G.e31; avec[0] = GlobalC::ucell.latvec.e31;
		this->bvec[1] = GlobalC::ucell.G.e32; avec[1] = GlobalC::ucell.latvec.e32;
		this->bvec[2] = GlobalC::ucell.G.e33; avec[2] = GlobalC::ucell.latvec.e33;
	}
	else
	{
		throw range_error("Efield::add_efield, edir is < 1 or > 3. "+ModuleBase::GlobalFunc::TO_STRING(__FILE__)+" line "+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
		//ModuleBase::WARNING_QUIT("Efield::add_efield","edir is < 1 or > 3.");
	}
//	std::cout << " bvec=" << bvec[0] << " " << bvec[1] << " " << bvec[2] << std::endl;

	this->bmod = sqrt(pow(bvec[0],2) + pow(bvec[1],2) + pow(bvec[2],2));

	const double debye = 2.54176;

	double tot_dipole = 0.0;
	double e_dipole = 0.0;
	double ion_dipole = 0.0;

	// ModuleBase::e2 is 2.0, means the square of the electron charge
	// calculation of dipole

	ModuleBase::matrix fdip(GlobalC::ucell.nat, 3);

	if(GlobalV::DIPOLE)
	{
		// id dipole correction is active
		this->compute_el_dip(emaxpos, eopreg, edir, rho, e_dipole);
		this->compute_ion_dip(emaxpos, eopreg, edir, ion_dipole, this->bmod, bvec);
		tot_dipole = ion_dipole - e_dipole;
		Parallel_Reduce::reduce_double_all(tot_dipole);
		// calculate the correction to total energy
		etotefield = -ModuleBase::e2*(eamp-tot_dipole/2.0)*tot_dipole*GlobalC::ucell.omega/ModuleBase::FOUR_PI;
		// calculate the force
		if(GlobalV::CAL_FORCE)
		{
//			std::cout << "\n dipole force: " << std::endl;
			int iat = 0;
			for(int it=0; it<GlobalC::ucell.ntype; ++it)
			{
				for(int ia=0; ia<GlobalC::ucell.atoms[it].na; ++ia)
				{
					for(int jj=0; jj<3; ++jj)
					{
						fdip(iat, jj) = ModuleBase::e2 *(eamp - tot_dipole)*GlobalC::ucell.atoms[it].zv*bvec[jj]/bmod;
					}
//					std::cout << std::setw(15) << fdip(iat, 0) 
//					<< std::setw(15) << fdip(iat, 1)
//					<< std::setw(15) << fdip(iat, 2) << std::endl;
					++iat;
				}
			}
		}
	}
	else
	{
		// dipole correction is not active
		this->compute_ion_dip(emaxpos, eopreg, edir, ion_dipole, bmod, bvec);
		this->etotefield=-ModuleBase::e2*eamp*ion_dipole*GlobalC::ucell.omega/ModuleBase::FOUR_PI;
		// calculate the force
		if(GlobalV::CAL_FORCE)
		{
			assert(bmod>0);
//			std::cout << "\n dipole force: " << std::endl;
			int iat = 0;
			for(int it=0; it<GlobalC::ucell.ntype; ++it)
			{
				for(int ia=0; ia<GlobalC::ucell.atoms[it].na; ++ia)
				{
					for(int jj=0; jj<3; ++jj)
					{
						fdip(iat, jj) = ModuleBase::e2 *eamp * GlobalC::ucell.atoms[it].zv * bvec[jj]/bmod;
					}
//					std::cout << std::setw(15) << fdip(iat, 0) 
//					<< std::setw(15) << fdip(iat, 1)
//					<< std::setw(15) << fdip(iat, 2) << std::endl;
					++iat;
				}
			}
		}
	}


	// calculation of potentials.
	const double length = (1.0-eopreg)*(GlobalC::ucell.lat0*
		std::sqrt(avec[0]*avec[0]+avec[1]*avec[1]+avec[2]*avec[2]));

	const double vamp = ModuleBase::e2*(eamp-tot_dipole)*length;
	//const double fac = GlobalC::ucell.omega/ModuleBase::FOUR_PI;

//	std::cout << " Elec. dipole = " << e_dipole << " (Ry au), also = " 
//	<< e_dipole * debye << " (debye)" << std::endl; 
	
	GlobalV::ofs_running << " Ion dipole = " << ion_dipole << " (Ry au), also = " 
	<< ion_dipole * debye << " (debye)" << std::endl;

//	std::cout << " Dipole = " << (tot_dipole* fac) << " (Ry au), also = " 
//	<< (tot_dipole * fac) * debye << " (debye)" << std::endl;


	if( abs(eamp) > 0.0) 
	{
		ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"Amplitute of Efield (Hartree)",eamp);
		ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"Potential amplitute is (Ry)",vamp);
		ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"Total length is (Bohr)",length);
	}
	
	//-----------------------------------------------------------
	// add potential
	// V\left(ijk\right) = e^{2} \left( eamp - dip \right) z_{v}
	// Saw\left( \frac{k}{nr3} \right) \frac{alat}{bmod}
	//-----------------------------------------------------------

	int i,j,k,index;	//index0;
	double value;

	if(GlobalV::MY_RANK==0)
	{
		std::stringstream ss;
		ss << GlobalV::global_out_dir << "EFIELD.dat";
		std::ofstream ofs(ss.str().c_str());
		
		int npoi;
		if(edir==1) npoi = GlobalC::pw.ncx;
		else if(edir == 2) npoi = GlobalC::pw.ncy;
		else if(edir == 3) npoi = GlobalC::pw.ncz;
		else throw range_error("Efield::add_efield, edir is < 1 or > 3. "+ModuleBase::GlobalFunc::TO_STRING(__FILE__)+" line "+ModuleBase::GlobalFunc::TO_STRING(__LINE__));

		for(int ip=0; ip<npoi; ++ip)
		{
			const double sawarg = (double)ip/(double)npoi;
			value = ModuleBase::e2 * (eamp - tot_dipole) *
				saw(emaxpos, eopreg, sawarg) * (GlobalC::ucell.lat0/bmod);
			ofs << ip << " " << value << std::endl;
		}
		ofs.close();
	}
	
	// compared this to gint_gamma_vl.cpp	
	const int yz = GlobalC::pw.ncy*GlobalC::pw.nczp;
	for(int ir=0; ir<GlobalC::pw.nrxx; ++ir)
	{
		index = ir;
		i     = index / yz; // get the z, z is the fastest
		index = index - yz * i;// get (x,y)
		j     = index / GlobalC::pw.nczp;// get y
		k 	  = index - GlobalC::pw.nczp*j + GlobalC::pw.nczp_start;// get x

		double sawarg;
		if (edir == 1) sawarg = (double)i/(double)GlobalC::pw.ncx;
		else if (edir == 2) sawarg = (double)j/(double)GlobalC::pw.ncy;
		else if (edir == 3) sawarg = (double)k/(double)GlobalC::pw.ncz;
		else throw range_error("Efield::add_efield, edir is < 1 or > 3. "+ModuleBase::GlobalFunc::TO_STRING(__FILE__)+" line "+ModuleBase::GlobalFunc::TO_STRING(__LINE__));

     	value = ModuleBase::e2 * (eamp - tot_dipole) * 
			saw(emaxpos, eopreg, sawarg) * (GlobalC::ucell.lat0/bmod);

		/*
		std::cout << " eamp=" << eamp << std::endl;
		std::cout << " tot_dipole=" << tot_dipole << std::endl;
		std::cout << " emaxpos=" << emaxpos << std::endl;
		std::cout << " eopreg=" << eopreg << std::endl;
		std::cout << " sawarg=" << sawarg << std::endl;
		std::cout << " bmod=" << bmod << std::endl;
		*/

     	v_in[ir] += value;
	}

	return;
}

void Efield::compute_el_dip(const double &emaxpos, const double &eopreg,
	const int &edir, const double*const rho, double &e_dipole)const
{
	ModuleBase::TITLE("Efield","compute_el_dip");


	return;
}

void Efield::compute_ion_dip(const double &emaxpos, const double &eopreg,
	const int &edir, double &ion_dipole, const double &bmod_in, const double *bvec_in)const
{
	ModuleBase::TITLE("Efield","compute_ion_dip");

	// Calculate IONIC dipole
	// 
	// P_{ion} = \sum^{nat}_{s} z_{v} Saw\left( \vec{t_{s}}\cdot\vec{b_{edir}}}
	//                             \right) \frac{alat}{bmod} \frac{4\pi}{\Omega}
	// 

	double tvectb = 0.0;
	ion_dipole = 0.0;
	for(int it=0; it<GlobalC::ucell.ntype; ++it)
	{
		Atom* atom = &GlobalC::ucell.atoms[it];
		for(int ia=0; ia< atom->na; ++ia)
		{
			tvectb = atom->tau[ia].x*bvec_in[0] 
				+ atom->tau[ia].y * bvec_in[1] 
				+ atom->tau[ia].z * bvec_in[2];

			ion_dipole += atom->zv * this->saw(emaxpos,eopreg, tvectb )
			 * (GlobalC::ucell.lat0/bmod_in) * (ModuleBase::FOUR_PI/GlobalC::ucell.omega);

//			std::cout << "--------------------------------------" << std::endl;
//			std::cout << " tvectb=" << tvectb << std::endl;
//			std::cout << " tau=" << atom->tau[ia].x << " " << atom->tau[ia].y << " " << atom->tau[ia].z << std::endl;
//			std::cout << " bvec=" << bvec[0] << " " << bvec[1] << " " << bvec[2] << std::endl;
//			std::cout << " saw=" << this->saw(emaxpos,eopreg, tvectb ) << std::endl;
//			std::cout << " ion_dipole=" << ion_dipole << std::endl;
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

void Efield::compute_force(ModuleBase::matrix &fdip)
{
	assert(bmod>0);
	int iat = 0;
	for(int it=0; it<GlobalC::ucell.ntype; ++it)
	{
		for(int ia=0; ia<GlobalC::ucell.atoms[it].na; ++ia)
		{
			for(int jj=0; jj<3; ++jj)
			{
				fdip(iat, jj) = ModuleBase::e2 * Efield::eamp * GlobalC::ucell.atoms[it].zv * Efield::bvec[jj]/Efield::bmod;
			}
			++iat;
		}
	}
}
