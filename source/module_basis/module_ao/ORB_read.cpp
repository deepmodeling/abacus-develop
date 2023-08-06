#include "ORB_read.h"
#include <cstring>		// Peize Lin fix bug about strcmp 2016-08-02
#include <cassert>
#include "module_base/math_integral.h"
#include "module_base/tool_check.h"
#include "module_base/parallel_common.h"
#include <algorithm>
#include "module_base/timer.h"

//==============================
// Define an object here! 
//==============================
/// PLEASE avoid using 'ORB' as global variable 
// mohan note 2021-03-23
namespace GlobalC
{
LCAO_Orbitals ORB;
}

LCAO_Orbitals::LCAO_Orbitals()
{
	this->nchimax = 0;// this initialzied must specified
	this->Phi = new Numerical_Orbital[1];	
	this->Alpha = new Numerical_Orbital[1];

	this->read_in_flag = false;	

	this->dr_uniform = 0.001;

	this->lmax_d = 0;
    this->nchimax_d = 0;
    this->rcutmax_Phi = 0.0;
}

LCAO_Orbitals::~LCAO_Orbitals()
{
	delete[] Phi;
	delete[] Alpha;
}

const LCAO_Orbitals& LCAO_Orbitals::get_const_instance()
{
	return GlobalC::ORB;
}

#ifdef __MPI
// be called in UnitCell.
void LCAO_Orbitals::bcast_files(
	const int &ntype_in, 
	const int &my_rank)
{
	ModuleBase::TITLE("LCAO_Orbitals","bcast_files");

	// 'read_in_flag' is true when there is a
	// block "NUMERICAL_ORBITAL" in structure
	// file.
	Parallel_Common::bcast_bool(read_in_flag);
    Parallel_Common::bcast_string(descriptor_file);
	if(!read_in_flag)
	{
		return;
	}

	assert(ntype_in > 0 );

	GlobalV::ofs_running << "\n READING ORBITAL FILE NAMES FOR LCAO" << std::endl;
	for(int it=0; it<ntype_in; it++)
	{
		std::string ofile;
		std::string nfile;

		if(my_rank==0)
		{
			ofile = orbital_file[it];
			//-----------------------------------
			// Turn off the read in NONLOCAL file
			// function since 2013-08-02 by mohan
			//-----------------------------------
//			nfile = nonlocal_file[it];
		}

// PLEASE avoid using 'bcast_string' as global variable 
// mohan note 2021-03-23
		Parallel_Common::bcast_string(ofile);
		//-----------------------------------
		// Turn off the read in NONLOCAL file
		// function since 2013-08-02 by mohan
		//-----------------------------------
//		Parallel_Common::bcast_string(nfile);

		if(my_rank!=0)
		{
			orbital_file.push_back( ofile );
			//-----------------------------------
			// Turn off the read in NONLOCAL file
			// function since 2013-08-02 by mohan
			//-----------------------------------
//			nonlocal_file.push_back ( nfile );
		}

		GlobalV::ofs_running << " orbital file: " << orbital_file[it] << std::endl;
//		GlobalV::ofs_running << " nonlocal file: " << nonlocal_file[it] << std::endl;
	}
	return;
}
#endif


void LCAO_Orbitals::Read_Orbitals(
	std::ofstream &ofs_in,
	const int &ntype_in, 
	const int &lmax_in,
	const bool &deepks_setorb,
	const int &out_mat_r,
	const bool &force_flag, // mohan add 2021-05-07
	const int &my_rank) // mohan add 2021-04-26
{
	ModuleBase::TITLE("LCAO_Orbitals", "Read_Orbitals");
	ModuleBase::timer::tick("LCAO_Orbitals","Read_Orbitals");

	ofs_in << "\n\n\n\n";
	ofs_in << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
	ofs_in << " |                                                                    |" << std::endl;
	ofs_in << " | Setup numerical orbitals:                                          |" << std::endl;
	ofs_in << " | This part setup: numerical atomic orbitals, non-local projectors   |" << std::endl;
	ofs_in << " | and neutral potential (1D). The atomic orbitals information        |" << std::endl;
	ofs_in << " | including the radius, angular momentum and zeta number.            |" << std::endl;
	ofs_in << " | The neutral potential is the sum of local part of pseudopotential  |" << std::endl;
	ofs_in << " | and potential given by atomic charge, they will cancel out beyond  |" << std::endl;
	ofs_in << " | a certain radius cutoff, because the Z/r character.                |" << std::endl;
	ofs_in << " |                                                                    |" << std::endl;
	ofs_in << " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
	ofs_in << "\n\n\n\n";	

	//////////////////////
	/// (1) check dk, dR, Rmax.
	//////////////////////

	ofs_in << "\n SETUP ONE DIMENSIONAL ORBITALS/POTENTIAL" << std::endl;

	if(!read_in_flag)
	{
		ModuleBase::WARNING_QUIT("LCAO_Orbitals::Read_Orbitals","Set the NUMERICAL_ORBITAL block in structure file.");
	}


	//OUT(ofs_in,"ecutwfc for kmesh",ecutwfc);
	ModuleBase::GlobalFunc::OUT(ofs_in,"delta k  (1/Bohr)",dk);
	ModuleBase::GlobalFunc::OUT(ofs_in,"delta r    (Bohr)",dR);
	ModuleBase::GlobalFunc::OUT(ofs_in,"dr_uniform (Bohr)",dr_uniform);
	ModuleBase::GlobalFunc::OUT(ofs_in,"rmax       (Bohr)",Rmax);

	// check the read in data.
    assert(dk > 0.0);
    assert(ecutwfc > 0.0);
    assert(dR > 0.0);
    assert(Rmax > 0.0);

	/// ntype: number of atom species
	this->ntype = ntype_in; 
	assert(ntype>0);

	/// lmax: lmax used in local orbitals as basis sets
	assert(lmax_in>=0); // mohan add 2021-04-16
	this->lmax = lmax_in;

	//////////////////////////////////////////////////////////
	/// (2) set the kmesh according to ecutwfc and dk. 
	//////////////////////////////////////////////////////////

	//-----------------------------------------------------------------
	/// calculate number of k mesh according to energy cutoff.
	/// Mohan choose ecutwfc according to interpolation requirement.
	//	std::cout << " ecutwfc=" << ecutwfc << std::endl;
	//LiuXh modified 2016-01-25, 2016-07-20
	if(ecutwfc< 20)
	{
		this->kmesh = static_cast<int>( 2 * sqrt(ecutwfc) / dk )  + 4;
	}
	else
	{
		this->kmesh = static_cast<int>( sqrt(ecutwfc) / dk )  + 4;
	}

	// jingan add for calculate r(R) matrix
	// if(out_mat_r) 
	// {
	// 	kmesh = kmesh * 4;
	// }

	//	this->kmesh = static_cast<int> (PI / 0.01 / 4 / this->dk);
	if(kmesh%2==0) kmesh++;
	ModuleBase::GlobalFunc::OUT(ofs_in,"kmesh",kmesh);
	//-----------------------------------------------------------------


	//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	//~~~~~~~~~~~~~~~~~~~~~~   1    ~~~~~~~~~~~~~~~~~~~~~~~~~
	// Read in numerical atomic orbitals for each atom type.
	//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	delete[] this->Phi;

	this->Phi = new Numerical_Orbital[ntype];
	for(int it=0; it<ntype; it++)
	{
		this->Read_PAO(ofs_in, it, force_flag, my_rank);
		//caoyu add 2021-05-24	to reconstruct atom_arrange::set_sr_NL
		this->rcutmax_Phi = std::max(this->rcutmax_Phi, this->Phi[it].getRcut());
	}

	

	//caoyu add 2021-3-16
	///>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	///~~~~~~~~~~~~~~~~~~~~~~   3    ~~~~~~~~~~~~~~~~~~~~~~~~~
	/// Read in numerical basis for descriptor.
	///>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


	if (deepks_setorb)	//condition: descriptor in lcao line
	{
		
		delete[] this->Alpha;
		this->Alpha = new Numerical_Orbital[1];	//not related to atom type -- remain to be discussed

		this->Read_Descriptor(ofs_in, force_flag, my_rank);

	}

	ModuleBase::timer::tick("LCAO_Orbitals","Read_Orbitals");
	return;
}





//-------------------------------------------------------
// mohan note 2021-04-26
// to_caoyu: 
// 1. read in lmaxt and nchi directly from orbital files
// 2. pass nchi to phi via this->Phi[it].set_orbital_info 
// be careful! nchi[l] may be different for differnt phi
//-------------------------------------------------------
void LCAO_Orbitals::Read_PAO(
	std::ofstream &ofs_in,
	const int& it, 
	const bool &force_flag, // mohan add 2021-05-07
	const int &my_rank) // mohan add 2021-04-26
{
	ModuleBase::TITLE("LCAO_Orbitals","Read_PAO");

	std::ifstream in_ao;
	bool open=false;
	if(my_rank==0)
	{
		in_ao.open(this->orbital_file[it].c_str());
		if(in_ao)
		{
			open=true;
		}
	}
#ifdef __MPI
	Parallel_Common::bcast_bool( open );
#endif
	if(!open)
	{
		std::cout << " Orbital file : " << this->orbital_file[it] << std::endl;
		ModuleBase::WARNING_QUIT("LCAO_Orbitals::Read_PAO","Couldn't find orbital files");
	}

	ofs_in << " " << std::setw(12) << "ORBITAL" << std::setw(3) << "L" 
	<< std::setw(3) << "N" << std::setw(8) << "nr" << std::setw(8) << "dr"
	<< std::setw(8) << "RCUT" << std::setw(12) << "CHECK_UNIT"
		<< std::setw(12) << "NEW_UNIT" << std::endl;
	
	//lmax and nchimax for type it
	int lmaxt=0;
	int nchimaxt=0;

	this->read_orb_file(ofs_in, in_ao, it, lmaxt, nchimaxt, this->Phi, force_flag, my_rank);

	//lmax and nchimax for all types
	this->lmax = std::max(this->lmax, lmaxt);
	this->nchimax = std::max(this->nchimax, nchimaxt);
	
	in_ao.close();
	return;
}


//caoyu add 2021-3-16
void LCAO_Orbitals::Read_Descriptor(
	std::ofstream &ofs_in,
	const bool &force_flag, // mohan add 2021-05-07
	const int &my_rank)	//read descriptor basis
{
	ModuleBase::TITLE("LCAO_Orbitals", "Read_Descriptor");

	std::ifstream in_de;
	ofs_in << " " << std::setw(12) << "DESCRIPTOR" << std::setw(3) << "L"
		<< std::setw(3) << "N" << std::setw(8) << "nr" << std::setw(8) << "dr"
		<< std::setw(8) << "RCUT" << std::setw(12) << "CHECK_UNIT"
		<< std::setw(12) << "NEW_UNIT" << std::endl;

	// check if the descriptor file exists.
	bool open = false;
	if (my_rank == 0)
	{
		in_de.open(this->descriptor_file.c_str());
		if (in_de)
		{
			open = true;
		}
	}
#ifdef __MPI
	Parallel_Common::bcast_bool(open);
#endif
	if (!open)
	{
		std::cout << " Orbital file : " << this->descriptor_file << std::endl;
		ModuleBase::WARNING_QUIT("LCAO_Orbitals::Read_Descriptor", "Couldn't find orbital files for descriptor");
	}

	this->lmax_d = 0;
	this->nchimax_d = 0;

	this->read_orb_file(ofs_in, in_de, 0, this->lmax_d, this->nchimax_d, this->Alpha, force_flag, my_rank);

	in_de.close();

	return;
}


void LCAO_Orbitals::read_orb_file(
	std::ofstream &ofs_in, // GlobalV::ofs_running
	std::ifstream &ifs,
	const int &it, 
	int &lmax, 
	int &nchimax, 
	Numerical_Orbital* ao,
	const bool &force_flag,
	const int &my_rank)
{
	ModuleBase::TITLE("LCAO_Orbitals","read_orb_file");
	char word[80];
	std::string orb_label;
	if (my_rank == 0)
	{
		while (ifs.good())
		{
			ifs >> word;
			if (std::strcmp(word, "Element") == 0)
			{
				ifs >> orb_label;
				continue;
			}
			if (std::strcmp(word, "Lmax") == 0)
			{
				ifs >> lmax;
				break;
			}
		}
	}
#ifdef __MPI
	Parallel_Common::bcast_int(lmax);
#endif	
	
	int* nchi = new int[lmax+1];		// allocate space: number of chi for each L.
	
	if (my_rank == 0)
	{	
		for (int l = 0; l <= lmax; l++)
		{
			ifs >> word >> word >> word >> nchi[l];
			nchimax = std::max(nchimax, nchi[l]);
		}
	}

#ifdef __MPI
	Parallel_Common::bcast_int(nchimax);
	Parallel_Common::bcast_int(nchi, lmax + 1);
#endif		

	// calculate total number of chi
	int total_nchi = 0;
	for (int l = 0; l <= lmax; l++)
	{
		total_nchi += nchi[l];
	}

	//OUT(GlobalV::ofs_running,"Total number of chi(l,n)",total_nchi);
	delete[] ao[it].phiLN;
	ao[it].phiLN = new Numerical_Orbital_Lm[total_nchi];

	int meshr=0; // number of mesh points
	int meshr_read=0;
	double dr=0.0; 

	if (my_rank == 0)
	{
		while (ifs.good())
		{
			ifs >> word;
			if (std::strcmp(word, "END") == 0)		// Peize Lin fix bug about strcmp 2016-08-02
			{
				break;
			}
		}
		ModuleBase::CHECK_NAME(ifs, "Mesh");
		ifs >> meshr;
		meshr_read = meshr;
		if (meshr % 2 == 0)
		{
			++meshr;
		}
		ModuleBase::CHECK_NAME(ifs, "dr");
		ifs >> dr;
	}

#ifdef __MPI
	Parallel_Common::bcast_int(meshr);
	Parallel_Common::bcast_int(meshr_read);
	Parallel_Common::bcast_double(dr);
#endif		

	int count = 0;
	std::string name1;
	std::string name2;
	std::string name3;
	int tmp_it=0;
	int tmp_l=0;
	int tmp_n=0;

	for (int L = 0; L <= lmax; L++)
	{
		for (int N = 0; N < nchi[L]; N++)
		{
			ofs_in << " " << std::setw(12) << count + 1 << std::setw(3) << L << std::setw(3) << N;

			double* radial; // radial mesh
			double* psi; // radial local orbital
			double* psir;// psi * r
			double* rab;// dr

			// set the number of mesh and the interval distance.
			ofs_in << std::setw(8) << meshr << std::setw(8) << dr;

			radial = new double[meshr];
			psi = new double[meshr];
			psir = new double[meshr];
			rab = new double[meshr];

			for(int im=0; im<meshr; ++im)
			{
				radial[im]=0.0;
				psi[im]=0.0;
				psir[im]=0.0;
				rab[im]=0.0;
			}

			for (int ir = 0; ir < meshr; ir++)
			{
				rab[ir] = dr;
				// plus one because we can't read in r = 0 term now.
				// change ir+1 to ir, because we need psi(r==0) information.
				radial[ir] = ir * dr; //mohan 2010-04-19
			}

			// set the length of orbital
			ofs_in << std::setw(8) << radial[meshr - 1];

			// mohan update 2010-09-07
			bool find = false;
			if (my_rank == 0)
			{
				while (!find)
				{
					if (ifs.eof())
					{
						std::cout << " Can't find l="
							<< L << " n=" << N << " orbital." << std::endl;
						break;
					}

					ifs >> name1 >> name2>> name3;
					ifs >> tmp_it >> tmp_l >> tmp_n;
					assert( name1 == "Type" );
					if (L == tmp_l && N == tmp_n)
					{
						// meshr_read is different from meshr if meshr is even number.
						for (int ir = 0; ir < meshr_read; ir++)
						{
							ifs >> psi[ir];
							psir[ir] = psi[ir] * radial[ir];
						}
						find = true;
					}
					else
					{
						double no_use;
						for (int ir = 0; ir < meshr_read; ir++)
						{
							ifs >> no_use;
						}
					}
				}//end find
			}

#ifdef __MPI
			Parallel_Common::bcast_bool(find);
#endif
			if (!find)
			{
				ModuleBase::WARNING_QUIT("LCAO_Orbitals::read_orb_file", "Can't find orbitals.");
			}

#ifdef __MPI
			Parallel_Common::bcast_double(psi, meshr_read);
			Parallel_Common::bcast_double(psir, meshr_read);
#endif

			// renormalize radial wave functions
			double* inner = new double[meshr]();
			for (int ir = 0; ir < meshr; ir++)
			{
				inner[ir] = psir[ir] * psir[ir];
			}
			double unit = 0.0;

			ModuleBase::Integral::Simpson_Integral(meshr, inner, rab, unit);

			assert(unit>0.0);

			// check unit: \sum ( psi[r] * r )^2 = 1
			ofs_in << std::setprecision(3) << std::setw(12) << unit;

			for (int ir = 0; ir < meshr; ir++)
			{
				psi[ir] /= sqrt(unit);
				psir[ir] /= sqrt(unit);
			}

			for (int ir = 0; ir < meshr; ir++)
			{
				inner[ir] = psir[ir] * psir[ir];
			}
			ModuleBase::Integral::Simpson_Integral(meshr, inner, rab, unit);
			delete[] inner;
			ofs_in << std::setw(12) << unit << std::endl;

			ao[it].phiLN[count].set_orbital_info(
                orb_label,
                it, //type
				L, //angular momentum L
				N, // number of orbitals of this L
				meshr, // number of radial mesh
				rab,
				radial,// radial mesh value(a.u.)
				Numerical_Orbital_Lm::Psi_Type::Psi,// psi type next
				psi, // radial wave function
				this->kmesh,
				this->dk,
				this->dr_uniform,
				GlobalV::out_element_info,
				true,
				force_flag); // delta k mesh in reciprocal space

			delete[] radial;
			delete[] rab;
			delete[] psir;
			delete[] psi;

			++count;
		}
	}
	ao[it].set_orbital_info(
        it, // type
        orb_label, // label	
		lmax,
		nchi,
		total_nchi); //copy twice !
	
	delete[] nchi;
	return;
}
