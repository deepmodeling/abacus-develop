#include "module_base/global_variable.h"
#include "module_orbital/ORB_read.h"
#include "src_pw/occupy.h"
#include "src_pw/klist.h"
#include "src_pw/magnetism.h"
#include "src_pw/wf_atomic.h"
#include "src_pw/wavefunc.h"
#include "src_pw/charge_mixing.h"
#include "module_elecstate/potentials/potential_new.h"
#include "module_cell/atom_pseudo.h"
#include "module_cell/atom_spec.h"
#include "module_cell/unitcell.h"
#include "module_cell/pseudo_nc.h"
#include "module_symmetry/symmetry_basic.h"
#include "module_symmetry/symmetry.h"
#include "src_parallel/parallel_grid.h"
#include "src_parallel/parallel_kpoints.h"
#include "src_pw/pw_complement.h"
#include "src_pw/structure_factor.h"
#include "src_pw/VNL_in_pw.h"
#include "input.h"
#include "src_pw/energy.h"
#include "module_xc/xc_functional.h"
#include "module_pw/pw_basis_k.h"
#include "src_io/restart.h"

int ModuleSymmetry::Symmetry::symm_flag;

LCAO_Orbitals::LCAO_Orbitals(){}
LCAO_Orbitals::~LCAO_Orbitals(){}

ModuleSymmetry::Symmetry::Symmetry(){}
ModuleSymmetry::Symmetry::~Symmetry(){}
ModuleSymmetry::Symmetry_Basic::Symmetry_Basic(){}
ModuleSymmetry::Symmetry_Basic::~Symmetry_Basic(){}

pseudo_nc::pseudo_nc(){}
pseudo_nc::~pseudo_nc(){}
Atom::Atom(){}
Atom::~Atom(){}
Atom_pseudo::Atom_pseudo(){}
Atom_pseudo::~Atom_pseudo(){}
InfoNonlocal::InfoNonlocal(){}
InfoNonlocal::~InfoNonlocal(){}
UnitCell::UnitCell(){}
UnitCell::~UnitCell(){}
Parallel_Grid::Parallel_Grid(){}
Parallel_Grid::~Parallel_Grid(){}
WF_igk::WF_igk(){}
WF_igk::~WF_igk(){}
WF_atomic::WF_atomic(){}
WF_atomic::~WF_atomic(){}
wavefunc::wavefunc(){}
wavefunc::~wavefunc(){}
Magnetism::Magnetism(){}
Magnetism::~Magnetism(){}
ORB_gaunt_table::ORB_gaunt_table(){}
ORB_gaunt_table::~ORB_gaunt_table(){}
pseudopot_cell_vl::pseudopot_cell_vl(){}
pseudopot_cell_vl::~pseudopot_cell_vl(){}
pseudopot_cell_vnl::pseudopot_cell_vnl(){}
pseudopot_cell_vnl::~pseudopot_cell_vnl(){}
energy::energy(){}
energy::~energy(){}


XC_Functional::XC_Functional(){}

XC_Functional::~XC_Functional(){}

int XC_Functional::func_type = 0;

int XC_Functional::get_func_type()
{
    return func_type;
}


namespace GlobalC
{
K_Vectors kv;
wavefunc wf;
Charge CHR;
UnitCell ucell;
ModuleSymmetry::Symmetry symm;
Parallel_Grid Pgrid;
Structure_Factor sf;
ModulePW::PW_Basis* rhopw;
ModulePW::PW_Basis_K* wfcpw;
pseudopot_cell_vnl ppcell;
energy en;
Parallel_Kpoints Pkpoints;
Restart restart;
} // namespace GlobalC
Input INPUT;


void Restart::load_disk(const std::string mode, const int i, double** rho) const {}


psi::Psi<complex<double>>* wavefunc::allocate(const int nks)
{
	this->npwx = GlobalC::wfcpw->npwk_max;
	psi::Psi<std::complex<double>>* psi = new psi::Psi<std::complex<double>>(nks, GlobalV::NBANDS,npwx, nullptr);
	return psi;
}

bool Charge::read_rho(const int &is, const std::string &fn, double* rho) //add by dwan
{
	std::ifstream ifs(fn.c_str());

	std::string name;
	ifs >> name;
    
	bool quit=false;
	ModuleBase::CHECK_DOUBLE(ifs,GlobalC::ucell.lat0 * 0.529177,quit);
	ModuleBase::CHECK_DOUBLE(ifs,GlobalC::ucell.latvec.e11,quit);
	ModuleBase::CHECK_DOUBLE(ifs,GlobalC::ucell.latvec.e12,quit);
	ModuleBase::CHECK_DOUBLE(ifs,GlobalC::ucell.latvec.e13,quit);
	ModuleBase::CHECK_DOUBLE(ifs,GlobalC::ucell.latvec.e21,quit);
	ModuleBase::CHECK_DOUBLE(ifs,GlobalC::ucell.latvec.e22,quit);
	ModuleBase::CHECK_DOUBLE(ifs,GlobalC::ucell.latvec.e23,quit);
	ModuleBase::CHECK_DOUBLE(ifs,GlobalC::ucell.latvec.e31,quit);
	ModuleBase::CHECK_DOUBLE(ifs,GlobalC::ucell.latvec.e32,quit);
	ModuleBase::CHECK_DOUBLE(ifs,GlobalC::ucell.latvec.e33,quit);
	for(int it=0; it<GlobalC::ucell.ntype; it++)
	{
		ModuleBase::CHECK_STRING(ifs,GlobalC::ucell.atoms[it].label,quit);
	}

	for(int it=0; it<GlobalC::ucell.ntype; it++)
	{
		ModuleBase::CHECK_DOUBLE(ifs,GlobalC::ucell.atoms[it].na,quit);
	}

	std::string coordinate;
	ifs >> coordinate;
	double tau;
	for(int it=0; it<GlobalC::ucell.ntype; it++)
	{
		for(int ia=0; ia<GlobalC::ucell.atoms[it].na; ia++)
		{
			ifs>>tau;
			ifs>>tau;
			ifs>>tau;
			/*
			ModuleBase::CHECK_DOUBLE(ifs,GlobalC::ucell.atoms[it].taud[ia].x,quit);
			ModuleBase::CHECK_DOUBLE(ifs,GlobalC::ucell.atoms[it].taud[ia].y,quit);
			ModuleBase::CHECK_DOUBLE(ifs,GlobalC::ucell.atoms[it].taud[ia].z,quit);
			*/
		}
	}

	ModuleBase::CHECK_INT(ifs, GlobalV::NSPIN);
	ModuleBase::GlobalFunc::READ_VALUE(ifs, GlobalC::en.ef);

	ModuleBase::CHECK_INT(ifs, GlobalC::rhopw->nx);	
	ModuleBase::CHECK_INT(ifs, GlobalC::rhopw->ny);	
	ModuleBase::CHECK_INT(ifs, GlobalC::rhopw->nz);	

	const int nxy = GlobalC::rhopw->nx * GlobalC::rhopw->ny;
	double *zpiece = new double[nxy];
	for(int iz=0; iz<GlobalC::rhopw->nz; iz++)
	{
		ModuleBase::GlobalFunc::ZEROS(zpiece, nxy);
		for(int j=0; j<GlobalC::rhopw->ny; j++)
		{
			for(int i=0; i<GlobalC::rhopw->nx; i++)
			{
				ifs >> zpiece[ i*GlobalC::rhopw->ny + j ];
			}
		}

		for(int ir=0; ir<nxy; ir++)
		{
			rho[ir*GlobalC::rhopw->nplane+iz] = zpiece[ir];
		}
	}// iz
	delete[] zpiece;

	ifs.close();
	return true;
}

//bool Occupy::use_gaussian_broadening=false;

bool ModuleSymmetry::Symmetry_Basic::equal(double const&m, double const&n) const{return false;}

void UnitCell::setup_cell(
#ifdef __LCAO
		LCAO_Orbitals &orb,
#endif
		const std::string &s_pseudopot_dir,
		const std::string &fn,
		std::ofstream &log)
{
	// (1) init mag
	assert(ntype>0);
	// (2) init *Atom class array.
	this->atoms = new Atom[this->ntype]; // atom species.
	this->set_atom_flag = true;


	bool ok = true;
	bool ok2 = true;

	// (3) read in atom information
		// open "atom_unitcell" file.
	std::ifstream ifa(fn.c_str(), ios::in);

	//========================
	// call read_atom_species
	//========================
	const int error = this->read_atom_species(orb, ifa, log);


	//==========================
	// call read_atom_positions
	//==========================
	ok2 = this->read_atom_positions(orb, ifa, log, GlobalV::ofs_warning);
	//========================================================
	// Calculate unit cell volume
	// the reason to calculate volume here is 
	// Firstly, latvec must be read in.
	//========================================================
	assert(lat0 > 0.0);
	this->omega = abs( latvec.Det() ) * this->lat0 * lat0 * lat0 ;
		
	//==========================================================
	// Calculate recip. lattice vectors and dot products
	// latvec have the unit of lat0, but G has the unit 2Pi/lat0
	//==========================================================
	this->GT = latvec.Inverse();
	this->G  = GT.Transpose();
	this->GGT = G * GT;
	this->invGGT = GGT.Inverse();

	//this->cal_meshx();

	return;
}

int UnitCell::read_atom_species(LCAO_Orbitals &orb, std::ifstream &ifa, std::ofstream &ofs_running)
{
	delete[] atom_label;
	delete[] atom_mass;
	delete[] pseudo_fn;
	this->atom_mass  = new double[ntype]; //atom masses
	this->atom_label = new std::string[ntype]; //atom labels
	this->pseudo_fn  = new std::string[ntype]; //file name of pseudopotential
	if( ModuleBase::GlobalFunc::SCAN_BEGIN(ifa, "ATOMIC_SPECIES") )
	{	
		ModuleBase::GlobalFunc::OUT(ofs_running,"ntype",ntype);
		for (int i = 0;i < ntype;i++)
		{
			ifa >> this->atom_label[i] >> this->atom_mass[i];
		}
	}
	if( ModuleBase::GlobalFunc::SCAN_BEGIN(ifa, "LATTICE_CONSTANT") )
	{
		ModuleBase::GlobalFunc::READ_VALUE(ifa, lat0);
	}
	lat0_angstrom = lat0 * 0.529177 ;
	this->tpiba  = ModuleBase::TWO_PI / lat0;
	this->tpiba2 = tpiba * tpiba;
	if(latName=="sc"){//simple-cubic
		latvec.e11 = 1.0; latvec.e12 = 0.0; latvec.e13 = 0.0;
		latvec.e21 = 0.0; latvec.e22 = 1.0;	latvec.e23 = 0.0;
		latvec.e31 = 0.0; latvec.e32 = 0.0;	latvec.e33 = 1.0;
	}
	a1.x = latvec.e11;
	a1.y = latvec.e12;
	a1.z = latvec.e13;

	a2.x = latvec.e21;
	a2.y = latvec.e22;
	a2.z = latvec.e23;

	a3.x = latvec.e31;
	a3.y = latvec.e32;
	a3.z = latvec.e33;
	return 0;
}

bool UnitCell::read_atom_positions(LCAO_Orbitals &orb, std::ifstream &ifpos, std::ofstream &ofs_running, std::ofstream &ofs_warning)
{
	if( ModuleBase::GlobalFunc::SCAN_BEGIN(ifpos, "ATOMIC_POSITIONS"))
	{
		ModuleBase::GlobalFunc::READ_VALUE( ifpos, Coordinate);

		ModuleBase::Vector3<double> v;
		ModuleBase::Vector3<int> mv;
		int na = 0;
		this->nat = 0;
		assert(ntype>0);
		for (int it = 0;it < ntype; it++)
		{
			//=======================================
			// (1) read in atom label
			// start magnetization
			//=======================================
			ModuleBase::GlobalFunc::READ_VALUE(ifpos, atoms[it].label);
			bool found = false;
			for(int it2=0; it2<ntype; it2++)
			{
				if( this->atoms[it].label == this->atom_label[it] )
				{
					found = true;
				}
			}
			double magnet;
			ModuleBase::GlobalFunc::READ_VALUE(ifpos, magnet );
			//=========================
			// (3) read in atom number
			//=========================
			ModuleBase::GlobalFunc::READ_VALUE(ifpos, na);
			this->atoms[it].na = na;
			this->nat += na;
		}
	}
        return true;
}

void elecstate::Potential::init_pot(int istep, const Charge* chg)
{
	return;
}