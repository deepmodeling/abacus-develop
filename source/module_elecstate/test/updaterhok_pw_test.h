#include "module_base/global_variable.h"
#include "module_basis/module_ao/ORB_read.h"
#include "module_elecstate/occupy.h"
#include "module_cell/klist.h"
#include "module_elecstate/magnetism.h"
#include "module_hamilt_pw/hamilt_pwdft/wavefunc.h"
#include "module_elecstate/module_charge/charge_mixing.h"
#include "module_elecstate/potentials/potential_new.h"
#include "module_cell/atom_pseudo.h"
#include "module_cell/atom_spec.h"
#include "module_cell/unitcell.h"
#include "module_cell/pseudo_nc.h"
#include "module_cell/module_symmetry/symmetry_basic.h"
#include "module_cell/module_symmetry/symmetry.h"
#include "module_hamilt_pw/hamilt_pwdft/parallel_grid.h"
#include "module_cell/parallel_kpoints.h"
#include "module_hamilt_pw/hamilt_pwdft/structure_factor.h"
#include "module_hamilt_pw/hamilt_pwdft/VNL_in_pw.h"
#include "module_io/input.h"
#include "module_elecstate/energy.h"
#include "module_hamilt_general/module_xc/xc_functional.h"
#include "module_basis/module_pw/pw_basis_k.h"
#include "module_io/restart.h"
#include "module_io/rho_io.h"

LCAO_Orbitals::LCAO_Orbitals(){}
LCAO_Orbitals::~LCAO_Orbitals(){}

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
double& energy::get_ef(const int&is, const bool& two_efermi){return this->ef;} //just mock


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

bool ModuleIO::read_rho(const int &is,
	const int &nspin,
	const std::string &fn,
	double* rho,
	int& nx,
	int& ny,
	int& nz,
	double& ef,
	const UnitCell* ucell,
	int &prenspin)
{
	std::ifstream ifs(fn.c_str());
	bool quit=false;

	ifs.ignore(300, '\n'); // skip the header

	ModuleBase::CHECK_INT(ifs, nspin);
	ifs.ignore(150, ')');
	ifs >> ef;
	ifs.ignore(150, '\n');

	ModuleBase::CHECK_INT(ifs,ucell->nat,quit);
	ifs.ignore(150, '\n');

	double fac=ucell->lat0;
	ModuleBase::CHECK_INT(ifs,nx);	
	ModuleBase::CHECK_DOUBLE(ifs, fac*ucell->latvec.e11/double(nx), quit);
	ModuleBase::CHECK_DOUBLE(ifs, fac*ucell->latvec.e12/double(nx), quit);
	ModuleBase::CHECK_DOUBLE(ifs, fac*ucell->latvec.e13/double(nx), quit);
	ModuleBase::CHECK_INT(ifs, ny);	
	ModuleBase::CHECK_DOUBLE(ifs, fac*ucell->latvec.e21/double(ny), quit);
	ModuleBase::CHECK_DOUBLE(ifs, fac*ucell->latvec.e22/double(ny), quit);
	ModuleBase::CHECK_DOUBLE(ifs, fac*ucell->latvec.e23/double(ny), quit);
	ModuleBase::CHECK_INT(ifs, nz);	
	ModuleBase::CHECK_DOUBLE(ifs, fac*ucell->latvec.e31/double(nz), quit);
	ModuleBase::CHECK_DOUBLE(ifs, fac*ucell->latvec.e32/double(nz), quit);
	ModuleBase::CHECK_DOUBLE(ifs, fac*ucell->latvec.e33/double(nz), quit);

	int temp = 0;
	for(int it=0; it<ucell->ntype; it++)
	{
		for(int ia=0; ia<ucell->atoms[it].na; ia++)
		{
			ifs >> temp;
			ifs >> temp; 
			ifs >> temp; 
			ifs >> temp; 
			ifs >> temp; 
		}
	}

	for(int i=0; i<nx; i++)
	{
		for(int j=0; j<ny; j++)
		{
			for(int k=0; k<nz; k++)
			{
				ifs >> rho[k*nx*ny+i*ny+j];
			}
		}
	}

	ifs.close();
    return true;

}

//bool Occupy::use_gaussian_broadening=false;

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
