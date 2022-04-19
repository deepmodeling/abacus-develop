#include "module_base/global_variable.h"
#include "module_orbital/ORB_read.h"
#include "src_pw/occupy.h"
#include "src_pw/klist.h"
#include "src_pw/magnetism.h"
#include "src_pw/wf_atomic.h"
#include "src_pw/wavefunc.h"
#include "src_pw/charge_broyden.h"
#include "src_pw/potential.h"
#include "module_cell/atom_pseudo.h"
#include "module_cell/atom_spec.h"
#include "module_cell/unitcell_pseudo.h"
#include "module_cell/pseudo_nc.h"
#include "module_symmetry/symmetry_basic.h"
#include "module_symmetry/symmetry.h"
#include "src_parallel/parallel_grid.h"
#include "src_parallel/parallel_kpoints.h"
#include "src_parallel/ft.h"
#include "src_pw/use_fft.h"
#include "src_pw/pw_complement.h"
#include "src_pw/pw_basis.h"
#include "src_pw/VNL_in_pw.h"
#include "src_pw/hamilt.h"
#include "src_parallel/ft.h"
#include "src_parallel/parallel_pw.h"
#include "input.h"
#include "src_pw/energy.h"
#include "module_xc/xc_functional.h"

bool ModuleSymmetry::Symmetry::symm_flag;

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
Charge_Mixing::Charge_Mixing(){}
Charge_Mixing::~Charge_Mixing(){}
Charge_Pulay::Charge_Pulay(){}
Charge_Pulay::~Charge_Pulay(){}
Charge_Broyden::Charge_Broyden(){}
Charge_Broyden::~Charge_Broyden(){}
Potential::Potential(){}
Potential::~Potential(){}
InfoNonlocal::InfoNonlocal(){}
InfoNonlocal::~InfoNonlocal(){}
UnitCell::UnitCell(){}
UnitCell::~UnitCell(){}
UnitCell_pseudo::UnitCell_pseudo(){}
UnitCell_pseudo::~UnitCell_pseudo(){}
Parallel_Grid::Parallel_Grid(){}
Parallel_Grid::~Parallel_Grid(){}
Use_FFT::Use_FFT(){}
Use_FFT::~Use_FFT(){}
WF_igk::WF_igk(){}
WF_igk::~WF_igk(){}
WF_atomic::WF_atomic(){}
WF_atomic::~WF_atomic(){}
wavefunc::wavefunc(){}
wavefunc::~wavefunc(){}
Parallel_PW::Parallel_PW(){}
Parallel_PW::~Parallel_PW(){}
Magnetism::Magnetism(){}
Magnetism::~Magnetism(){}
ORB_gaunt_table::ORB_gaunt_table(){}
ORB_gaunt_table::~ORB_gaunt_table(){}
pseudopot_cell_vl::pseudopot_cell_vl(){}
pseudopot_cell_vl::~pseudopot_cell_vl(){}
pseudopot_cell_vnl::pseudopot_cell_vnl(){}
pseudopot_cell_vnl::~pseudopot_cell_vnl(){}
Hamilt_PW::Hamilt_PW(){}
Hamilt_PW::~Hamilt_PW(){}
Hamilt::Hamilt(){}
Hamilt::~Hamilt(){}
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
Charge_Broyden CHR;
Potential pot;
UnitCell_pseudo ucell;
ModuleSymmetry::Symmetry symm;
Parallel_Grid Pgrid;
Use_FFT UFFT;
PW_Basis pw;
pseudopot_cell_vnl ppcell;
Hamilt hm;
energy en;
Parallel_Kpoints Pkpoints;
} // namespace GlobalC
Input INPUT;

void Occupy::calculate_weights()
{
	GlobalC::wf.wg(0,0)=2.0;
	GlobalC::wf.wg(0,1)=0.0;
	GlobalC::wf.wg(0,2)=0.0;
	GlobalC::wf.wg(0,3)=0.0;
}

void Use_FFT::allocate()
{
    delete[] porter;
    porter = new std::complex<double>[GlobalC::pw.nrxx];
    return;
}

void wavefunc::allocate(const int nks)
{
	this->npwx = GlobalC::pw.setupIndGk(this->igk, GlobalC::kv.ngk);
	this->wg.create(nks,GlobalV::NBANDS);
	this->evc=new ModuleBase::ComplexMatrix[nks];
	this->ekb = new double*[nks];
	for (int ik=0;ik<nks;ik++)
	{
		this->ekb[ik] = new double[GlobalV::NBANDS];
		this->evc[ik].create(GlobalV::NBANDS,npwx);
	}
	return;
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

	ModuleBase::CHECK_INT(ifs, GlobalC::pw.ncx);	
	ModuleBase::CHECK_INT(ifs, GlobalC::pw.ncy);	
	ModuleBase::CHECK_INT(ifs, GlobalC::pw.ncz);	

	const int nxy = GlobalC::pw.ncx * GlobalC::pw.ncy;
	double *zpiece = new double[nxy];
	for(int iz=0; iz<GlobalC::pw.ncz; iz++)
	{
		ModuleBase::GlobalFunc::ZEROS(zpiece, nxy);
		for(int j=0; j<GlobalC::pw.ncy; j++)
		{
			for(int i=0; i<GlobalC::pw.ncx; i++)
			{
				ifs >> zpiece[ i*GlobalC::pw.ncy + j ];
			}
		}

		for(int ir=0; ir<nxy; ir++)
		{
			rho[ir*GlobalC::pw.nczp+iz] = zpiece[ir];
		}
	}// iz
	delete[] zpiece;

	ifs.close();
	return true;
}

void Use_FFT::ToRealSpace(int const&is, ModuleBase::ComplexMatrix const&vg, double*vr){}
void Use_FFT::ToRealSpace(std::complex<double> const*vg, double*vr){}
bool Occupy::use_gaussian_broadening=false;
bool Occupy::use_tetrahedron_method = false;
double Magnetism::get_nelup(){return 0;}
double Magnetism::get_neldw(){return 0;}
void PW_Basis::bspline_sf(const int norder){}

bool ModuleSymmetry::Symmetry_Basic::equal(double const&m, double const&n) const{return false;}

void UnitCell_pseudo::setup_cell(
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

int UnitCell_pseudo::read_atom_species(LCAO_Orbitals &orb, std::ifstream &ifa, std::ofstream &ofs_running)
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

bool UnitCell_pseudo::read_atom_positions(LCAO_Orbitals &orb, std::ifstream &ifpos, std::ofstream &ofs_running, std::ofstream &ofs_warning)
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
