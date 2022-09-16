#include "input.h"
#include "module_base/global_variable.h"
#include "module_base/mathzone.h"
#include "module_cell/atom_pseudo.h"
#include "module_cell/atom_spec.h"
#include "module_cell/pseudo_nc.h"
#include "module_cell/unitcell_pseudo.h"
#include "module_orbital/ORB_read.h"
#include "module_symmetry/symmetry.h"
#include "module_symmetry/symmetry_basic.h"
#include "module_xc/xc_functional.h"
#include "src_parallel/parallel_grid.h"
#include "src_parallel/parallel_kpoints.h"
#include "src_parallel/parallel_pw.h"
#include "src_pw/VNL_in_pw.h"
#include "src_pw/charge_broyden.h"
#include "src_pw/energy.h"
#include "src_pw/hamilt.h"
#include "src_pw/klist.h"
#include "src_pw/magnetism.h"
#include "src_pw/occupy.h"
#include "src_pw/potential.h"
#include "src_pw/structure_factor.h"
#include "src_pw/pw_complement.h"
#include "src_pw/use_fft.h"
#include "src_pw/wavefunc.h"
#include "src_pw/wf_atomic.h"
#include "module_pw/pw_basis_k.h"

#ifdef __LCAO
#include "module_neighbor/sltk_atom_arrange.h"
#include "module_orbital/ORB_atomic.h"
#include "module_orbital/ORB_atomic_lm.h"
#include "module_orbital/ORB_control.h"
#include "module_orbital/ORB_gaunt_table.h"
#include "module_orbital/ORB_gen_tables.h"
#include "module_orbital/ORB_nonlocal.h"
#include "module_orbital/ORB_nonlocal_lm.h"
#include "module_orbital/ORB_table_alpha.h"
#include "module_orbital/ORB_table_beta.h"
#include "module_orbital/ORB_table_phi.h"
#include "module_orbital/parallel_orbitals.h"
#include "driver.h"
#include "src_io/wf_local.h"
#include "module_gint/gint_gamma.h"
#include "module_gint/gint_tools.h"
#include "module_gint/grid_bigcell.h"
#include "module_gint/grid_meshball.h"
#include "module_gint/grid_meshcell.h"
#include "module_gint/grid_meshk.h"
#include "module_gint/grid_technique.h"
#include "src_lcao/local_orbital_charge.h"
#include "src_lcao/local_orbital_wfc.h"
#include "src_parallel/parallel_global.h"
#include "src_io/output.h"
#include "module_cell/unitcell_pseudo.h"
#include "module_cell/read_pp.h"

Parallel_Orbitals::Parallel_Orbitals()
{
}
Parallel_Orbitals::~Parallel_Orbitals()
{
}
Local_Orbital_Charge::Local_Orbital_Charge()
{
}
Local_Orbital_Charge::~Local_Orbital_Charge()
{
}
Local_Orbital_wfc::Local_Orbital_wfc()
{
}
Local_Orbital_wfc::~Local_Orbital_wfc()
{
}

#endif

int ModuleSymmetry::Symmetry::symm_flag;

ModuleSymmetry::Symmetry::Symmetry()
{
}
ModuleSymmetry::Symmetry::~Symmetry()
{
}
ModuleSymmetry::Symmetry_Basic::Symmetry_Basic()
{
}
ModuleSymmetry::Symmetry_Basic::~Symmetry_Basic()
{
}

Charge::Charge()
{
}
Charge::~Charge()
{
}
Charge_Mixing::Charge_Mixing()
{
}
Charge_Mixing::~Charge_Mixing()
{
}
Charge_Pulay::Charge_Pulay()
{
}
Charge_Pulay::~Charge_Pulay()
{
}
Charge_Broyden::Charge_Broyden()
{
}
Charge_Broyden::~Charge_Broyden()
{
}
Potential::Potential()
{
}
Potential::~Potential()
{
}
UnitCell::UnitCell()
{
}
UnitCell::~UnitCell()
{
}
UnitCell_pseudo::UnitCell_pseudo()
{
}
UnitCell_pseudo::~UnitCell_pseudo()
{
}
Parallel_Grid::Parallel_Grid()
{
}
Parallel_Grid::~Parallel_Grid()
{
}
Use_FFT::Use_FFT()
{
}
Use_FFT::~Use_FFT()
{
}
WF_igk::WF_igk()
{
}
WF_igk::~WF_igk()
{
}
WF_atomic::WF_atomic()
{
}
WF_atomic::~WF_atomic()
{
}
wavefunc::wavefunc()
{
}
wavefunc::~wavefunc()
{
}
Parallel_PW::Parallel_PW()
{
}
Parallel_PW::~Parallel_PW()
{
}
Magnetism::Magnetism()
{
}
Magnetism::~Magnetism()
{
}
pseudopot_cell_vl::pseudopot_cell_vl()
{
}
pseudopot_cell_vl::~pseudopot_cell_vl()
{
}
pseudopot_cell_vnl::pseudopot_cell_vnl()
{
}
pseudopot_cell_vnl::~pseudopot_cell_vnl()
{
}
Hamilt_PW::Hamilt_PW()
{
}
Hamilt_PW::~Hamilt_PW()
{
}
Hamilt::Hamilt()
{
}
Hamilt::~Hamilt()
{
}
energy::energy()
{
}
energy::~energy()
{
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
Structure_Factor sf;
ModulePW::PW_Basis* rhopw;
ModulePW::PW_Basis_Big *bigpw = static_cast<ModulePW::PW_Basis_Big*>(rhopw);
ModulePW::PW_Basis_K* wfcpw;
pseudopot_cell_vnl ppcell;
Hamilt hm;
energy en;
Parallel_Kpoints Pkpoints;
// Grid_Technique GridT; already defined in ../grid_technique.cpp
#ifdef __LCAO
LCAO_Orbitals ORB;
Grid_Driver GridD(GlobalV::test_deconstructor, GlobalV::test_grid_driver, GlobalV::test_grid);
ORB_gen_tables UOT;
#endif
} // namespace GlobalC
Input INPUT;

void Occupy::calculate_weights()
{
    GlobalC::wf.wg(0, 0) = 2.0;
    GlobalC::wf.wg(0, 1) = 0.0;
    GlobalC::wf.wg(0, 2) = 0.0;
    GlobalC::wf.wg(0, 3) = 0.0;
}

void wavefunc::allocate_ekb_wg(const int nks)
{
    ModuleBase::TITLE("wavefunc", "init_local");
    this->npwx = GlobalC::wfcpw->npwk_max;
    this->ekb = new double *[nks];
    for (int ik = 0; ik < nks; ik++)
    {
        ekb[ik] = new double[GlobalV::NBANDS];
        ModuleBase::GlobalFunc::ZEROS(ekb[ik], GlobalV::NBANDS);
    }
    this->wg.create(nks, GlobalV::NBANDS);
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
			ModuleBase::CHECK_DOUBLE(ifs,GlobalC::ucell.atoms[it].taud[ia].x,quit);
			ModuleBase::CHECK_DOUBLE(ifs,GlobalC::ucell.atoms[it].taud[ia].y,quit);
			ModuleBase::CHECK_DOUBLE(ifs,GlobalC::ucell.atoms[it].taud[ia].z,quit);
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

void Charge::allocate(const int &nspin_in, const int &nrxx_in, const int &ngmc_in)
{

    this->nspin = nspin_in;
    this->nrxx = nrxx_in;
    this->ngmc = ngmc_in;

    rho = new double *[nspin];
    rhog = new std::complex<double> *[nspin];
    rho_save = new double *[nspin];
    rhog_save = new std::complex<double> *[nspin];
    for (int is = 0; is < nspin; is++)
    {
        rho[is] = new double[nrxx];
        rhog[is] = new std::complex<double>[ngmc];
        rho_save[is] = new double[nrxx];
        rhog_save[is] = new std::complex<double>[ngmc];
        ModuleBase::GlobalFunc::ZEROS(rho[is], nrxx);
        ModuleBase::GlobalFunc::ZEROS(rhog[is], ngmc);
        ModuleBase::GlobalFunc::ZEROS(rho_save[is], nrxx);
        ModuleBase::GlobalFunc::ZEROS(rhog_save[is], ngmc);
    }

    this->rho_core = new double[nrxx]; // core charge in real space
    ModuleBase::GlobalFunc::ZEROS(rho_core, nrxx);

    this->rhog_core = new std::complex<double>[ngmc]; // reciprocal core charge
    ModuleBase::GlobalFunc::ZEROS(rhog_core, ngmc);

    return;
}

void Charge::renormalize_rho(void)
{
    const double sr = this->sum_rho();
//    std::cout<<"charge before normalized "<<sr<<std::endl;
    const double normalize_factor = nelec / sr;

	for(int is=0; is<nspin; is++)
	{
	    	for(int ir=0; ir<nrxx; ir++)
		{
			rho[is][ir] *= normalize_factor;
		}
	}

//	double tmpdb = this->sum_rho();
//   std::cout<<"charge after normalized "<<tmpdb<<std::endl;
    return;
}

double Charge::sum_rho(void) const
{
    double sum_rho = 0.0;
    int nspin0 = 1;
    if(nspin==2)
    {
	nspin0 = 2;
    }

	for(int is=0; is<nspin0; is++)
	{
		for(int ir=0; ir<nrxx; ir++)
		{
			sum_rho += this->rho[is][ir];
		}
	}

	// multiply the sum of charge density by a factor
    sum_rho *= GlobalC::ucell.omega / static_cast<double>( GlobalC::rhopw->nxyz );

	// mohan fixed bug 2010-01-18,
	// sum_rho may be smaller than 1, like Na bcc.
    if (sum_rho <= 0.1)
    {
		GlobalV::ofs_warning << " sum_rho=" << sum_rho << std::endl;
        ModuleBase::WARNING_QUIT("Charge::renormalize_rho","Can't find even an electron!");
    }

    return sum_rho;
}



void Use_FFT::ToRealSpace(int const &is, const ModuleBase::ComplexMatrix &vg, double *vr, ModulePW::PW_Basis* rho_basis)
{
}
void Use_FFT::ToRealSpace(const std::complex<double> *vg, double *vr, ModulePW::PW_Basis* rho_basis)
{
}
bool Occupy::use_gaussian_broadening = false;
bool Occupy::use_tetrahedron_method = false;
double Magnetism::get_nelup()
{
    return 0;
}
double Magnetism::get_neldw()
{
    return 0;
}

bool ModuleSymmetry::Symmetry_Basic::equal(double const &m, double const &n) const
{
    return false;
}

void UnitCell_pseudo::setup_cell(
#ifdef __LCAO
    LCAO_Orbitals &orb,
#endif
    const std::string &s_pseudopot_dir,
    const std::string &fn,
    std::ofstream &log)
{
    // (1) init mag
    assert(ntype > 0);
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
    this->omega = abs(latvec.Det()) * this->lat0 * lat0 * lat0;

    //==========================================================
    // Calculate recip. lattice vectors and dot products
    // latvec have the unit of lat0, but G has the unit 2Pi/lat0
    //==========================================================
    this->GT = latvec.Inverse();
    this->G = GT.Transpose();
    this->GGT = G * GT;
    this->invGGT = GGT.Inverse();

    this->read_cell_pseudopots(s_pseudopot_dir,log);

#ifdef __LCAO
    // setup GlobalV::NLOCAL
    this->cal_nwfc(log);
#endif

    return;
}

int UnitCell_pseudo::read_atom_species(LCAO_Orbitals &orb, std::ifstream &ifa, std::ofstream &ofs_running)
{
    delete[] atom_label;
    delete[] atom_mass;
    delete[] pseudo_fn;
    this->atom_mass = new double[ntype]; // atom masses
    this->atom_label = new std::string[ntype]; // atom labels
    this->pseudo_fn = new std::string[ntype]; // file name of pseudopotential
    if (ModuleBase::GlobalFunc::SCAN_BEGIN(ifa, "ATOMIC_SPECIES"))
    {
        ModuleBase::GlobalFunc::OUT(ofs_running, "ntype", ntype);
        for (int i = 0; i < ntype; i++)
        {
            ifa >> this->atom_label[i] >> this->atom_mass[i];
	    ModuleBase::GlobalFunc::READ_VALUE(ifa, pseudo_fn[i]);
        }
    }
#ifdef __LCAO
    if (GlobalV::BASIS_TYPE == "lcao")
    {
        if (ModuleBase::GlobalFunc::SCAN_BEGIN(ifa, "NUMERICAL_ORBITAL"))
        {
            orb.read_in_flag = true;
            for (int i = 0; i < ntype; i++)
            {
                std::string ofile;
                ifa >> ofile;
                ofile = GlobalV::global_orbital_dir + ofile;
                orb.orbital_file.push_back(ofile);
            }
        }
    }
#endif
    if (ModuleBase::GlobalFunc::SCAN_BEGIN(ifa, "LATTICE_CONSTANT"))
    {
        ModuleBase::GlobalFunc::READ_VALUE(ifa, lat0);
    }
    lat0_angstrom = lat0 * 0.529177;
    this->tpiba = ModuleBase::TWO_PI / lat0;
    this->tpiba2 = tpiba * tpiba;
    if (latName == "sc")
    { // simple-cubic
        latvec.e11 = 1.0;
        latvec.e12 = 0.0;
        latvec.e13 = 0.0;
        latvec.e21 = 0.0;
        latvec.e22 = 1.0;
        latvec.e23 = 0.0;
        latvec.e31 = 0.0;
        latvec.e32 = 0.0;
        latvec.e33 = 1.0;
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

bool UnitCell_pseudo::read_atom_positions(LCAO_Orbitals &orb,
                                          std::ifstream &ifpos,
                                          std::ofstream &ofs_running,
                                          std::ofstream &ofs_warning)
{
    if (ModuleBase::GlobalFunc::SCAN_BEGIN(ifpos, "ATOMIC_POSITIONS"))
    {
        ModuleBase::GlobalFunc::READ_VALUE(ifpos, Coordinate);

        ModuleBase::Vector3<double> v;
        ModuleBase::Vector3<int> mv;
        int na = 0;
        this->nat = 0;
        assert(ntype > 0);
        for (int it = 0; it < ntype; it++)
        {
            //=======================================
            // (1) read in atom label
            // start magnetization
            //=======================================
            ModuleBase::GlobalFunc::READ_VALUE(ifpos, atoms[it].label);
            bool found = false;
            for (int it2 = 0; it2 < ntype; it2++)
            {
                if (this->atoms[it].label == this->atom_label[it])
                {
                    found = true;
                }
            }
            double magnet;
            ModuleBase::GlobalFunc::READ_VALUE(ifpos, magnet);
            //=========================
            // (3) read in atom number
            //=========================
            ModuleBase::GlobalFunc::READ_VALUE(ifpos, na);
            this->atoms[it].na = na;
            this->nat += na;
            if (na > 0)
            {
                atoms[it].tau = new ModuleBase::Vector3<double>[na];
				atoms[it].tau_original = new ModuleBase::Vector3<double>[na];
                atoms[it].taud = new ModuleBase::Vector3<double>[na];
                atoms[it].vel = new ModuleBase::Vector3<double>[na];
                for (int ia = 0; ia < na; ia++)
                {
                    ifpos >> v.x >> v.y >> v.z;
                    ifpos >> mv.x >> mv.y >> mv.z;
                    string tmpid;
                    tmpid = ifpos.get();
                    while (tmpid != "\n" && (ifpos.eof() == false))
                    {
                        tmpid = ifpos.get();
                    }
                    atoms[it].vel[ia].set(0, 0, 0);
                    if (Coordinate == "Cartesian")
                    {
                        atoms[it].tau[ia] = v;
                        // from Cartesian to Direct
                        double dx, dy, dz;
                        ModuleBase::Mathzone::Cartesian_to_Direct(atoms[it].tau[ia].x,
                                                                  atoms[it].tau[ia].y,
                                                                  atoms[it].tau[ia].z,
                                                                  latvec.e11,
                                                                  latvec.e12,
                                                                  latvec.e13,
                                                                  latvec.e21,
                                                                  latvec.e22,
                                                                  latvec.e23,
                                                                  latvec.e31,
                                                                  latvec.e32,
                                                                  latvec.e33,
                                                                  dx,
                                                                  dy,
                                                                  dz);

                        atoms[it].taud[ia].x = dx;
                        atoms[it].taud[ia].y = dy;
                        atoms[it].taud[ia].z = dz;
                        /*
                        std::cout<<"ia "<<ia <<" "
                            <<atoms[it].taud[ia].x
                            <<atoms[it].taud[ia].y
                            << atoms[it].taud[ia].z<<std::endl;
                        */
                    }
                }
            }
        }
    }

    // move negative coordinates
    this->check_dtau();
#ifdef __LCAO
    if (GlobalV::BASIS_TYPE == "lcao")
    {
        for (int it = 0; it < ntype; it++)
        {
            std::ifstream ifs(orb.orbital_file[it].c_str(), ios::in);
            char word[80];
            this->atoms[it].nw = 0;
            int L = 0;
            while (ifs.good())
            {
                ifs >> word;
                if (strcmp("Lmax", word) == 0)
                {
                    ModuleBase::GlobalFunc::READ_VALUE(ifs, this->atoms[it].nwl);
                    delete[] this->atoms[it].l_nchi;
                    this->atoms[it].l_nchi = new int[this->atoms[it].nwl + 1];
                }
                assert(this->atoms[it].nwl < 10);

                if (strcmp("Cutoff(a.u.)", word) == 0)
                {
                    ModuleBase::GlobalFunc::READ_VALUE(ifs, this->atoms[it].Rcut);
                }

                if (strcmp("Sorbital-->", word) == 0)
                {
                    ModuleBase::GlobalFunc::READ_VALUE(ifs, this->atoms[it].l_nchi[L]);
                    this->atoms[it].nw += (2 * L + 1) * this->atoms[it].l_nchi[L];
                    std::stringstream ss;
                    ss << "L=" << L << ", number of zeta";
                    ModuleBase::GlobalFunc::OUT(ofs_running, ss.str(), atoms[it].l_nchi[L]);
                    L++;
                }
                if (strcmp("Porbital-->", word) == 0)
                {
                    ModuleBase::GlobalFunc::READ_VALUE(ifs, this->atoms[it].l_nchi[L]);
                    this->atoms[it].nw += (2 * L + 1) * this->atoms[it].l_nchi[L];
                    std::stringstream ss;
                    ss << "L=" << L << ", number of zeta";
                    ModuleBase::GlobalFunc::OUT(ofs_running, ss.str(), atoms[it].l_nchi[L]);
                    L++;
                }
                if (strcmp("Dorbital-->", word) == 0)
                {
                    ModuleBase::GlobalFunc::READ_VALUE(ifs, this->atoms[it].l_nchi[L]);
                    this->atoms[it].nw += (2 * L + 1) * this->atoms[it].l_nchi[L];
                    std::stringstream ss;
                    ss << "L=" << L << ", number of zeta";
                    ModuleBase::GlobalFunc::OUT(ofs_running, ss.str(), atoms[it].l_nchi[L]);
                    L++;
                }
            }
            ifs.close();
        }
    }
#endif
    return true;
}

void UnitCell_pseudo::check_dtau(void)
{
    for (int it = 0; it < ntype; it++)
    {
        Atom *atom1 = &atoms[it];
        for (int ia = 0; ia < atoms[it].na; ia++)
        {
            double dx2 = (atom1->taud[ia].x + 10000) - int(atom1->taud[ia].x + 10000);
            double dy2 = (atom1->taud[ia].y + 10000) - int(atom1->taud[ia].y + 10000);
            double dz2 = (atom1->taud[ia].z + 10000) - int(atom1->taud[ia].z + 10000);

            // mohan add 2011-04-07
            while (dx2 >= 1)
            {
                GlobalV::ofs_warning << " dx2 is >=1 " << std::endl;
                dx2 -= 1.0;
            }
            while (dy2 >= 1)
            {
                GlobalV::ofs_warning << " dy2 is >=1 " << std::endl;
                dy2 -= 1.0;
            }
            while (dz2 >= 1)
            {
                GlobalV::ofs_warning << " dz2 is >=1 " << std::endl;
                dz2 -= 1.0;
            }
            // mohan add 2011-04-07
            while (dx2 < 0)
            {
                GlobalV::ofs_warning << " dx2 is <0 " << std::endl;
                dx2 += 1.0;
            }
            while (dy2 < 0)
            {
                GlobalV::ofs_warning << " dy2 is <0 " << std::endl;
                dy2 += 1.0;
            }
            while (dz2 < 0)
            {
                GlobalV::ofs_warning << " dz2 is <0 " << std::endl;
                dz2 += 1.0;
            }

            atom1->taud[ia].x = dx2;
            atom1->taud[ia].y = dy2;
            atom1->taud[ia].z = dz2;

            double cx2, cy2, cz2;

            ModuleBase::Mathzone::Direct_to_Cartesian(atom1->taud[ia].x,
                                                      atom1->taud[ia].y,
                                                      atom1->taud[ia].z,
                                                      latvec.e11,
                                                      latvec.e12,
                                                      latvec.e13,
                                                      latvec.e21,
                                                      latvec.e22,
                                                      latvec.e23,
                                                      latvec.e31,
                                                      latvec.e32,
                                                      latvec.e33,
                                                      cx2,
                                                      cy2,
                                                      cz2);

            atom1->tau[ia].x = cx2;
            atom1->tau[ia].y = cy2;
            atom1->tau[ia].z = cz2;
        }
    }
    return;
}

//#ifdef _LCAO

void Run_lcao::lcao_line(ModuleESolver::ESolver *p_esolver)
{
    GlobalC::sf.set(INPUT.nbspline);
    GlobalC::ucell.setup_cell(GlobalC::ORB, GlobalV::global_pseudo_dir, GlobalV::stru_file, GlobalV::ofs_running);
    GlobalC::kv.set(GlobalC::symm,
                    GlobalV::global_kpoint_card,
                    GlobalV::NSPIN,
                    GlobalC::ucell.G,
                    GlobalC::ucell.latvec);

	// pw_rho = new ModuleBase::PW_Basis();
    //temporary, it will be removed
    GlobalC::rhopw = new ModulePW::PW_Basis_Big();
    GlobalC::bigpw = static_cast<ModulePW::PW_Basis_Big*>(GlobalC::rhopw);
	GlobalC::bigpw->setbxyz(INPUT.bx,INPUT.by,INPUT.bz);
    GlobalC::wfcpw = new ModulePW::PW_Basis_K_Big(); 
    ModulePW::PW_Basis_K_Big* tmp2 = static_cast<ModulePW::PW_Basis_K_Big*>(GlobalC::wfcpw);
    tmp2->setbxyz(INPUT.bx,INPUT.by,INPUT.bz);
#ifdef __MPI
    GlobalC::rhopw->initmpi(1, 0 ,POOL_WORLD);
    GlobalC::wfcpw->initmpi(1, 0 ,POOL_WORLD);
#endif
    GlobalC::rhopw->initgrids(GlobalC::ucell.lat0, GlobalC::ucell.latvec, 4 * INPUT.ecutwfc);
    GlobalC::rhopw->initparameters(false, 4 * INPUT.ecutwfc);
	GlobalC::rhopw->setuptransform();
	GlobalC::wfcpw->initgrids(GlobalC::ucell.lat0, GlobalC::ucell.latvec, GlobalC::rhopw->nx, GlobalC::rhopw->ny, GlobalC::rhopw->nz);
    GlobalC::wfcpw->initparameters(false, INPUT.ecutwfc, GlobalC::kv.nks, GlobalC::kv.kvec_d.data());
	
    GlobalC::CHR.allocate(GlobalV::NSPIN, GlobalC::rhopw->nrxx, GlobalC::rhopw->npw);

    ORB_control orb_con(GlobalV::GAMMA_ONLY_LOCAL,
                        GlobalV::NLOCAL,
                        GlobalV::NBANDS,
                        GlobalV::NSPIN,
                        GlobalV::DSIZE,
                        GlobalV::NB2D,
                        GlobalV::DCOLOR,
                        GlobalV::DRANK,
                        GlobalV::MY_RANK,
                        GlobalV::CALCULATION,
                        GlobalV::KS_SOLVER);

    // void Run_lcao::Init_Basis_lcao
    orb_con.read_orb_first(GlobalV::ofs_running,
                           GlobalC::ORB,
                           GlobalC::ucell.ntype,
                           GlobalC::ucell.lmax,
                           INPUT.lcao_ecut,
                           INPUT.lcao_dk,
                           INPUT.lcao_dr,
                           INPUT.lcao_rmax,
                           GlobalV::deepks_setorb,
                           INPUT.out_mat_r,
                           GlobalV::CAL_FORCE,
                           GlobalV::MY_RANK);

    GlobalC::ucell.infoNL.setupNonlocal(GlobalC::ucell.ntype, GlobalC::ucell.atoms, GlobalV::ofs_running, GlobalC::ORB);

    int Lmax = 0;
    orb_con.set_orb_tables(GlobalV::ofs_running,
                           GlobalC::UOT,
                           GlobalC::ORB,
                           GlobalC::ucell.lat0,
                           GlobalV::deepks_setorb,
                           Lmax,
                           GlobalC::ucell.infoNL.nprojmax,
                           GlobalC::ucell.infoNL.nproj,
                           GlobalC::ucell.infoNL.Beta);
}

void UnitCell_pseudo::cal_nwfc(std::ofstream &log)
{
    assert(ntype > 0);
    assert(nat > 0);

    //===========================
    // (1) set iw2l, iw2n, iw2m
    //===========================
    for (int it = 0; it < ntype; it++)
    {
        this->atoms[it].set_index();
    }

    //===========================
    // (2) set namax and nwmax
    //===========================
    this->namax = 0;
    this->nwmax = 0;
    for (int it = 0; it < ntype; it++)
    {
        this->namax = std::max(atoms[it].na, namax);
        this->nwmax = std::max(atoms[it].nw, nwmax);
    }
    assert(namax > 0);
    //===========================
    // (3) set nwfc and stapos_wf
    //===========================
    GlobalV::NLOCAL = 0;
    for (int it = 0; it < ntype; it++)
    {
        // std::cout<<" it "<<it<<" na "<<atoms[it].na<<" nw "<<atoms[it].nw<<std::endl;
        atoms[it].stapos_wf = GlobalV::NLOCAL;
        const int nlocal_it = atoms[it].nw * atoms[it].na;
        GlobalV::NLOCAL += nlocal_it;
    }
    this->set_iat2itia(); // used in gint_tools::get_block_info()
    assert(GlobalV::NLOCAL > 0);
    delete[] iwt2iat;
    delete[] iwt2iw;
    this->iwt2iat = new int[GlobalV::NLOCAL];
    this->iwt2iw = new int[GlobalV::NLOCAL];

    this->itia2iat.create(ntype, namax);
    this->itiaiw2iwt.create(ntype, namax, nwmax * GlobalV::NPOL);
    int iat = 0;
    int iwt = 0;
    for (int it = 0; it < ntype; it++)
    {
        for (int ia = 0; ia < atoms[it].na; ia++)
        {
            this->itia2iat(it, ia) = iat;
            // this->iat2ia[iat] = ia;
            for (int iw = 0; iw < atoms[it].nw * GlobalV::NPOL; iw++)
            {
                this->itiaiw2iwt(it, ia, iw) = iwt;
                this->iwt2iat[iwt] = iat;
                this->iwt2iw[iwt] = iw;
                ++iwt;
            }
            ++iat;
        }
    }

    //========================
    // (5) set lmax and nmax
    //========================
    this->lmax = 0;
    this->nmax = 0;
    for (int it = 0; it < ntype; it++)
    {
        lmax = std::max(lmax, atoms[it].nwl);
        for (int l = 0; l < atoms[it].nwl + 1; l++)
        {
            nmax = std::max(nmax, atoms[it].l_nchi[l]);
        }

        int nchi = 0;
        for (int l = 0; l < atoms[it].nwl + 1; l++)
        {
            nchi += atoms[it].l_nchi[l];
        }
        this->nmax_total = std::max(nmax_total, nchi);
    }

    //=======================
    // (6) set lmax_ppwf
    //=======================
    this->lmax_ppwf = 0;
    for (int it = 0; it < ntype; it++)
    {
        for (int ic = 0; ic < atoms[it].nchi; ic++)
        {
            if (lmax_ppwf < atoms[it].lchi[ic])
            {
                this->lmax_ppwf = atoms[it].lchi[ic];
            }
        }
    }
}

void UnitCell::set_iat2itia(void)
{
    assert(nat > 0);
    delete[] iat2it;
    delete[] iat2ia;
    this->iat2it = new int[nat];
    this->iat2ia = new int[nat];
    int iat = 0;
    for (int it = 0; it < ntype; it++)
    {
        for (int ia = 0; ia < atoms[it].na; ia++)
        {
            this->iat2it[iat] = it;
            this->iat2ia[iat] = ia;
            ++iat;
        }
    }
    return;
}

/******************************
 * Read LOWF_GAMMA_S1.dat
 * mock void Local_Orbital_Charge::gamma_file
 * because of lack of codes without mpi
 * in WF_Local::distri_lowf_new() called by WF_Local::read_lowf()
 ******************************/
void Local_Orbital_Charge::gamma_file(psi::Psi<double>* psid, Local_Orbital_wfc &lowf)
{
    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        lowf.wfc_gamma[is].create(this->ParaV->ncol, this->ParaV->nrow);
        lowf.wfc_gamma[is].zero_out();
        std::stringstream ss;
        ss << GlobalV::global_out_dir << "LOWF_GAMMA_S" << is + 1;
        // std::cout << " name is = " << ss.str() << std::endl;
        std::ifstream ifs;
        if (GlobalV::DRANK == 0)
        {
            ifs.open(ss.str().c_str());
            if (!ifs)
            {
                GlobalV::ofs_warning << " Can't open file:" << ss.str() << std::endl;
            }
        }
        if (GlobalV::MY_RANK == 0)
        {
            int nbands, nlocal;
            ModuleBase::GlobalFunc::READ_VALUE(ifs, nbands);
            ModuleBase::GlobalFunc::READ_VALUE(ifs, nlocal);

            if (nbands != GlobalV::NBANDS)
            {
                GlobalV::ofs_warning << " read in nbands=" << nbands;
                GlobalV::ofs_warning << " NBANDS=" << GlobalV::NBANDS << std::endl;
            }
            else if (nlocal != GlobalV::NLOCAL)
            {
                GlobalV::ofs_warning << " read in nlocal=" << nlocal;
                GlobalV::ofs_warning << " NLOCAL=" << GlobalV::NLOCAL << std::endl;
            }
        }
        for (int i = 0; i < GlobalV::NBANDS; i++)
        {
            int ib;
            ModuleBase::GlobalFunc::READ_VALUE(ifs, ib);
            ModuleBase::GlobalFunc::READ_VALUE(ifs, GlobalC::wf.ekb[GlobalV::CURRENT_SPIN][i]);
            ModuleBase::GlobalFunc::READ_VALUE(ifs, GlobalC::wf.wg(GlobalV::CURRENT_SPIN, i));
            assert((i + 1) == ib);
            for (int j = 0; j < GlobalV::NLOCAL; j++)
            {
                ifs >> lowf.wfc_gamma[is](i, j);
                // std::cout <<lowf.wfc_gamma[is](i,j)  << " ";
            }
            // std::cout << std::endl;
        }
    }
}

void Local_Orbital_Charge::allocate_gamma(const int& lgd, psi::Psi<double>* psid)
{
    this->DM = new double **[GlobalV::NSPIN];
    this->DM_pool = new double *[GlobalV::NSPIN];
    for (int is = 0; is < GlobalV::NSPIN; is++)
    {
        this->DM_pool[is] = new double[lgd * lgd];
        ModuleBase::GlobalFunc::ZEROS(DM_pool[is], lgd * lgd);
        this->DM[is] = new double *[lgd];
        for (int i = 0; i < lgd; i++)
        {
            DM[is][i] = &DM_pool[is][i * lgd];
        }
    }
}

#include "../../module_base/blas_connector.h"
//calculate density matrix from wfc_gamma in gamma_only case
void Local_Orbital_Charge::cal_dm(const ModuleBase::matrix& wg,
    std::vector<ModuleBase::matrix>& wfc_gamma,
    std::vector<ModuleBase::matrix>& dm_gamma)
{
	
	// dm = wfc.T * wg * wfc.conj()
	// dm[is](iw1,iw2) = \sum_{ib} wfc[is](ib,iw1).T * wg(is,ib) * wfc[is](ib,iw2).conj()
	assert(wg.nc<=GlobalV::NLOCAL);
    assert(wg.nr==GlobalV::NSPIN);
    for(int is=0; is!=GlobalV::NSPIN; ++is)
    {
        std::vector<double> wg_local(this->ParaV->ncol,0.0);
        for(int ib_global=0; ib_global!=wg.nc; ++ib_global)
        {
            const int ib_local = this->ParaV->trace_loc_col[ib_global];
            if(ib_local>=0)
            {
                wg_local[ib_local] = wg(is,ib_global);
            }
        }
        
        // wg_wfc(ib,iw) = wg[ib] * wfc(ib,iw);
        ModuleBase::matrix wg_wfc(wfc_gamma[is]);
        for(int ir=0; ir!=wg_wfc.nr; ++ir)
        {
            BlasConnector::scal( wg_wfc.nc, wg_local[ir], wg_wfc.c+ir*wg_wfc.nc, 1 );
        }

        // C++: dm(iw1,iw2) = wfc(ib,iw1).T * wg_wfc(ib,iw2)
        const double one_float=1.0, zero_float=0.0;
        const int one_int=1;
        const char N_char='N', T_char='T';
        dm_gamma[is].create( this->ParaV->ncol, this->ParaV->nrow );
        const int lda=GlobalV::NLOCAL;
        dgemm_(
            &N_char, &T_char, 
            &GlobalV::NLOCAL, &GlobalV::NLOCAL, &wg.nc,
            &one_float,
            wg_wfc.c, &lda,
            wfc_gamma[is].c, &lda,
            &zero_float,
            dm_gamma[is].c, &lda);
    }

	return;
}

void Local_Orbital_Charge::read_dm(const int &is, const std::string &fn)
{
    std::ifstream ifs;
    ifs.open(fn.c_str());
    std::string name;
    ifs >> name;

    bool quit = false;

    // check lattice constant, unit is Angstrom
    ModuleBase::CHECK_DOUBLE(ifs, GlobalC::ucell.lat0 * ModuleBase::BOHR_TO_A, quit);
    ModuleBase::CHECK_DOUBLE(ifs, GlobalC::ucell.latvec.e11, quit);
    ModuleBase::CHECK_DOUBLE(ifs, GlobalC::ucell.latvec.e12, quit);
    ModuleBase::CHECK_DOUBLE(ifs, GlobalC::ucell.latvec.e13, quit);
    ModuleBase::CHECK_DOUBLE(ifs, GlobalC::ucell.latvec.e21, quit);
    ModuleBase::CHECK_DOUBLE(ifs, GlobalC::ucell.latvec.e22, quit);
    ModuleBase::CHECK_DOUBLE(ifs, GlobalC::ucell.latvec.e23, quit);
    ModuleBase::CHECK_DOUBLE(ifs, GlobalC::ucell.latvec.e31, quit);
    ModuleBase::CHECK_DOUBLE(ifs, GlobalC::ucell.latvec.e32, quit);
    ModuleBase::CHECK_DOUBLE(ifs, GlobalC::ucell.latvec.e33, quit);

    for (int it = 0; it < GlobalC::ucell.ntype; it++)
    {
        ModuleBase::CHECK_STRING(ifs, GlobalC::ucell.atoms[it].label, quit);
    }

    for (int it = 0; it < GlobalC::ucell.ntype; it++)
    {
        ModuleBase::CHECK_DOUBLE(ifs, GlobalC::ucell.atoms[it].na, quit);
    }
    std::string coordinate;
    ifs >> coordinate;
    double tau;

    for (int it = 0; it < GlobalC::ucell.ntype; it++)
    {
        for (int ia = 0; ia < GlobalC::ucell.atoms[it].na; ia++)
        {
            ModuleBase::CHECK_DOUBLE(ifs, GlobalC::ucell.atoms[it].taud[ia].x, quit);
            ModuleBase::CHECK_DOUBLE(ifs, GlobalC::ucell.atoms[it].taud[ia].y, quit);
            ModuleBase::CHECK_DOUBLE(ifs, GlobalC::ucell.atoms[it].taud[ia].z, quit);
        }
    }
    ModuleBase::CHECK_INT(ifs, GlobalV::NSPIN);
    if (GlobalV::NSPIN == 1 || GlobalV::NSPIN == 4)
    {
        ModuleBase::GlobalFunc::READ_VALUE(ifs, GlobalC::en.ef);
    }
    ModuleBase::CHECK_INT(ifs, GlobalV::NLOCAL);
    ModuleBase::CHECK_INT(ifs, GlobalV::NLOCAL);
    if (GlobalV::GAMMA_ONLY_LOCAL)
    {
        for (int i = 0; i < GlobalV::NLOCAL; ++i)
        {
            for (int j = 0; j < GlobalV::NLOCAL; ++j)
            {
                ifs >> DM[is][i][j];
                // std::cout<<DM[is][i][j]<<" ";
            }
            // std::cout<<endl;
        }
    }
    ifs.close();
}

void set_matrix_grid()
{
    GlobalV::SEARCH_RADIUS = atom_arrange::set_sr_NL(GlobalV::ofs_running,
                                                     GlobalV::OUT_LEVEL,
                                                     GlobalC::ORB.get_rcutmax_Phi(),
                                                     GlobalC::ucell.infoNL.get_rcutmax_Beta(),
                                                     GlobalV::GAMMA_ONLY_LOCAL);
    atom_arrange::search(GlobalV::SEARCH_PBC,
                         GlobalV::ofs_running,
                         GlobalC::GridD,
                         GlobalC::ucell,
                         GlobalV::SEARCH_RADIUS,
                         GlobalV::test_atom_input);
    GlobalC::GridT.set_pbc_grid(GlobalC::rhopw->nx,
                                GlobalC::rhopw->ny,
                                GlobalC::rhopw->nz,
                                GlobalC::bigpw->bx,
                                GlobalC::bigpw->by,
                                GlobalC::bigpw->bz,
                                GlobalC::bigpw->nbx,
                                GlobalC::bigpw->nby,
                                GlobalC::bigpw->nbz,
                                GlobalC::bigpw->nbxx,
                                GlobalC::bigpw->nbzp_start,
                                GlobalC::bigpw->nbzp);
    // std::cout << "GridT.max_atom in set_matrix_grid " << GlobalC::GridT.max_atom << std::endl;
}
//#endif
//====================================================
//mock of  module_orbital/ORB_control.cpp
#include "module_orbital/ORB_control.h"
#include "module_orbital/ORB_gen_tables.h"
#include "module_base/timer.h"
#include "src_io/wf_local.h"
#include "module_base/lapack_connector.h"
#include "module_base/memory.h"

ORB_control::ORB_control(
    const bool& gamma_only_in,
    const int& nlocal_in,
    const int& nbands_in,
    const int& nspin_in,
    const int& dsize_in,
    const int& nb2d_in,
    const int& dcolor_in,
    const int& drank_in,
    const int& myrank_in,
    const std::string& calculation_in,
    const std::string& ks_solver_in) :
    gamma_only(gamma_only_in),
    nlocal(nlocal_in),
    nbands(nbands_in),
    nspin(nspin_in),
    dsize(dsize_in),
    nb2d(nb2d_in),
    dcolor(dcolor_in),
    drank(drank_in),
    myrank(myrank_in),
    calculation(calculation_in),
    ks_solver(ks_solver_in),
    setup_2d(true)
{
    this->ParaV.nspin = nspin_in;
}

ORB_control::ORB_control() :
    setup_2d(false)
{}
ORB_control::~ORB_control()
{}

void ORB_control::read_orb_first(
	std::ofstream &ofs_in, 
	LCAO_Orbitals &orb,
	const int &ntype, // mohan add 2021-04-26
	const int &lmax, // mohan add 2021-04-26 
	const double &lcao_ecut_in, // mohan add 2021-04-16
	const double &lcao_dk_in, // mohan add 2021-04-16
	const double &lcao_dr_in, // mohan add 2021-04-16
	const double &lcao_rmax_in, // mohan add 2021-04-16
	const bool &deepks_setorb,
	const int &out_mat_r,
	const bool &force_flag, // mohan add 2021-05-07
	const int &my_rank // mohan add 2021-04-26
) 
{
    ModuleBase::TITLE("ORB_control","read_orb_first");
	ModuleBase::timer::tick("ORB_control","read_orb_first");
    
	/////////////////////////////////////////////////////////////////
	/// (1) FUNCTION : use 'info' to generate 'Numerical Orbital'
	///
	/// (1) RESULT : We have 'Numerical Orbital' for calculate S-table and T-table.
	/////////////////////////////////////////////////////////////////

	// mohan add 2021-04-16
	assert(ntype>0);
	assert(lmax>=0);
	assert(lcao_ecut_in>0.0);
	assert(lcao_dk_in>0.0);
	assert(lcao_dr_in>0.0);
	assert(lcao_rmax_in>0.0);

	// mohan add 2021-04-16
	orb.ecutwfc = lcao_ecut_in;
	orb.dk = lcao_dk_in;
	orb.dR = lcao_dr_in;
	orb.Rmax = lcao_rmax_in;
	
    orb.Read_Orbitals(
		ofs_in,
		ntype, 
		lmax, 
		deepks_setorb, 
		out_mat_r, 
		force_flag,
		my_rank);

	ModuleBase::timer::tick("ORB_control","read_orb_first");
	return;
}

void ORB_control::set_orb_tables(
	std::ofstream &ofs_in,
	ORB_gen_tables &OGT, 
	LCAO_Orbitals &orb,
	const double &lat0,
	const bool &deepks_setorb,
	const int &Lmax_exx,
	const int &nprojmax, 
	const int* nproj,
	const Numerical_Nonlocal* beta_) 
{
    ModuleBase::TITLE("ORB_control","set_orb_tables");
	ModuleBase::timer::tick("ORB_control","set_orb_tables");


#ifdef __NORMAL

#else
	if(calculation=="test")
	{
		ModuleBase::timer::tick("ORB_control","set_orb_tables");
		return;
	}
#endif


    ///////////////////////////////////////////////////////////////////
    /// (2) FUNCTION : Generate Gaunt_Coefficients and S-table using OGT.init
	/// 	   Must have 'Numerical Orbital' infomation
	///
	/// (2) RESULT : we have tabulated S table for use.
    ///////////////////////////////////////////////////////////////////
    const int job0 = 3;
    /// job0 :
    /// 1. generate overlap table
    /// 2. generate kinetic table
    /// 3. generate overlap & kinetic table
    OGT.gen_tables(ofs_in, job0, orb, Lmax_exx, deepks_setorb, nprojmax, nproj, beta_);
    // init lat0, in order to interpolated value from this table.

	assert(lat0>0.0);
    OGT.set_unit(lat0);

	ModuleBase::timer::tick("ORB_control","set_orb_tables");
    return;
}
//mock of  module_orbital/ORB_control.cpp
//====================================================
//====================================================
//mock of  module_orbital/ORB_gen_tables.cpp
#include "module_orbital/ORB_read.h"
#include "module_orbital/ORB_gen_tables.h"
#include "module_base/ylm.h"
#include "module_base/math_polyint.h"
#include "module_base/timer.h"


ORB_gen_tables::ORB_gen_tables() {}
ORB_gen_tables::~ORB_gen_tables() {}

/// call in hamilt_linear::init_before_ions.
void ORB_gen_tables::gen_tables(
	std::ofstream &ofs_in,
	const int &job0,
	LCAO_Orbitals &orb,
	const int &Lmax_exx,
	const bool &deepks_setorb,
	const int &nprojmax, 
	const int* nproj,
	const Numerical_Nonlocal* beta_)
{
	ModuleBase::TITLE("ORB_gen_tables", "gen_tables");
	ModuleBase::timer::tick("ORB_gen_tables", "gen_tables");

	ofs_in << "\n SETUP THE TWO-CENTER INTEGRATION TABLES" << std::endl;

	//////////////////////////////
	/// (1) MOT: make overlap table.
	//////////////////////////////
	MOT.allocate(
		orb.get_ntype(),
		orb.get_lmax(),	 
		orb.get_kmesh(), 
		orb.get_Rmax(),	
		orb.get_dR(),	 
		orb.get_dk());	 

	tbeta.allocate(
		orb.get_ntype(),
		orb.get_lmax(),	
		orb.get_kmesh(), 
		orb.get_Rmax(),	
		orb.get_dR(),
		orb.get_dk());

	//caoyu add 2021-03-18
	//mohan update 2021-04-22
	if (deepks_setorb)
	{
		talpha.allocate(
			orb.get_ntype(), 
			orb.get_lmax(),	
			orb.get_kmesh(),
			orb.get_Rmax(),
			orb.get_dR(),
			orb.get_dk());
	}

	// OV: overlap
	MOT.init_OV_Tpair(orb);
	MOT.init_OV_Opair(orb);

	// NL: nonlocal
	tbeta.init_NL_Tpair(orb.Phi, beta_);
	tbeta.init_NL_Opair(orb, nprojmax, nproj); // add 2009-5-8

	//caoyu add 2021-03-18
	// DS: Descriptor
	if (deepks_setorb)
	{
		talpha.init_DS_Opair(orb);
		talpha.init_DS_2Lplus1(orb);
	}

	//////////////////////////////
	/// (2) init Ylm Coef
	//////////////////////////////
	//liaochen add 2010/4/29
	ModuleBase::Ylm::set_coefficients();

	// PLEASE add explanations for all options of 'orb_num' and 'mode'
	// mohan add 2021-04-03
	// Peize Lin update 2016-01-26
#ifdef __ORBITAL
	int orb_num = 4;
#else
	int orb_num = 2; //
#endif
	int mode = 1;	 // 1: <phi|phi> and <phi|beta>
	int Lmax_used = 0;
	int Lmax = 0;


	MOT.init_Table_Spherical_Bessel(orb_num, mode, Lmax_used, Lmax, Lmax_exx, orb, beta_);

	//calculate S(R) for interpolation
	MOT.init_Table(job0, orb);
	tbeta.init_Table_Beta(MOT.pSB, orb.Phi, beta_, nproj); // add 2009-5-8

	//caoyu add 2021-03-18
	if (deepks_setorb)
	{
		talpha.init_Table_Alpha(MOT.pSB, orb);
		if(GlobalV::deepks_out_unittest) talpha.print_Table_DSR(orb);
	}

	/////////////////////////////
	/// (3) make Gaunt coefficients table
	/////////////////////////////

	const int lmax = (Lmax_used - 1) / 2;
	//MGT.init_Ylm_Gaunt(orb.get_lmax()+1, 0.0,PI,0.0,ModuleBase::TWO_PI);
	MGT.init_Gaunt_CH(lmax);
	//MGT.init_Gaunt(orb.get_lmax()+1);
	MGT.init_Gaunt(lmax);

	ModuleBase::timer::tick("ORB_gen_tables", "gen_tables");
	return;
}
//mock of  module_orbital/ORB_gen_tables.cpp
//====================================================
//====================================================
//mock of module_orbital/ORB_read.cpp
#include "module_orbital/ORB_read.h"
#include <cstring>		// Peize Lin fix bug about strcmp 2016-08-02
#include <cassert>
#include "module_base/math_integral.h"
#include "module_base/tool_check.h"
#include <algorithm>
#include "module_base/timer.h"

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
	if(out_mat_r) 
	{
		kmesh = kmesh * 4;
	}

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
	
	int* nchi = new int[lmax+1];		// allocate space: number of chi for each L.
	
	if (my_rank == 0)
	{	
		for (int l = 0; l <= lmax; l++)
		{
			ifs >> word >> word >> word >> nchi[l];
			nchimax = std::max(nchimax, nchi[l]);
		}
	}


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

			if (!find)
			{
				ModuleBase::WARNING_QUIT("LCAO_Orbitals::read_orb_file", "Can't find orbitals.");
			}


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
//mock of module_orbital/ORB_read.cpp
//====================================================
//====================================================
//mock of module_cell/setup_nonlocal.cpp
#include "../../module_cell/setup_nonlocal.h"
#include "../../module_base/complexmatrix.h"

#ifdef __LCAO
//#include "../src_pw/global.h"
//#include "../src_pw/soc.h"
// mohan add 2013-08-02
// In order to get rid of the read in file .NONLOCAL.

InfoNonlocal::InfoNonlocal()
{
    this->Beta = new Numerical_Nonlocal[1];
	this->nproj = nullptr;
    this->nprojmax = 0;
    this->rcutmax_Beta = 0.0;
}
InfoNonlocal::~InfoNonlocal()
{
    delete[] Beta;
	delete[] nproj;
}


void InfoNonlocal::setupNonlocal(
    const int& ntype,
    Atom* atoms,
    std::ofstream &log,
    LCAO_Orbitals &orb
)
{
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	//~~~~~~~~~~~~~~~~~~~~~~   2    ~~~~~~~~~~~~~~~~~~~~~~~~~
	// Read in non-local projector for each atom type.
	// In fact this should be improved,
	// the non-local projector should be transferred
	// from .UPF file directly.
	// mohan note 2011-03-04
	//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	if(GlobalV::BASIS_TYPE=="lcao" || GlobalV::BASIS_TYPE=="lcao_in_pw")
	{
		delete[] this->Beta;
		this->Beta = new Numerical_Nonlocal[ntype];

		delete[] this->nproj;
		this->nproj = new int[ntype];
		ModuleBase::GlobalFunc::ZEROS(this->nproj, ntype);
		
		this->nprojmax = 0;
		

		// if true: read in the nonlocal file from file.
		// if false: get nonlocal information from .upf or .vwr directly
		bool readin_nonlocal = false;

		for(int it=0; it<ntype; it++)
		{
			Atom* atom = &atoms[it];
			if(readin_nonlocal)
			{
				this->Read_NonLocal(
					it, 
					atom, 
					this->nproj[it], 
					GlobalV::MY_RANK, 
					orb.get_kmesh(),
					orb.get_dk(),
					orb.get_dr_uniform(),
					orb.orbital_file[it] );	
			}
			else
			{
				this->Set_NonLocal(
					it, 
					atom, 
					this->nproj[it],
					orb.get_kmesh(),
					orb.get_dk(),
					orb.get_dr_uniform(),
					log);
			}
			this->nprojmax = std::max(this->nprojmax, this->nproj[it]);
			//caoyu add 2021-05-24 to reconstruct atom_arrange::set_sr_NL
			this->rcutmax_Beta = std::max(this->rcutmax_Beta, this->Beta[it].get_rcut_max());
		}

		log << " max number of nonlocal projetors among all species is " << this->nprojmax << std::endl; 
	}
    return;
}

void InfoNonlocal::Set_NonLocal(
    const int &it, 
    Atom* atom, 
    int &n_projectors,
    const int& kmesh,
    const double& dk,
    const double& dr_uniform,
	std::ofstream &log)
{
	ModuleBase::TITLE("InfoNonlocal","Set_NonLocal");

	// set a pointer
	//Atom* atom = &GlobalC::ucell.atoms[it];

	// get the number of non-local projectors
	n_projectors = atom->nbeta;

	//const int nh = atom->nh;//zhengdy-soc

	// set the nonlocal projector objects
	//Numerical_Nonlocal_Lm* tmpBeta_lm = new Numerical_Nonlocal_Lm[n_projectors];
	Numerical_Nonlocal_Lm* tmpBeta_lm = new Numerical_Nonlocal_Lm[n_projectors];

	//ModuleBase::ComplexMatrix coefficient_D_nc_in(nh*2, nh*2);//zhengdy-soc

	if(!atom->has_so)
	{
		for(int p1 = 0; p1<n_projectors; p1++)
		{
			const int lnow = atom->lll[p1];

			// only keep the nonzero part.
			int cut_mesh = atom->mesh; 
			for(int ir=atom->mesh-1; ir>=0; --ir)
			{
				if( abs( atom->betar(p1,ir) ) > 1.0e-10 )
				{
					cut_mesh = ir; 
					break;
				}
			}
			if(cut_mesh %2 == 0) ++cut_mesh;

//		std::cout << " cut_mesh=" << cut_mesh << std::endl;
			double* beta_r = new double[cut_mesh];
			ModuleBase::GlobalFunc::ZEROS(beta_r, cut_mesh);
			for(int ir=0; ir<cut_mesh; ++ir)
			{
				beta_r[ir] = atom->betar(p1,ir);
			}

			tmpBeta_lm[p1].set_NL_proj(
					atom->label,
					it, //type
					lnow, // angular momentum L
					cut_mesh, // number of radial mesh
					atom->rab,
					atom->r, // radial mesh value (a.u.)
					beta_r,
					kmesh,
					dk,
					dr_uniform); // delta k mesh in reciprocal space

			if(GlobalV::out_element_info)tmpBeta_lm[p1].plot(GlobalV::MY_RANK);

			delete[] beta_r;
				
		}
		
		// mohan comment out 2021-04-26
		//ModuleBase::WARNING("InfoNonlocal::Set_NonLocal","bug in line "+TO_STRING(__LINE__)+", matrix ic>=nc");		


		// Peize Lin add 2019-01-23
		this->Beta[it].set_type_info(
			it, 
			atom->label, 
			atom->pp_type, 
			atom->lmax, 
			n_projectors, 
			tmpBeta_lm);//LiuXh 2016-01-14, 2016-07-19

		// mohan add 2021-05-07
		//atom->set_d_so(coefficient_D_nc_in,n_projectors,0,0);
	}

	delete[] tmpBeta_lm;

	log << " SET NONLOCAL PSEUDOPOTENTIAL PROJECTORS" << std::endl;
	return;
}

#endif
//mock of module_cell/setup_nonlocal.cpp
//====================================================
//====================================================
//mock of module_cell/read_cell_pseudopots.cpp
#include "module_cell/unitcell_pseudo.h"
#include <cstring>

//==========================================================
// Read pseudopotential according to the dir
//==========================================================
void UnitCell_pseudo::read_cell_pseudopots(const std::string &pp_dir, std::ofstream &log)
{
	ModuleBase::TITLE("UnitCell_pseudo","read_cell_pseudopots");
	// setup reading log for pseudopot_upf
	std::stringstream ss;
	ss << GlobalV::global_out_dir << "atom_pseudo.log";
	
	// Read in the atomic pseudo potentials
	std::string pp_address;
	for (int i = 0;i < ntype;i++)
	{
		Pseudopot_upf upf;
	
		// mohan update 2010-09-12	
		int error = 0;
		int error_ap = 0;
		
		if(GlobalV::MY_RANK==0)
		{
			pp_address = pp_dir + this->pseudo_fn[i];
			error = upf.init_pseudo_reader( pp_address ); //xiaohui add 2013-06-23

			if(error==0) // mohan add 2021-04-16
			{
				if(this->atoms[i].flag_empty_element)	// Peize Lin add for bsse 2021.04.07
				{
					upf.set_empty_element();			
				}
				//average pseudopotential if needed
				error_ap = upf.average_p(GlobalV::soc_lambda); //added by zhengdy 2020-10-20
			}
		}


		if(error_ap) 
		{
			ModuleBase::WARNING_QUIT("UnitCell_pseudo::read_pseudopot","error when average the pseudopotential.");
		}

		if(error==1)
		{
			std::cout << " Pseudopotential directory now is : " << pp_address << std::endl;
			GlobalV::ofs_warning << " Pseudopotential directory now is : " << pp_address << std::endl;
			ModuleBase::WARNING_QUIT("read_pseudopot","Couldn't find pseudopotential file.");
		}
		else if(error==2)
		{
			ModuleBase::WARNING_QUIT("read_pseudopot","Pseudopotential data do not match.");
		}
		else if(error==3)
		{
			ModuleBase::WARNING_QUIT("read_pseudopot","Check the reference states in pseudopotential .vwr file.\n Also the norm of the read in pseudo wave functions\n explicitly please check S, P and D channels.\n If the norm of the wave function is \n unreasonable large (should be near 1.0), ABACUS would quit. \n The solution is to turn off the wave functions  \n and the corresponding non-local projectors together\n in .vwr pseudopotential file.");
		}

//xiaohui add 2015-03-24
		//xiaohui add 2015-03-24

		if(GlobalV::MY_RANK==0)
		{
//			upf.print_pseudo_upf( ofs );
			atoms[i].set_pseudo_nc( upf );

			log << "\n Read in pseudopotential file is " << pseudo_fn[i] << std::endl;
			ModuleBase::GlobalFunc::OUT(log,"pseudopotential type",atoms[i].pp_type);
			ModuleBase::GlobalFunc::OUT(log,"exchange-correlation functional", atoms[i].xc_func);
			ModuleBase::GlobalFunc::OUT(log,"nonlocal core correction", atoms[i].nlcc);
//			ModuleBase::GlobalFunc::OUT(log,"spin orbital",atoms[i].has_so);
			ModuleBase::GlobalFunc::OUT(log,"valence electrons", atoms[i].zv);
			ModuleBase::GlobalFunc::OUT(log,"lmax", atoms[i].lmax);
			ModuleBase::GlobalFunc::OUT(log,"number of zeta", atoms[i].nchi);
			ModuleBase::GlobalFunc::OUT(log,"number of projectors", atoms[i].nbeta);
			for(int ib=0; ib<atoms[i].nbeta; ib++)
			{
				ModuleBase::GlobalFunc::OUT(log,"L of projector", atoms[i].lll[ib]);
			}
//			ModuleBase::GlobalFunc::OUT(log,"Grid Mesh Number", atoms[i].mesh);
		}
		if(upf.functional_error == 1)
		{
			std::cout << "In Pseudopot_upf::read_pseudo_header : dft_functional from INPUT does not match that in pseudopot file" << std::endl;
			std::cout << "Please make sure this is what you need" << std::endl;
			atoms[i].xc_func = GlobalV::DFT_FUNCTIONAL;
			transform(atoms[i].xc_func.begin(), atoms[i].xc_func.end(), atoms[i].xc_func.begin(), (::toupper));
			if(GlobalV::MY_RANK==0)
			{
				log << "\n In Pseudopot_upf::read_pseudo_header : dft_functional from INPUT does not match that in pseudopot file" << std::endl;
				log << " Please make sure this is what you need" << std::endl;
				log << " XC functional updated to : " << GlobalV::DFT_FUNCTIONAL << std::endl;
				ModuleBase::GlobalFunc::OUT(log,"exchange-correlation functional", atoms[i].xc_func);
			}
		}
			
		//atoms[i].print_pseudo_us(ofs);
	}

//	if(GlobalV::MY_RANK==0)
//	{
//		ofs.close();
//	}
	return;
}
//mock of module_cell/read_cell_pseudopots.cpp
//====================================================
//====================================================
//mock of module_cell/read_pp.cpp
#include "module_cell/read_pp.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <sstream>
#include <cstring> // Peize Lin fix bug about strcpy 2016-08-02


using namespace std;

Pseudopot_upf::Pseudopot_upf()
{
	this->els = new std::string[1];
	this->lchi = nullptr;
	this->oc = nullptr;

	this->r = nullptr;
	this->rab = nullptr;
	this->vloc = nullptr;

	this->kkbeta = nullptr;
	this->lll = nullptr;

	this->rho_at = nullptr;
	this->rho_atc = nullptr;

	this->nn = nullptr;//zhengdy-soc
	this->jchi = nullptr;
	this->jjj = nullptr;

	functional_error = 0;//xiaohui add 2015-03-24
}

Pseudopot_upf::~Pseudopot_upf()
{
	delete [] els; 
	delete [] lchi;
	delete [] oc;

	delete [] r;    //mesh_1
	delete [] rab;  //mesh_2
	delete [] vloc;  //local_1

	delete [] kkbeta; // nl_1
	delete [] lll; // nl_2

	delete [] rho_at;// psrhoatom_1
	delete [] rho_atc;

	delete [] nn;
	delete [] jjj;
	delete [] jchi;
}

int Pseudopot_upf::init_pseudo_reader(const std::string &fn)
{
    ModuleBase::TITLE("Pseudopot_upf","init");
    // First check if this pseudo-potential has spin-orbit information
    std::ifstream ifs(fn.c_str(), ios::in);

	// can't find the file.
	if (!ifs)
    {
        return 1;
    }


    if(GlobalV::global_pseudo_type=="auto") //zws
	{
		set_pseudo_type(fn);
	}

	if(GlobalV::global_pseudo_type=="upf201")
	{
		int info = read_pseudo_upf201(ifs);
		return info;
	}

	return 0;
}


//----------------------------------------------------------
// setting the type of the pseudopotential file
//----------------------------------------------------------
int Pseudopot_upf::set_pseudo_type(const std::string &fn) //zws add
{
    std::ifstream pptype_ifs(fn.c_str(), ios::in);
    std::string dummy;
	std::string strversion;

	if (pptype_ifs.good())
	{
		getline(pptype_ifs,dummy);

		std::stringstream wdsstream(dummy);
		getline(wdsstream,strversion,'"');
		getline(wdsstream,strversion,'"');

		if ( trim(strversion) == "2.0.1" )
		{
			GlobalV::global_pseudo_type = "upf201";
		}
		else
		{
			GlobalV::global_pseudo_type = "upf";
		}
	}
	return 0;
}

std::string& Pseudopot_upf::trim(std::string &in_str)
{
    static const std::string deltri = " \t" ; // delete tab or space
    std::string::size_type position = in_str.find_first_of(deltri, 0);
    if (position == std::string::npos)
	{
        return in_str;
	}
    return trim(in_str.erase(position, 1) );
}

std::string Pseudopot_upf::trimend(std::string &in_str)
{
    const std::string &deltri =" \t" ;
    std::string::size_type position = in_str.find_last_not_of(deltri)+1;
    std::string tmpstr=in_str.erase(position);
    return tmpstr.erase(0,tmpstr.find_first_not_of(deltri));
} //zws


int Pseudopot_upf::average_p(const double& lambda)
{
	int error = 0;
	double lambda_ = lambda;
	if(!GlobalV::LSPINORB) lambda_ = 0.0;
	if(!this->has_so && GlobalV::LSPINORB) 
	{
		error++; 
		std::cout<<"warning_quit! no soc upf used for lspinorb calculation, error!"<<std::endl; 
		return error;
	}
	//ModuleBase::WARNING_QUIT("average_p", "no soc upf used for lspinorb calculation, error!");

	if(!this->has_so || (GlobalV::LSPINORB && abs(lambda_ - 1.0) < 1.0e-8) )
	{
		return error; 
	}

	//if(abs(lambda_)<1.0e-8)
	if(!GlobalV::LSPINORB)
	{
		int new_nbeta = 0; //calculate the new nbeta
		for(int nb=0; nb< this->nbeta; nb++)
		{
			new_nbeta++;
			if(this->lll[nb] != 0 && abs(this->jjj[nb] - this->lll[nb] - 0.5) < 1e-6) //two J = l +- 0.5 average to one
			{
				new_nbeta--;
			}
		}

		this->nbeta = new_nbeta;
		ModuleBase::matrix dion_new;
		dion_new.create(this->nbeta, this->nbeta);

		int old_nbeta=-1;
		for(int nb=0; nb<this->nbeta; nb++)
		{
			old_nbeta++;
			int l = this->lll[old_nbeta];
			int ind=0, ind1=0;
			if(l != 0)
			{
				if(abs(this->jjj[old_nbeta] - this->lll[old_nbeta] + 0.5) < 1e-6)
				{
					if(abs(this->jjj[old_nbeta+1]-this->lll[old_nbeta+1]-0.5)>1e-6) 
					{
						error = 1;
						std::cout<<"warning_quit! error beta function 1 !" <<endl;
						return error;
					}
					ind = old_nbeta +1;
					ind1 = old_nbeta;
				}
				else
				{
					if(abs(this->jjj[old_nbeta+1]-this->lll[old_nbeta+1]+0.5)>1e-6)
					{
						error = 1;
						std::cout<<"warning_quit! error beta function 2 !" <<endl;
						return error;
					}
					ind = old_nbeta;
					ind1 = old_nbeta +1;
				}
				double vion1 = ((l+1.0) * this->dion(ind,ind) + l * this->dion(ind1,ind1)) / (2.0*l+1.0);
				if(std::abs(vion1)<1.0e-8) vion1 = 0.1;
				//average beta (betar)
				for(int ir = 0; ir<this->mesh;ir++)
				{
					this->beta(nb, ir) = 1.0 / (2.0 * l + 1.0) * 
							( (l + 1.0) * sqrt(std::abs(this->dion(ind,ind) / vion1)) *
							this->beta(ind, ir) + 
							l * sqrt(std::abs(this->dion(ind1,ind1) / vion1)) *
							this->beta(ind1, ir) ) ;
				}
				//average the dion matrix
				this->dion(nb, nb) = vion1;
				old_nbeta++;	
			}
			else
			{
				for(int ir = 0; ir<this->mesh;ir++)
					this->beta(nb, ir) = this->beta(old_nbeta, ir);
				this->dion(nb, nb) = this->dion(old_nbeta, old_nbeta);
			}
			this->lll[nb] = this->lll[old_nbeta]; //reset the lll index, ignore jjj index
		}

		//store the old dion and then recreate dion 
		for(int i=0;i<this->nbeta; i++)
		{
			for(int j=0;j<this->nbeta;j++)
			{
				dion_new(i,j) = this->dion(i,j);
			}
		}

		this->dion = dion_new;
		
		int new_nwfc = 0;
		for(int nb=0; nb<this->nwfc; nb++)
		{
			new_nwfc++;
			if(this->lchi[nb] != 0 && abs(this->jchi[nb] - this->lchi[nb] - 0.5)<1e-6)
			{
				new_nwfc--;
			}
		}

		this->nwfc = new_nwfc;
		int old_nwfc=-1;
		for(int nb=0; nb<this->nwfc; nb++)
		{
			old_nwfc++;
			int l = this->lchi[old_nwfc];
			int ind=0, ind1=0;
			if(l!=0)
			{
				if(abs(this->jchi[old_nwfc] - this->lchi[old_nwfc] + 0.5) < 1e-6)
				{
					if(abs(this->jchi[old_nwfc+1]-this->lchi[old_nwfc+1]-0.5)>1e-6) 
					{error++; std::cout<<"warning_quit! error chi function 1 !"<<endl; return error;}
	//					ModuleBase::WARNING_QUIT("average_p", "error chi function 1 !");
					ind = old_nwfc +1;
					ind1 = old_nwfc;
				}
				else
				{
					if(abs(this->jchi[old_nwfc+1]-this->lchi[old_nwfc+1]+0.5)>1e-6)
					{error++; std::cout<<"warning_quit! error chi function 2 !"<<endl; return error;}
	//					ModuleBase::WARNING_QUIT("average_p", "error chi function 2 !");
					ind = old_nwfc;
					ind1 = old_nwfc +1;
				}
				//average chi
				for(int ir = 0; ir<this->mesh;ir++)
				{
					this->chi(nb, ir) = 1.0 / (2.0 * l + 1.0) *
						( (l+1.0)*this->chi(ind,ir) + (l*this->chi(ind1,ir)) );
				}
				old_nwfc++;
			}
			else{
				for(int ir = 0; ir<this->mesh;ir++)
					this->chi(nb, ir) = this->chi(old_nwfc, ir);
			}
			this->lchi[nb] = this->lchi[old_nwfc]; //reset lchi index
		}
		this->has_so = 0;	
		return error;
	}
	else//lambda_ != 0, modulate the soc effect in pseudopotential
	{
		for(int nb=0; nb<this->nbeta; nb++)
		{
			int l = this->lll[nb];
			int ind=0, ind1=0;
			if(l != 0)
			{
				if(abs(this->jjj[nb] - this->lll[nb] + 0.5) < 1e-6)
				{
					if(abs(this->jjj[nb+1]-this->lll[nb+1]-0.5)>1e-6) 
					{
						error = 1;
						std::cout<<"warning_quit! error beta function 1 !" <<endl;
						return error;
					}
					ind = nb +1;
					ind1 = nb;
				}
				else
				{
					if(abs(this->jjj[nb+1]-this->lll[nb+1]+0.5)>1e-6)
					{
						error = 1;
						std::cout<<"warning_quit! error beta function 2 !" <<endl;
						return error;
					}
					ind = nb;
					ind1 = nb +1;
				}
				double vion1 = ((l+1.0) * this->dion(ind,ind) + l * this->dion(ind1,ind1)) / (2.0*l+1.0);
				if(abs(vion1)<1.0e-10) vion1 = 0.1;
				//average beta (betar)
				const double sqrtDplus = sqrt(std::abs(this->dion(ind,ind) / vion1));
				const double sqrtDminus = sqrt(std::abs(this->dion(ind1,ind1) / vion1));
				this->dion(ind, ind) = vion1;
				this->dion(ind1, ind1) = vion1;
				for(int ir = 0; ir<this->mesh;ir++)
				{
					double avera = 1.0 / (2.0 * l + 1.0) * 
							( (l + 1.0) * sqrtDplus *
							this->beta(ind, ir) + 
							l * sqrtDminus *
							this->beta(ind1, ir) ) ;
					double delta = 1.0 / (2.0 * l + 1.0) * 
							( sqrtDplus *
							this->beta(ind, ir) - 
							sqrtDminus *
							this->beta(ind1, ir) ) ;
					this->beta(ind, ir) = (avera + l * delta * lambda_) ;
					this->beta(ind1, ir) = (avera - (l + 1) * delta * lambda_); 
				}
				nb++;
			}
		}

		for(int nb=0; nb<this->nwfc; nb++)
		{
			int l = this->lchi[nb];
			int ind=0, ind1=0;
			if(l!=0)
			{
				if(abs(this->jchi[nb] - this->lchi[nb] + 0.5) < 1e-6)
				{
					if(abs(this->jchi[nb+1]-this->lchi[nb+1]-0.5)>1e-6) 
					{error++; std::cout<<"warning_quit! error chi function 1 !"<<endl; return error;}
					ind = nb +1;
					ind1 = nb;
				}
				else
				{
					if(abs(this->jchi[nb+1]-this->lchi[nb+1]+0.5)>1e-6)
					{error++; std::cout<<"warning_quit! error chi function 2 !"<<endl; return error;}
					ind = nb;
					ind1 = nb +1;
				}
				//average chi
				for(int ir = 0; ir<this->mesh;ir++)
				{
					double avera = 0.5 * 
						( this->chi(ind,ir) + this->chi(ind1,ir) );
					double delta = 0.5 * 
						( this->chi(ind,ir) - this->chi(ind1,ir) );
					this->chi(ind, ir) = avera + delta * lambda_ ; 
					this->chi(ind1, ir) = avera - delta * lambda_ ; 
				}
				nb++;
			}
		}
		return error;
	}
}

// Peize Lin add for bsse 2021.04.07
void Pseudopot_upf::set_empty_element(void)
{
	this->zp = 0;
	for(int ir=0; ir<mesh; ++ir)
	{
		this->vloc[ir] = 0;
	}
	for(int i=0; i<nbeta; ++i)
	{
		for(int j=0; j<nbeta; ++j)
		{
			this->dion(i,j) = 0;
		}
	}
	for(int ir=0; ir<mesh; ++ir)
	{
		this->rho_at[ir] = 0;
	}
	return;
}
//mock of module_cell/read_pp.cpp
//====================================================
//====================================================
//mock of module_cell/read_pp_upf201_mock.cpp
#include "module_cell/read_pp.h"
int Pseudopot_upf::read_pseudo_upf201(std::ifstream &ifs)
{

    std::string word;
    int ONCVPSP;
	//--------------------------------------
	//-              PP_HEADER             - 
	//--------------------------------------
	if(!ModuleBase::GlobalFunc::SCAN_BEGIN(ifs,"<PP_HEADER"))	ModuleBase::WARNING_QUIT("read_pseudo_upf201","Found no PP_HEADER");
	std::string *name=new std::string[50];
	std::string *val=new std::string[50];
	int nparameter;
	this->getnameval(ifs, nparameter, name, val);
	ONCVPSP = 1;
	
	for(int ip = 0 ; ip < nparameter; ++ip)
	{
		if(name[ip]=="generated")
		{
			//add something//
		}
		else if(name[ip]=="author"){}
		else if(name[ip]=="date"){}
		else if(name[ip]=="comment"){}
		else if(name[ip]=="element"){
			psd = val[ip];
		}
		else if(name[ip]=="pseudo_type"){
			pp_type = val[ip];
			if(pp_type!="NC") 
			{
				ModuleBase::WARNING_QUIT("Pseudopot_upf::read_pseudo_header","unknown pseudo type");
			}
		}
		else if(name[ip]=="relativistic"){}
		else if(name[ip]=="is_ultrasoft"){
			if(val[ip]=="T" || val[ip]=="TRUE" || val[ip]=="True" || val[ip]=="true")
			{
				ModuleBase::WARNING_QUIT("Pseudopot_upf::read_pseudo_header","ULTRASOFT PSEUDOPOTENTIAL IS NOT SUPPORTED !!!");
			}
		}
		else if(name[ip]=="is_paw"){
			if(val[ip]=="T" || val[ip]=="TRUE" || val[ip]=="True" || val[ip]=="true")
			{
				ModuleBase::WARNING_QUIT("Pseudopot_upf::read_pseudo_header","PAW PSEUDOPOTENTIAL IS NOT SUPPORTED !!!");
			}
		}
		else if(name[ip]=="is_coulomb"){}
		else if(name[ip]=="has_so"){
			if(val[ip]=="T" || val[ip]=="TRUE" || val[ip]=="True" || val[ip]=="true")
				has_so = true;
			else
				has_so = false;
		}
		else if(name[ip]=="has_wfc"){}
		else if(name[ip]=="has_gipaw"){}
		else if(name[ip]=="paw_as_gipaw"){
			ONCVPSP = 0;
		}
		else if(name[ip]=="core_correction"){
			if(val[ip]=="T" || val[ip]=="TRUE" || val[ip]=="True" || val[ip]=="true")
				nlcc = true;
			else
				nlcc = false;
		}
		else if(name[ip]=="functional"){
			xc_func = val[ip];
		}
		else if(name[ip]=="z_valence"){
			zp = atoi(val[ip].c_str());
		}
		else if(name[ip]=="total_psenergy"){
			etotps = atof(val[ip].c_str());
		}
		else if(name[ip]=="wfc_cutoff"){}
		else if(name[ip]=="rho_cutoff"){}
		else if(name[ip]=="l_max"){
			lmax = atoi(val[ip].c_str());
		}
		else if(name[ip]=="l_max_rho"){}
		else if(name[ip]=="l_local"){}
		else if(name[ip]=="mesh_size"){
			mesh = atoi(val[ip].c_str());
		}
		else if(name[ip]=="number_of_wfc"){
			nwfc = atoi(val[ip].c_str());
		}
		else if(name[ip]=="number_of_proj"){
			nbeta = atoi(val[ip].c_str());
		}
		else
		{
			std::string warningstr = name[ip] + " is not read in. Please add this parameter in read_pp_upf201.cpp if needed.";
			ModuleBase::WARNING("PP_HEADRER reading", warningstr);
		}
	}
			
	//--------------------------------------
	//-              PP_MESH               - 
	//--------------------------------------
	if(ONCVPSP == 0)
	{
		ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_MESH");
		this->getnameval(ifs, nparameter, name, val);
		for(int ip = 0 ; ip < nparameter; ++ip)
		{
			if(name[ip]=="dx"){}
			else if(name[ip]=="mesh"){}
			else if(name[ip]=="xmin"){}
			else if(name[ip]=="rmax"){}
			else if(name[ip]=="zmesh"){}
			else
			{
				std::string warningstr = name[ip] + " is not read in. Please add this parameter in read_pp_upf201.cpp if needed.";
				ModuleBase::WARNING("PP_MESH reading", warningstr);
			}

		}
	}
	else
	{
		ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_MESH>");
	}

	if (ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_R", true, false))
	{
		ModuleBase::GlobalFunc::READ_VALUE(ifs, word); // type size columns
		this->read_pseudo_upf201_r(ifs);
	}
	else if (ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_R>"))
	{
		this->read_pseudo_upf201_r(ifs);
	}
	ModuleBase::GlobalFunc::SCAN_END(ifs, "</PP_R>");

    if (ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_RAB", true, false))
	{
		ModuleBase::GlobalFunc::READ_VALUE(ifs, word); // type size columns
		this->read_pseudo_upf201_rab(ifs);
	}
	else if (ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_RAB>"))
	{
		this->read_pseudo_upf201_rab(ifs);
	}
	ModuleBase::GlobalFunc::SCAN_END(ifs, "</PP_RAB>");
	ModuleBase::GlobalFunc::SCAN_END(ifs, "</PP_MESH>");

	//--------------------------------------
	//-              PP_NLCC               - 
	//--------------------------------------
	if (this->nlcc)
	{
		ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_NLCC");
		ModuleBase::GlobalFunc::READ_VALUE(ifs, word);    // type size columns
		delete[] rho_atc;
		this->rho_atc = new double[mesh];
		ModuleBase::GlobalFunc::ZEROS(rho_atc, mesh);
		for (int ir = 0;ir < mesh;ir++)
		{
			ifs >> this->rho_atc[ir];
		}
		ModuleBase::GlobalFunc::SCAN_END(ifs, "</PP_NLCC>");
	}

	//--------------------------------------
	//-              PP_LOCAL              - 
	//--------------------------------------
	ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_LOCAL");
	ModuleBase::GlobalFunc::READ_VALUE(ifs, word);    // type size columns
	delete[] vloc;
	this->vloc = new double[mesh];
	ModuleBase::GlobalFunc::ZEROS(vloc, mesh);
	for (int ir = 0;ir < mesh;ir++)
	{
		ifs >> this->vloc[ir];
	}
	ModuleBase::GlobalFunc::SCAN_END(ifs, "</PP_LOCAL>");

	//--------------------------------------
	//-            PP_NONLOCAL             - 
	//--------------------------------------
	ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_NONLOCAL>");
	delete[] kkbeta;
	delete[] lll;
	this->kkbeta = new int[nbeta];
	this->lll = new int[nbeta];
	this->beta.create(nbeta , mesh);
	this->dion.create(nbeta , nbeta);
	for(int ib=0;ib<nbeta;ib++)
	{
		ifs >> word; //number of beta
		this->getnameval(ifs, nparameter, name, val);
		for(int ip = 0 ; ip < nparameter; ++ip)
		{
			if(name[ip]=="type"){}
			else if(name[ip]=="size"){}
			else if(name[ip]=="columns"){}
			else if(name[ip]=="index"){}
			else if(name[ip]=="label"){}
			else if(name[ip]=="angular_momentum"){
				lll[ib] = atoi(val[ip].c_str());
			}
			else if(name[ip]=="cutoff_radius_index"){
				kkbeta[ib] = atoi (val[ip].c_str());
			}
			else if(name[ip]=="cutoff_radius"){}
			else if(name[ip]=="ultrasoft_cutoff_radius"){}
			else
			{
				std::string warningstr = name[ip] + " is not read in. Please add this parameter in read_pp_upf201.cpp if needed.";
				ModuleBase::WARNING("PP_BETA reading", warningstr);
			}
		}
		for (int ir=0;ir<mesh;ir++)
		{
			ifs >> this->beta(ib, ir);
		}
		ifs >> word; //number of beta
	}

	if (ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_DIJ", true, false))
	{
		ModuleBase::GlobalFunc::READ_VALUE(ifs, word);  // type size columns
		this->read_pseudo_upf201_dij(ifs);
	}
	else if ( ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_DIJ>"))
	{
		this->read_pseudo_upf201_dij(ifs);
	}
	ModuleBase::GlobalFunc::SCAN_END(ifs, "</PP_DIJ>");
	ModuleBase::GlobalFunc::SCAN_END(ifs, "</PP_NONLOCAL>");

	//--------------------------------------
	//-            PP_PSWFC                - 
	//--------------------------------------
	ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_PSWFC>");
	delete[] els;
	delete[] lchi;
	delete[] oc;
	this->els = new std::string[nwfc];
	this->lchi = new int[nwfc];
	this->oc = new double[nwfc];
	ModuleBase::GlobalFunc::ZEROS(lchi, nwfc); // angular momentum of each orbital
	ModuleBase::GlobalFunc::ZEROS(oc, nwfc);//occupation of each orbital
	this->chi.create(this->nwfc, this->mesh);
	for(int iw=0;iw<nwfc;iw++)
	{
		ifs >> word; //number of chi
		this->getnameval(ifs, nparameter, name, val);
		for(int ip = 0 ; ip < nparameter; ++ip)
		{
			if(name[ip]=="type"){}
			else if(name[ip]=="size"){}
			else if(name[ip]=="columns"){}
			else if(name[ip]=="index"){}
			else if(name[ip]=="label"){
				els[iw] = val[ip];
			}
			else if(name[ip]=="l"){
				lchi[iw] = atoi(val[ip].c_str());
			}
			else if(name[ip]=="occupation"){
				oc[iw] = atof(val[ip].c_str());
			}
			else if(name[ip]=="n"){}
			else if(name[ip]=="pseudo_energy"){}
			else if(name[ip]=="cutoff_radius"){}
			else if(name[ip]=="ultrasoft_cutoff_radius"){}
			else
			{
				std::string warningstr = name[ip] + " is not read in. Please add this parameter in read_pp_upf201.cpp if needed.";
				ModuleBase::WARNING("PP_CHI reading", warningstr);
			}
		}
		for (int ir=0;ir<mesh;ir++)
		{
			ifs >> this->chi(iw, ir);
		}
		ifs >> word; //number of chi
	}
	ModuleBase::GlobalFunc::SCAN_END(ifs, "</PP_PSWFC>");

	//--------------------------------------
	//-          PP_RHOATOM                - 
	//--------------------------------------
	if (ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_RHOATOM", true, false))
	{
		ModuleBase::GlobalFunc::READ_VALUE(ifs, word); // type size columns
		this->read_pseudo_upf201_rhoatom(ifs);
	}
	else if (ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_RHOATOM>"))
	{
		this->read_pseudo_upf201_rhoatom(ifs);
	}
	ModuleBase::GlobalFunc::SCAN_END(ifs, "</PP_RHOATOM>");

	//--------------------------------------
	//-          PP_SPIN_ORB               - 
	//--------------------------------------
	ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_SPIN_ORB>", true, false);
	//added by zhengdy-soc
	delete[] this->jchi;
	delete[] this->jjj;
	delete[] this->nn;
	this->jchi = new double [nwfc];
	this->jjj = new double [nbeta];
	this->nn = new int [nwfc];
	ModuleBase::GlobalFunc::ZEROS(jchi,nwfc);
	ModuleBase::GlobalFunc::ZEROS(jjj,nbeta);
	ModuleBase::GlobalFunc::ZEROS(nn,nwfc);

	for(int round=0;round<2;round++)
	{
		ifs>>word;
		if(word=="<PP_RELBETA.1")
		{
			for(int nb = 0;nb<nbeta;nb++)
			{
				if(nb!=0) ifs>>word; //RELBETA
				this->getnameval(ifs, nparameter, name, val);
				for(int ip = 0 ; ip < nparameter; ++ip)
				{
					if(name[ip]=="index"){}
					else if(name[ip]=="lll"){
						lll[nb] = atoi(val[ip].c_str());
					}
					else if(name[ip]=="jjj"){
						jjj[nb] = atof(val[ip].c_str());
					}
					else
					{
						std::string warningstr = name[ip] + " is not read in. Please add this parameter in read_pp_upf201.cpp if needed.";
						ModuleBase::WARNING("PP_RELBETA reading", warningstr);
					}
				}
			}
		}
		else if(word=="<PP_RELWFC.1")
		{
			for(int nw = 0;nw<nwfc;nw++)
			{
				if(nw!=0) ifs>>word;     //RELWFC
				this->getnameval(ifs, nparameter, name, val);
				for(int ip = 0 ; ip < nparameter; ++ip)
				{
					if(name[ip]=="index"){}
					else if(name[ip]=="els"){
						els[nw] = val[ip];
					}
					else if(name[ip]=="nn"){
						nn[nw] = atoi(val[ip].c_str());
					}
					else if(name[ip]=="lchi"){
						lchi[nw]=atoi(val[ip].c_str());
					}
					else if(name[ip]=="jchi"){
						jchi[nw]=atof(val[ip].c_str());
					}
					else if(name[ip]=="oc"){
						oc[nw] = atof(val[ip].c_str());
					}
					else
					{
						std::string warningstr = name[ip] + " is not read in. Please add this parameter in read_pp_upf201.cpp if needed.";
						ModuleBase::WARNING("PP_RELWFC reading", warningstr);
					}
				}
			}
		}
		else if(round==0)
		{
			this->has_so = 0;
			break;
		}
	}
	ModuleBase::GlobalFunc::SCAN_END(ifs, "</PP_SPIN_ORB>", false);
	if (mesh%2 == 0)
	{
		mesh -= 1;
	}
	
	ModuleBase::GlobalFunc::SCAN_END(ifs, "</UPF>", false);
	delete []name;
	delete []val;
	
	if(GlobalV::DFT_FUNCTIONAL!="default")
	{
		string xc_func1 = GlobalV::DFT_FUNCTIONAL;
		transform(xc_func1.begin(), xc_func1.end(), xc_func1.begin(), (::toupper));
		if(xc_func1 != xc_func)
		{
			functional_error = 1;

			std::cout << " dft_functional readin is: " << GlobalV::DFT_FUNCTIONAL << std::endl;
			std::cout << " dft_functional in pseudopot file is: " << xc_func << std::endl;
			GlobalV::ofs_warning << " dft_functional readin is: " << GlobalV::DFT_FUNCTIONAL << std::endl;
			GlobalV::ofs_warning << " dft_functional in pseudopot file is: " << xc_func << std::endl;
		}
	}
	return 0;

}
void Pseudopot_upf:: getnameval(std::ifstream& ifs,int &n, std::string * name, std::string *val)
{
	std::string txt,word;
	//get long txt
	ifs>>txt;
	while(ifs>>word)
	{
		size_t wl= word.length()-1;
		txt=txt+" "+word;
		if(word.substr(wl,1)==">")
		{
			break;
		}
	}

	//count number of parameters according to "="
	size_t pos=0;
	n=0;
	while(1)
	{
		pos = txt.find("=",pos);
		if(pos == std::string::npos) break;
		pos++;
		n++;
	}

	//get name & value
	pos=0;
	size_t pos2, ll;
	for(int i = 0; i < n; ++i)
	{
		pos2 = txt.find("=",pos);
		for(; pos2 > pos ; --pos2)//There may be a space before "=";
		{
			if(txt.substr(pos2-1,1) != " ") break;
		}
		ll=pos2-pos;
		name[i] = txt.substr(pos,ll);
		std::string mark;
		bool findmark=false;
		for(int j = 0; j < 100; ++j)//The mark can be ' or " or .
		{
			mark = txt.substr(pos2,1);
			pos2++;
			if(mark=="\""||mark=="\'"||mark==".")
			{
				findmark = true;
				break;
			}
		}
		if(!findmark) ModuleBase::WARNING_QUIT("read_upf201",
		"The values are not in \' or \". Please improve the program in read_pp_upf201.cpp");
		pos = pos2;
		pos2 = txt.find(mark,pos);
		ll=pos2-pos;
		std::string tmpval = txt.substr(pos,ll);
		tmpval = trim(tmpval);
		val[i]=tmpval;
		pos=pos2+1;
		for(int j = 0; j < 100; ++j)
		{
			if(txt.substr(pos,1)==" " || txt.substr(pos,1)==",")
				pos++;
			else
				break;
		}
	}
	return;
}

void Pseudopot_upf::read_pseudo_upf201_r(std::ifstream &ifs)
{
	delete[] r;
	delete[] rab;
	assert(mesh>0);
	this->r = new double[mesh];
	this->rab = new double[mesh];
	ModuleBase::GlobalFunc::ZEROS(r,mesh);
	ModuleBase::GlobalFunc::ZEROS(rab,mesh);
	for (int ir = 0;ir < mesh;ir++)
	{
		ifs >> this->r[ir];
	}
}
void Pseudopot_upf::read_pseudo_upf201_rab(std::ifstream &ifs)
{
	for (int ir = 0;ir < mesh;ir++)
	{
		ifs >> this->rab[ir];
	}
}
void Pseudopot_upf::read_pseudo_upf201_dij(std::ifstream &ifs)
{
	this->nd = nbeta * nbeta;
	for(int i=0;i<nbeta;i++)
	{
		for(int j=0;j<nbeta;j++)
		{
			ifs >> dion(i,j);
			if ( i != j  && dion(i,j) != 0.0 )
			{
				ModuleBase::WARNING_QUIT("read_pseudo_upf201","Error: for i != j, Dij of Pseudopotential must be 0.0");
			}
		}
	}
}
void Pseudopot_upf::read_pseudo_upf201_rhoatom(std::ifstream &ifs)
{
	delete[] rho_at;
	this->rho_at = new double[mesh];
	ModuleBase::GlobalFunc::ZEROS(rho_at, mesh);
	for (int ir = 0;ir < mesh;ir++)
	{
		ifs >> this->rho_at[ir];
	}
}

#include "src_io/berryphase.h"
bool berryphase::berry_phase_flag;
