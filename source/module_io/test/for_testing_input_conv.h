#define private public

#include "module_cell/module_symmetry/symmetry.h"
#include "module_cell/unitcell.h"
#include "module_elecstate/elecstate_lcao.h"
#include "module_elecstate/module_charge/charge_mixing.h"
#include "module_elecstate/occupy.h"
#include "module_elecstate/potentials/H_TDDFT_pw.h"
#include "module_elecstate/potentials/efield.h"
#include "module_elecstate/potentials/gatefield.h"
#include "module_hamilt_lcao/hamilt_lcaodft/FORCE_STRESS.h"
#include "module_hamilt_lcao/hamilt_lcaodft/local_orbital_charge.h"
#include "module_hamilt_lcao/module_dftu/dftu.h"
#include "module_hamilt_lcao/module_tddft/evolve_elec.h"
#include "module_hamilt_pw/hamilt_pwdft/VNL_in_pw.h"
#include "module_hamilt_pw/hamilt_pwdft/structure_factor.h"
#include "module_hamilt_pw/hamilt_pwdft/wavefunc.h"
#include "module_hsolver/hsolver_lcao.h"
#include "module_io/berryphase.h"
#include "module_io/restart.h"
#include "module_md/md_func.h"
#include "module_relax/relax_old/bfgs_basic.h"
#include "module_relax/relax_old/ions_move_basic.h"
#include "module_relax/relax_old/ions_move_cg.h"
#include "module_relax/relax_old/lattice_change_basic.h"

bool berryphase::berry_phase_flag = false;
int elecstate::ElecStateLCAO::out_wfc_lcao = 0;
bool elecstate::ElecStateLCAO::need_psi_grid = 1;
int hsolver::HSolverLCAO::out_mat_hs = 0;
int hsolver::HSolverLCAO::out_mat_hsR = 0;
int hsolver::HSolverLCAO::out_mat_t = 0;
int hsolver::HSolverLCAO::out_mat_dh = 0;
int Local_Orbital_Charge::out_dm = 0;
int Local_Orbital_Charge::out_dm1 = 0;
double module_tddft::Evolve_elec::td_force_dt;
bool module_tddft::Evolve_elec::td_vext;
std::vector<int> module_tddft::Evolve_elec::td_vext_dire_case;
bool module_tddft::Evolve_elec::out_dipole;
bool module_tddft::Evolve_elec::out_efield;
double module_tddft::Evolve_elec::td_print_eij;
int module_tddft::Evolve_elec::td_edm;
double elecstate::Gatefield::zgate = 0.5;
bool elecstate::Gatefield::relax = false;
bool elecstate::Gatefield::block = false;
double elecstate::Gatefield::block_down = 0.45;
double elecstate::Gatefield::block_up = 0.55;
double elecstate::Gatefield::block_height = 0.1;
int elecstate::Efield::efield_dir;
double elecstate::Efield::efield_pos_max;
double elecstate::Efield::efield_pos_dec;
double elecstate::Efield::efield_amp;

// parameters of electric field for tddft

int elecstate::H_TDDFT_pw::stype; // 0 : length gauge  1: velocity gauge

std::vector<int> elecstate::H_TDDFT_pw::ttype;
//  0  Gauss type function.
//  1  trapezoid type function.
//  2  Trigonometric functions, sin^2.
//  3  heaviside function.
//  4  HHG function.

int elecstate::H_TDDFT_pw::tstart;
int elecstate::H_TDDFT_pw::tend;
double elecstate::H_TDDFT_pw::dt;

// space domain parameters

// length gauge
double elecstate::H_TDDFT_pw::lcut1;
double elecstate::H_TDDFT_pw::lcut2;

// time domain parameters

// Gauss
int elecstate::H_TDDFT_pw::gauss_count;
std::vector<double> elecstate::H_TDDFT_pw::gauss_omega; // time(a.u.)^-1
std::vector<double> elecstate::H_TDDFT_pw::gauss_phase;
std::vector<double> elecstate::H_TDDFT_pw::gauss_sigma; // time(a.u.)
std::vector<double> elecstate::H_TDDFT_pw::gauss_t0;
std::vector<double> elecstate::H_TDDFT_pw::gauss_amp; // Ry/bohr

// trapezoid
int elecstate::H_TDDFT_pw::trape_count;
std::vector<double> elecstate::H_TDDFT_pw::trape_omega; // time(a.u.)^-1
std::vector<double> elecstate::H_TDDFT_pw::trape_phase;
std::vector<double> elecstate::H_TDDFT_pw::trape_t1;
std::vector<double> elecstate::H_TDDFT_pw::trape_t2;
std::vector<double> elecstate::H_TDDFT_pw::trape_t3;
std::vector<double> elecstate::H_TDDFT_pw::trape_amp; // Ry/bohr

// Trigonometric
int elecstate::H_TDDFT_pw::trigo_count;
std::vector<double> elecstate::H_TDDFT_pw::trigo_omega1; // time(a.u.)^-1
std::vector<double> elecstate::H_TDDFT_pw::trigo_omega2; // time(a.u.)^-1
std::vector<double> elecstate::H_TDDFT_pw::trigo_phase1;
std::vector<double> elecstate::H_TDDFT_pw::trigo_phase2;
std::vector<double> elecstate::H_TDDFT_pw::trigo_amp; // Ry/bohr

// Heaviside
int elecstate::H_TDDFT_pw::heavi_count;
std::vector<double> elecstate::H_TDDFT_pw::heavi_t0;
std::vector<double> elecstate::H_TDDFT_pw::heavi_amp; // Ry/bohr

double Force_Stress_LCAO::force_invalid_threshold_ev = 0.0;
double BFGS_Basic::relax_bfgs_w1 = -1.0;
double BFGS_Basic::relax_bfgs_w2 = -1.0;
double Ions_Move_Basic::relax_bfgs_rmax = -1.0;
double Ions_Move_Basic::relax_bfgs_rmin = -1.0;
double Ions_Move_Basic::relax_bfgs_init = -1.0;
int Ions_Move_Basic::out_stru = 0;
double Ions_Move_CG::RELAX_CG_THR = -1.0;
std::string Lattice_Change_Basic::fixed_axes = "None";
int ModuleSymmetry::Symmetry::symm_flag = 0;
bool ModuleSymmetry::Symmetry::symm_autoclose = false;

Charge_Mixing::Charge_Mixing()
{
}
Charge_Mixing::~Charge_Mixing()
{
}
pseudopot_cell_vnl::pseudopot_cell_vnl()
{
}
pseudopot_cell_vnl::~pseudopot_cell_vnl()
{
}
pseudopot_cell_vl::pseudopot_cell_vl()
{
}
pseudopot_cell_vl::~pseudopot_cell_vl()
{
}
ORB_gaunt_table::ORB_gaunt_table()
{
}
ORB_gaunt_table::~ORB_gaunt_table()
{
}
ModuleDFTU::DFTU::DFTU()
{
}
ModuleDFTU::DFTU::~DFTU()
{
}
Structure_Factor::Structure_Factor()
{
}
Structure_Factor::~Structure_Factor()
{
}
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
UnitCell::UnitCell()
{
    if (GlobalV::test_unitcell)
        ModuleBase::TITLE("unitcell", "Constructor");
    Coordinate = "Direct";
    latName = "none";
    lat0 = 0.0;
    lat0_angstrom = 0.0;

    bool init_vel;

    ntype = 0;
    nat = 0;
    namax = 0;
    nwmax = 0;

    iat2it = nullptr;
    iat2ia = nullptr;
    iwt2iat = nullptr;
    iwt2iw = nullptr;

    itia2iat.create(1, 1);
    lc = new int[3];

    latvec = ModuleBase::Matrix3();
    latvec_supercell = ModuleBase::Matrix3();
    G = ModuleBase::Matrix3();
    GT = ModuleBase::Matrix3();
    GGT = ModuleBase::Matrix3();
    invGGT = ModuleBase::Matrix3();

    tpiba = 0.0;
    tpiba2 = 0.0;
    omega = 0.0;

    atom_label = new std::string[1];
    atom_mass = nullptr;
    pseudo_fn = new std::string[1];
    pseudo_type = new std::string[1];
    orbital_fn = new std::string[1];

    set_atom_flag = false;
}
UnitCell::~UnitCell()
{
}
#ifdef __LCAO
InfoNonlocal::InfoNonlocal()
{
}
InfoNonlocal::~InfoNonlocal()
{
}
LCAO_Orbitals::LCAO_Orbitals()
{
}
LCAO_Orbitals::~LCAO_Orbitals()
{
}
#endif
Magnetism::Magnetism()
{
}
Magnetism::~Magnetism()
{
}

void Charge_Mixing::set_mixing(const std::string& mixing_mode_in,
                               const double& mixing_beta_in,
                               const int& mixing_ndim_in,
                               const double& mixing_gg0_in,
                               const bool& mixing_tau_in)
{
    return;
}
// void Charge_Mixing::need_auto_set(){}
void Charge_Mixing::need_auto_set()
{
    this->autoset = true;
}
void Occupy::decision(const std::string& name, const std::string& smearing_method, const double& smearing_sigma)
{
    return;
}
// void UnitCell::setup(const std::string&,const int&,const int&,const bool&,const std::string&){return;}
void UnitCell::setup(const std::string& latname_in,
                     const int& ntype_in,
                     const int& lmaxmax_in,
                     const bool& init_vel_in,
                     const std::string& fixed_axes_in)
{
    this->latName = latname_in;
    this->ntype = ntype_in;
    this->lmaxmax = lmaxmax_in;
    this->init_vel = init_vel_in;
    // pengfei Li add 2018-11-11
    if (fixed_axes_in == "None")
    {
        this->lc[0] = 1;
        this->lc[1] = 1;
        this->lc[2] = 1;
    }
    else if (fixed_axes_in == "volume")
    {
        this->lc[0] = 1;
        this->lc[1] = 1;
        this->lc[2] = 1;
        if (!GlobalV::relax_new)
        {
            ModuleBase::WARNING_QUIT(
                "Input",
                "there are bugs in the old implementation; set relax_new to be 1 for fixed_volume relaxation");
        }
    }
    else if (fixed_axes_in == "shape")
    {
        if (!GlobalV::relax_new)
        {
            ModuleBase::WARNING_QUIT("Input", "set relax_new to be 1 for fixed_shape relaxation");
        }
        this->lc[0] = 1;
        this->lc[1] = 1;
        this->lc[2] = 1;
    }
    else if (fixed_axes_in == "a")
    {
        this->lc[0] = 0;
        this->lc[1] = 1;
        this->lc[2] = 1;
    }
    else if (fixed_axes_in == "b")
    {
        this->lc[0] = 1;
        this->lc[1] = 0;
        this->lc[2] = 1;
    }
    else if (fixed_axes_in == "c")
    {
        this->lc[0] = 1;
        this->lc[1] = 1;
        this->lc[2] = 0;
    }
    else if (fixed_axes_in == "ab")
    {
        this->lc[0] = 0;
        this->lc[1] = 0;
        this->lc[2] = 1;
    }
    else if (fixed_axes_in == "ac")
    {
        this->lc[0] = 0;
        this->lc[1] = 1;
        this->lc[2] = 0;
    }
    else if (fixed_axes_in == "bc")
    {
        this->lc[0] = 1;
        this->lc[1] = 0;
        this->lc[2] = 0;
    }
    else if (fixed_axes_in == "abc")
    {
        this->lc[0] = 0;
        this->lc[1] = 0;
        this->lc[2] = 0;
    }
    else
    {
        ModuleBase::WARNING_QUIT("Input", "fixed_axes should be None,volume,shape,a,b,c,ab,ac,bc or abc!");
    }
    return;
}
void Structure_Factor::set(const int&)
{
    return;
}

namespace MD_func
{
void current_md_info(const int& my_rank, const std::string& file_dir, int& md_step, double& temperature)
{
    return;
}
} // namespace MD_func

namespace GlobalC
{
UnitCell ucell;
wavefunc wf;
ModuleSymmetry::Symmetry symm;
Structure_Factor sf;
ModuleDFTU::DFTU dftu;
Restart restart;
pseudopot_cell_vnl ppcell;
Charge_Mixing CHR_MIX;
} // namespace GlobalC

#undef private
