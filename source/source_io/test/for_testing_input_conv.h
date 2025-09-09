#ifndef INPUT_CONV_TEST_H
#define INPUT_CONV_TEST_H
#define private public
#include "source_io/module_parameter/parameter.h"
#include "source_cell/module_symmetry/symmetry.h"
#include "source_cell/unitcell.h"
#include "source_estate/elecstate_lcao.h"
#include "source_estate/module_charge/charge_mixing.h"
#include "source_estate/occupy.h"
#include "source_estate/module_pot/H_TDDFT_pw.h"
#include "source_estate/module_pot/efield.h"
#include "source_estate/module_pot/gatefield.h"
#include "source_lcao/FORCE_STRESS.h"
#include "source_lcao/module_dftu/dftu.h"
#include "source_lcao/module_rt/evolve_elec.h"
#include "source_lcao/module_rt/td_velocity.h"
#include "source_pw/module_pwdft/VNL_in_pw.h"
#include "source_pw/module_pwdft/structure_factor.h"
#include "source_hsolver/hsolver_lcao.h"
#include "source_io/berryphase.h"
#include "source_io/restart.h"
#include "source_md/md_func.h"
#include "source_relax/bfgs_basic.h"
#include "source_relax/ions_move_basic.h"
#include "source_relax/ions_move_cg.h"
#include "source_relax/lattice_change_basic.h"
#ifdef __PEXSI
#include "source_hsolver/module_pexsi/pexsi_solver.h"
#endif
#undef private
bool berryphase::berry_phase_flag = false;

bool TD_Velocity::out_current;
bool TD_Velocity::out_current_k;
bool TD_Velocity::out_vecpot;
bool TD_Velocity::init_vecpot_file;
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

// Parameters of electric field for RT-TDDFT

int elecstate::H_TDDFT_pw::stype;

std::vector<int> elecstate::H_TDDFT_pw::ttype;

int elecstate::H_TDDFT_pw::tstart;
int elecstate::H_TDDFT_pw::tend;
double elecstate::H_TDDFT_pw::dt;
double elecstate::H_TDDFT_pw::dt_int;

double elecstate::H_TDDFT_pw::lcut1;
double elecstate::H_TDDFT_pw::lcut2;

bool TD_Velocity::tddft_velocity;
bool TD_Velocity::out_mat_R;

// Gaussian
int elecstate::H_TDDFT_pw::gauss_count;
std::vector<double> elecstate::H_TDDFT_pw::gauss_omega; // time(a.u.)^-1
std::vector<double> elecstate::H_TDDFT_pw::gauss_phase;
std::vector<double> elecstate::H_TDDFT_pw::gauss_sigma; // time(a.u.)
std::vector<double> elecstate::H_TDDFT_pw::gauss_t0;
std::vector<double> elecstate::H_TDDFT_pw::gauss_amp; // Ry/bohr
std::vector<int> elecstate::H_TDDFT_pw::gauss_ncut;

// Trapezoid
int elecstate::H_TDDFT_pw::trape_count;
std::vector<double> elecstate::H_TDDFT_pw::trape_omega; // time(a.u.)^-1
std::vector<double> elecstate::H_TDDFT_pw::trape_phase;
std::vector<double> elecstate::H_TDDFT_pw::trape_t1;
std::vector<double> elecstate::H_TDDFT_pw::trape_t2;
std::vector<double> elecstate::H_TDDFT_pw::trape_t3;
std::vector<double> elecstate::H_TDDFT_pw::trape_amp; // Ry/bohr
std::vector<int> elecstate::H_TDDFT_pw::trape_ncut;

// Trigonometric
int elecstate::H_TDDFT_pw::trigo_count;
std::vector<double> elecstate::H_TDDFT_pw::trigo_omega1; // time(a.u.)^-1
std::vector<double> elecstate::H_TDDFT_pw::trigo_omega2; // time(a.u.)^-1
std::vector<double> elecstate::H_TDDFT_pw::trigo_phase1;
std::vector<double> elecstate::H_TDDFT_pw::trigo_phase2;
std::vector<double> elecstate::H_TDDFT_pw::trigo_amp; // Ry/bohr
std::vector<int> elecstate::H_TDDFT_pw::trigo_ncut;

// Heaviside
int elecstate::H_TDDFT_pw::heavi_count;
std::vector<double> elecstate::H_TDDFT_pw::heavi_t0;
std::vector<double> elecstate::H_TDDFT_pw::heavi_amp; // Ry/bohr

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
Soc::~Soc()
{
}
Fcoef::~Fcoef()
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
UnitCell::UnitCell()
{
    itia2iat.create(1, 1);
}
UnitCell::~UnitCell() {}
#ifdef __LCAO
InfoNonlocal::InfoNonlocal() {}
InfoNonlocal::~InfoNonlocal() {}
LCAO_Orbitals::LCAO_Orbitals() {}
LCAO_Orbitals::~LCAO_Orbitals() {}
#endif
Magnetism::Magnetism() {}
Magnetism::~Magnetism() {}
void Occupy::decision(const std::string& name,
                      const std::string& smearing_method,
                      const double& smearing_sigma) {
    return;
}
// void UnitCell::setup(const std::string&,const int&,const int&,const
// bool&,const std::string&){return;}
void UnitCell::setup(const std::string& latname_in,
                     const int& ntype_in,
                     const int& lmaxmax_in,
                     const bool& init_vel_in,
                     const std::string& fixed_axes_in) {
    this->latName = latname_in;
    this->ntype = ntype_in;
    this->lmaxmax = lmaxmax_in;
    this->init_vel = init_vel_in;
    // pengfei Li add 2018-11-11
    if (fixed_axes_in == "None") {
        this->lc[0] = 1;
        this->lc[1] = 1;
        this->lc[2] = 1;
    } else if (fixed_axes_in == "volume") {
        this->lc[0] = 1;
        this->lc[1] = 1;
        this->lc[2] = 1;
        if (!PARAM.input.relax_new) {
            ModuleBase::WARNING_QUIT(
                "Input",
                "there are bugs in the old implementation; set relax_new to be "
                "1 for fixed_volume relaxation");
        }
    } else if (fixed_axes_in == "shape") {
        if (!PARAM.input.relax_new) {
            ModuleBase::WARNING_QUIT(
                "Input",
                "set relax_new to be 1 for fixed_shape relaxation");
        }
        this->lc[0] = 1;
        this->lc[1] = 1;
        this->lc[2] = 1;
    } else if (fixed_axes_in == "a") {
        this->lc[0] = 0;
        this->lc[1] = 1;
        this->lc[2] = 1;
    } else if (fixed_axes_in == "b") {
        this->lc[0] = 1;
        this->lc[1] = 0;
        this->lc[2] = 1;
    } else if (fixed_axes_in == "c") {
        this->lc[0] = 1;
        this->lc[1] = 1;
        this->lc[2] = 0;
    } else if (fixed_axes_in == "ab") {
        this->lc[0] = 0;
        this->lc[1] = 0;
        this->lc[2] = 1;
    } else if (fixed_axes_in == "ac") {
        this->lc[0] = 0;
        this->lc[1] = 1;
        this->lc[2] = 0;
    } else if (fixed_axes_in == "bc") {
        this->lc[0] = 1;
        this->lc[1] = 0;
        this->lc[2] = 0;
    } else if (fixed_axes_in == "abc") {
        this->lc[0] = 0;
        this->lc[1] = 0;
        this->lc[2] = 0;
    } else {
        ModuleBase::WARNING_QUIT(
            "Input",
            "fixed_axes should be None,volume,shape,a,b,c,ab,ac,bc or abc!");
    }
    return;
}
// void Structure_Factor::set(const int&)
// {
//     return;
// }

namespace MD_func {
void current_md_info(const int& my_rank,
                     const std::string& file_dir,
                     int& md_step,
                     double& temperature) {
    return;
}
} // namespace MD_func

namespace GlobalC {
ModuleDFTU::DFTU dftu;
Restart restart;
} // namespace GlobalC

#ifdef __PEXSI
namespace pexsi {
int PEXSI_Solver::pexsi_npole = 0;
bool PEXSI_Solver::pexsi_inertia = 0;
int PEXSI_Solver::pexsi_nmax = 0;
// int PEXSI_Solver::pexsi_symbolic = 0;
bool PEXSI_Solver::pexsi_comm = 0;
bool PEXSI_Solver::pexsi_storage = 0;
int PEXSI_Solver::pexsi_ordering = 0;
int PEXSI_Solver::pexsi_row_ordering = 0;
int PEXSI_Solver::pexsi_nproc = 0;
bool PEXSI_Solver::pexsi_symm = 0;
bool PEXSI_Solver::pexsi_trans = 0;
int PEXSI_Solver::pexsi_method = 0;
int PEXSI_Solver::pexsi_nproc_pole = 0;
// double PEXSI_Solver::pexsi_spin = 2;
double PEXSI_Solver::pexsi_temp = 0.0;
double PEXSI_Solver::pexsi_gap = 0.0;
double PEXSI_Solver::pexsi_delta_e = 0.0;
double PEXSI_Solver::pexsi_mu_lower = 0.0;
double PEXSI_Solver::pexsi_mu_upper = 0.0;
double PEXSI_Solver::pexsi_mu = 0.0;
double PEXSI_Solver::pexsi_mu_thr = 0.0;
double PEXSI_Solver::pexsi_mu_expand = 0.0;
double PEXSI_Solver::pexsi_mu_guard = 0.0;
double PEXSI_Solver::pexsi_elec_thr = 0.0;
double PEXSI_Solver::pexsi_zero_thr = 0.0;
} // namespace pexsi
#endif

#endif
