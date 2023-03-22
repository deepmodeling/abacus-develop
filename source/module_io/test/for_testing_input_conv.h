#define private public

#include "module_cell/unitcell.h"
#include "module_hamilt_pw/hamilt_pwdft/wavefunc.h"
#include "module_hamilt_lcao/hamilt_lcaodft/FORCE_STRESS.h"
#include "module_relax/relax_old/bfgs_basic.h"
#include "module_relax/relax_old/ions_move_basic.h"
#include "module_relax/relax_old/lattice_change_basic.h"
#include "module_relax/relax_old/ions_move_cg.h"
#include "module_cell/module_symmetry/symmetry.h"
#include "module_hamilt_pw/hamilt_pwdft/structure_factor.h"
#include "module_elecstate/energy.h"
#include "module_hamilt_lcao/module_dftu/dftu.h"
#include "module_elecstate/potentials/efield.h"
#include "module_elecstate/potentials/gatefield.h"
#include "module_elecstate/potentials/gatefield.h"
#include "module_hamilt_lcao/module_tddft/ELEC_evolve.h"
#include "module_io/restart.h"
#include "module_hamilt_pw/hamilt_pwdft/VNL_in_pw.h"
#include "module_elecstate/occupy.h"
#include "module_elecstate/module_charge/charge_mixing.h"
#include "module_hamilt_lcao/hamilt_lcaodft/local_orbital_charge.h"
#include "module_hsolver/hsolver_lcao.h"
#include "module_elecstate/elecstate_lcao.h"
#include "module_io/berryphase.h"

bool berryphase::berry_phase_flag=false;
int elecstate::ElecStateLCAO::out_wfc_lcao = 0;
bool elecstate::ElecStateLCAO::need_psi_grid = 1;
int hsolver::HSolverLCAO::out_mat_hs = 0;
int hsolver::HSolverLCAO::out_mat_hsR = 0;
int hsolver::HSolverLCAO::out_hsR_interval = 1;
int hsolver::HSolverLCAO::out_mat_t = 0;
int hsolver::HSolverLCAO::out_mat_dh = 0;
int Local_Orbital_Charge::out_dm = 0;
int Local_Orbital_Charge::out_dm1 = 0;
double ELEC_evolve::td_force_dt;
int ELEC_evolve::td_val_elec_01;
int ELEC_evolve::td_val_elec_02;
int ELEC_evolve::td_val_elec_03;
bool ELEC_evolve::td_vext;
std::vector<int> ELEC_evolve::td_vext_dire_case;
bool ELEC_evolve::out_dipole;
bool ELEC_evolve::out_efield;
double ELEC_evolve::td_print_eij;
int ELEC_evolve::td_edm;
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
double Force_Stress_LCAO::force_invalid_threshold_ev = 0.0;
double BFGS_Basic::relax_bfgs_w1 = -1.0;
double BFGS_Basic::relax_bfgs_w2 = -1.0;
double Ions_Move_Basic::relax_bfgs_rmax = -1.0;
double Ions_Move_Basic::relax_bfgs_rmin = -1.0;
double Ions_Move_Basic::relax_bfgs_init = -1.0;
int Ions_Move_Basic::out_stru=0;
int Lattice_Change_Basic::out_stru = 0;
double Ions_Move_CG::RELAX_CG_THR =-1.0;
int ModuleSymmetry::Symmetry::symm_flag=0;

Charge_Mixing::Charge_Mixing(){}
Charge_Mixing::~Charge_Mixing(){}
pseudopot_cell_vnl::pseudopot_cell_vnl(){}
pseudopot_cell_vnl::~pseudopot_cell_vnl(){}
pseudopot_cell_vl::pseudopot_cell_vl(){}
pseudopot_cell_vl::~pseudopot_cell_vl(){}
ORB_gaunt_table::ORB_gaunt_table(){}
ORB_gaunt_table::~ORB_gaunt_table(){}
ModuleDFTU::DFTU::DFTU(){}
ModuleDFTU::DFTU::~DFTU(){}
energy::energy(){}
energy::~energy(){}
Structure_Factor::Structure_Factor(){}
Structure_Factor::~Structure_Factor(){}
ModuleSymmetry::Symmetry::Symmetry(){}
ModuleSymmetry::Symmetry::~Symmetry(){}
ModuleSymmetry::Symmetry_Basic::Symmetry_Basic(){}
ModuleSymmetry::Symmetry_Basic::~Symmetry_Basic(){}
WF_igk::WF_igk(){}
WF_igk::~WF_igk(){}
WF_atomic::WF_atomic(){}
WF_atomic::~WF_atomic(){}
wavefunc::wavefunc(){}
wavefunc::~wavefunc(){}
UnitCell::UnitCell(){}
UnitCell::~UnitCell(){}
#ifdef __LCAO
InfoNonlocal::InfoNonlocal(){}
InfoNonlocal::~InfoNonlocal(){}
LCAO_Orbitals::LCAO_Orbitals(){}
LCAO_Orbitals::~LCAO_Orbitals(){}
#endif
Magnetism::Magnetism(){}
Magnetism::~Magnetism(){}

void Charge_Mixing::set_mixing(
    const std::string &mixing_mode_in,
    const double &mixing_beta_in,
    const int &mixing_ndim_in,
    const double &mixing_gg0_in,
    const bool &mixing_tau_in
){return;}
void Charge_Mixing::need_auto_set(){}
void Occupy::decision(const std::string &name,const std::string &smearing_method,const double &smearing_sigma){return;}
void UnitCell::setup(const std::string&,const int&,const int&,const bool&,const std::string&){return;}
void Structure_Factor::set(const int&){return;}

namespace GlobalC
{
	UnitCell ucell;
	wavefunc wf;
	ModuleSymmetry::Symmetry symm;
	Structure_Factor sf;
	energy en;
	ModuleDFTU::DFTU dftu;
	Restart restart;
	pseudopot_cell_vnl ppcell;
	Charge_Mixing CHR_MIX;
}

#undef private
