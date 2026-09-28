#ifndef FORCE_STRESS_LCAO_H
#define FORCE_STRESS_LCAO_H

#include <string>

#include "edm.h"
#include "source_base/global_function.h"
#include "source_base/matrix.h"
#include "source_pw/module_pwdft/force_pw.h"
#include "source_pw/module_pwdft/stress_func.h"
#include "source_pw/module_pwdft/stru_fac.h"
#include "source_io/module_parameter/input_conv.h"
#include "source_psi/psi.h"
#ifdef __EXX
#include "source_lcao/module_ri/exx_lri_interface.h"
#endif
#include "force_stress_arrays.h"
#include "source_lcao/setup_exx.h" // for exx, mohan add 20251008
#include "source_lcao/setup_deepks.h" // for deepks, mohan add 20251010
#include "source_lcao/allocate_dm.h" // mohan add 2025-11-03
#include "source_pw/module_pwdft/dftu_base.h" // mohan add 2025-11-07
#include "source_hamilt/hamilt.h"
#include "source_hamilt/module_xc/exx_info.h"

namespace vdw
{
struct VdwResult;
}

class TwoCenterBundle;

// INPUT scalars that steer the force/stress kernels. Bundling them into one
// aggregate keeps getForceStress and its helpers from each re-reading the
// global PARAM object, and collapses five arguments into a single reference.
struct FSCalcConfig
{
    int nspin;
    int nbands;
    bool t_in_h;
    bool sc_mag_switch;
    std::string device;
};

// Force/stress component matrices assembled by getForceStress. Grouping them
// into a struct lets the assembly/print helpers take one reference instead of
// ~19 individual matrix arguments. Members are default-constructed and only
// created (allocated) when the corresponding term is active.
struct LCAOForceParts
{
    ModuleBase::matrix foverlap;
    ModuleBase::matrix ftvnl_dphi;
    ModuleBase::matrix fvnl_dbeta;
    ModuleBase::matrix fvl_dphi;
    ModuleBase::matrix fvl_dvl;
    ModuleBase::matrix fewalds;
    ModuleBase::matrix fcc;
    ModuleBase::matrix fscc;
    ModuleBase::matrix fvnl_dalpha; // deepks
    ModuleBase::matrix fpothybrid;
    ModuleBase::matrix force_u;
    ModuleBase::matrix force_dspin;
    ModuleBase::matrix force_exx;
    ModuleBase::matrix force_vdw;
    ModuleBase::matrix fefield;
    ModuleBase::matrix fefield_tddft;
    ModuleBase::matrix fgate;
    ModuleBase::matrix fsol;
};

struct LCAOStressParts
{
    ModuleBase::matrix soverlap;
    ModuleBase::matrix stvnl_dphi;
    ModuleBase::matrix svnl_dbeta;
    ModuleBase::matrix svl_dphi;
    ModuleBase::matrix sigmadvl;
    ModuleBase::matrix sigmahar;
    ModuleBase::matrix sigmaewa;
    ModuleBase::matrix sigmacc;
    ModuleBase::matrix sigmaxc;
    ModuleBase::matrix svnl_dalpha; // deepks
    ModuleBase::matrix stress_u;
    ModuleBase::matrix stress_dspin;
    ModuleBase::matrix stress_exx;
    ModuleBase::matrix stress_vdw;
};


template <typename T>
class Force_Stress_LCAO
{
    // mohan add 2021-02-09
    friend class md;
    friend void Input_Conv::Convert();
    friend class ions;

  public:
    Force_Stress_LCAO(Record_adj& ra, const int nat_in);
    ~Force_Stress_LCAO();

    void getForceStress(UnitCell& ucell,
                        const vdw::VdwResult* vdw_result,
                        const bool isforce,
                        const bool isstress,
                        const bool istestf,
                        const bool istests,
                        const Grid_Driver& gd,
                        Parallel_Orbitals& pv,
                        const elecstate::ElecState* pelec,
                        module_dm::Setup_DM<T> &dmat, // mohan add 2025-11-03
                        const psi::Psi<T>* psi,
                        const TwoCenterBundle& two_center_bundle,
                        const LCAO_Orbitals& orb,
                        ModuleBase::matrix& fcs,
                        ModuleBase::matrix& scs,
                        const pseudopot_cell_vl& locpp,
                        const Structure_Factor& sf,
                        const K_Vectors& kv,
                        ModulePW::PW_Basis* rhopw,
						surchem& solvent,
						Plus_U_Base &dftu, // mohan add 2025-11-07
                        Setup_DeePKS<T> &deepks,
                        Exx_NAO<T> &exx_nao,
                        ModuleSymmetry::Symmetry* symm,
                        const Exx_Info& exx_info,
                        const FSCalcConfig& cfg,
                        const int td_stype = 0,
                        hamilt::Hamilt<T>* p_hamilt = nullptr);

  private:
    int nat;
    Record_adj* RA = nullptr;
    CalEDM<T> edm_cal;
    Stress_Func<double> sc_pw;

    // Operator-based force/stress terms: kinetic, overlap, nonlocal,
    // rt-TDDFT hybrid gauge, local-potential Pulay term, and DeltaSpin.
    void cal_operator_fs(UnitCell& ucell,
                         const Grid_Driver& gd,
                         Parallel_Orbitals& pv,
                         const elecstate::ElecState* pelec,
                         module_dm::Setup_DM<T>& dmat,
                         const psi::Psi<T>* psi,
                         const TwoCenterBundle& two_center_bundle,
                         const LCAO_Orbitals& orb,
                         const K_Vectors& kv,
                         const bool isforce,
                         const bool isstress,
                         const FSCalcConfig& cfg,
                         const int td_stype,
                         hamilt::Hamilt<T>* p_hamilt,
                         LCAOForceParts& parts,
                         LCAOStressParts& sparts);


    // Local pseudopotential, Ewald, core-correction and self-consistent-field
    // force contributions, computed with the plane-wave Forces driver. Kept as
    // a member because it needs friend access to Forces::cal_force_*.
    void calForcePwPart(UnitCell& ucell,
                        ModuleBase::matrix& fvl_dvl,
                        ModuleBase::matrix& fewalds,
                        ModuleBase::matrix& fcc,
                        ModuleBase::matrix& fscc,
                        const double& etxc,
                        const ModuleBase::matrix& vnew,
                        const bool vnew_exist,
                        const Charge* const chr,
                        ModulePW::PW_Basis* rhopw,
                        const pseudopot_cell_vl& locpp,
                        const Structure_Factor& sf,
                        const std::string& device);

    static double force_invalid_threshold_ev;
};

template <typename T>
double Force_Stress_LCAO<T>::force_invalid_threshold_ev = 0.00;

// only for DFT+U, mohan add 2025-11-04
template <typename T>
void assign_dmk_ptr(
    module_dm::DensityMatrix<T,double>* dm,
    std::vector<std::vector<double>>*& dmk_d,
    std::vector<std::vector<std::complex<double>>>*& dmk_c
);

#endif
