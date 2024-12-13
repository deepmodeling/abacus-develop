#ifndef FORCE_STRESS_LCAO_H
#define FORCE_STRESS_LCAO_H

#include "FORCE.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/matrix.h"
#include "module_hamilt_pw/hamilt_pwdft/forces.h"
#include "module_hamilt_pw/hamilt_pwdft/stress_func.h"
#include "module_hamilt_pw/hamilt_pwdft/structure_factor.h"
#include "module_io/input_conv.h"
#include "module_psi/psi.h"
#ifdef __EXX
#include "module_ri/Exx_LRI.h"
#endif
#include "force_stress_arrays.h"
#include "module_hamilt_lcao/module_gint/gint_gamma.h"
#include "module_hamilt_lcao/module_gint/gint_k.h"

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

    void getForceStress(const bool isforce,
                        const bool isstress,
                        const bool istestf,
                        const bool istests,
                        const UnitCell& ucell,
                        const Grid_Driver& gd,
                        Parallel_Orbitals& pv,
                        const elecstate::ElecState* pelec,
                        const psi::Psi<T>* psi,
                        Gint_Gamma& gint_gamma, // mohan add 2024-04-01
                        Gint_k& gint_k,         // mohan add 2024-04-01
                        const TwoCenterBundle& two_center_bundle,
                        const LCAO_Orbitals& orb,
                        ModuleBase::matrix& fcs,
                        ModuleBase::matrix& scs,
                        const pseudopot_cell_vnl& nlpp,
                        const Structure_Factor& sf,
                        const K_Vectors& kv,
                        ModulePW::PW_Basis* rhopw,
#ifdef __EXX
                        Exx_LRI<double>& exx_lri_double,
                        Exx_LRI<std::complex<double>>& exx_lri_complex,
#endif
                        ModuleSymmetry::Symmetry* symm);

  private:
    int nat;
    Record_adj* RA;
    Force_LCAO<T> flk;
    Stress_Func<double> sc_pw;
    Forces<double> f_pw;

    void forceSymmetry(const UnitCell& ucell, 
                       ModuleBase::matrix& fcs, 
                       ModuleSymmetry::Symmetry* symm);

    void calForcePwPart(const UnitCell& ucell,
                        ModuleBase::matrix& fvl_dvl,
                        ModuleBase::matrix& fewalds,
                        ModuleBase::matrix& fcc,
                        ModuleBase::matrix& fscc,
                        const double& etxc,
                        const ModuleBase::matrix& vnew,
                        const bool vnew_exist,
                        const Charge* const chr,
                        ModulePW::PW_Basis* rhopw,
                        const pseudopot_cell_vnl& nlpp,
                        const Structure_Factor& sf);

    void integral_part(const bool isGammaOnly,
                       const bool isforce,
                       const bool isstress,
                       const UnitCell& ucell,
                       const Grid_Driver& gd,
                       ForceStressArrays& fsr, // mohan add 2024-06-15
                       const elecstate::ElecState* pelec,
                       const psi::Psi<T>* psi,
                       ModuleBase::matrix& foverlap,
                       ModuleBase::matrix& ftvnl_dphi,
                       ModuleBase::matrix& fvnl_dbeta,
                       ModuleBase::matrix& fvl_dphi,
                       ModuleBase::matrix& soverlap,
                       ModuleBase::matrix& stvnl_dphi,
                       ModuleBase::matrix& svnl_dbeta,
                       ModuleBase::matrix& svl_dphi,
#if __DEEPKS
                       ModuleBase::matrix& svnl_dalpha,
#endif
                       Gint_Gamma& gint_gamma,
                       Gint_k& gint_k,
                       const TwoCenterBundle& two_center_bundle,
                       const LCAO_Orbitals& orb,
                       const Parallel_Orbitals& pv,
                       const K_Vectors& kv);

    void calStressPwPart(const UnitCell& ucell,
                         ModuleBase::matrix& sigmadvl,
                         ModuleBase::matrix& sigmahar,
                         ModuleBase::matrix& sigmaewa,
                         ModuleBase::matrix& sigmacc,
                         ModuleBase::matrix& sigmaxc,
                         const double& etxc,
                         const Charge* const chr,
                         ModulePW::PW_Basis* rhopw,
                         const pseudopot_cell_vnl& nlpp,
                         const Structure_Factor& sf);

    static double force_invalid_threshold_ev;
};

template <typename T>
double Force_Stress_LCAO<T>::force_invalid_threshold_ev = 0.00;

#endif
