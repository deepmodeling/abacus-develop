#ifndef FORCE_STRESS_LCAO_H
#define FORCE_STRESS_LCAO_H

#include "FORCE.h"
#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "source_base/matrix.h"
#include "source_pw/module_pwdft/forces.h"
#include "source_pw/module_pwdft/stress_func.h"
#include "source_pw/module_pwdft/structure_factor.h"
#include "source_io/input_conv.h"
#include "source_psi/psi.h"
#ifdef __EXX
#include "source_lcao/module_ri/Exx_LRI_interface.h"
#endif
#include "force_stress_arrays.h"
#include "source_lcao/setup_exx.h" // for exx, mohan add 20251008
#include "source_lcao/setup_deepks.h" // for deepks, mohan add 20251010
#include "source_lcao/setup_dm.h" // mohan add 2025-11-03


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
                        const bool isforce,
                        const bool isstress,
                        const bool istestf,
                        const bool istests,
                        const Grid_Driver& gd,
                        Parallel_Orbitals& pv,
                        const elecstate::ElecState* pelec,
                        LCAO_domain::Setup_DM<T> &dmat, // mohan add 2025-11-03
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
                        Setup_DeePKS<T> &deepks,
                        Exx_NAO<T> &exx_nao,
                        ModuleSymmetry::Symmetry* symm);

  private:
    int nat;
    Record_adj* RA;
    Force_LCAO<T> flk;
    Stress_Func<double> sc_pw;

    void forceSymmetry(const UnitCell& ucell, ModuleBase::matrix& fcs, ModuleSymmetry::Symmetry* symm);

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
                        const Structure_Factor& sf);

    void integral_part(const bool isGammaOnly,
                       const bool isforce,
                       const bool isstress,
                       const UnitCell& ucell,
                       const Grid_Driver& gd,
                       ForceStressArrays& fsr, // mohan add 2024-06-15
					   const elecstate::ElecState* pelec,
					   const elecstate::DensityMatrix<T, double>* dm, // mohan add 2025-11-04
					   const psi::Psi<T>* psi,
                       ModuleBase::matrix& foverlap,
                       ModuleBase::matrix& ftvnl_dphi,
                       ModuleBase::matrix& fvnl_dbeta,
                       ModuleBase::matrix& fvl_dphi,
                       ModuleBase::matrix& soverlap,
                       ModuleBase::matrix& stvnl_dphi,
                       ModuleBase::matrix& svnl_dbeta,
                       ModuleBase::matrix& svl_dphi,
                       ModuleBase::matrix& fvnl_dalpha,
                       ModuleBase::matrix& svnl_dalpha,
                       Setup_DeePKS<T>& deepks,
                       const TwoCenterBundle& two_center_bundle,
                       const LCAO_Orbitals& orb,
                       const Parallel_Orbitals& pv,
                       const K_Vectors& kv);

    void calStressPwPart(UnitCell& ucell,
                         ModuleBase::matrix& sigmadvl,
                         ModuleBase::matrix& sigmahar,
                         ModuleBase::matrix& sigmaewa,
                         ModuleBase::matrix& sigmacc,
                         ModuleBase::matrix& sigmaxc,
                         const double& etxc,
                         const Charge* const chr,
                         ModulePW::PW_Basis* rhopw,
                         const pseudopot_cell_vl& locpp,
                         const Structure_Factor& sf);

    static double force_invalid_threshold_ev;
};

template <typename T>
double Force_Stress_LCAO<T>::force_invalid_threshold_ev = 0.00;

// only for DFT+U, mohan add 2025-11-04
template <typename T>
void assign_dmk_ptr(
    elecstate::DensityMatrix<T,double>* dm,
    std::vector<std::vector<double>>*& dmk_d,
    std::vector<std::vector<std::complex<double>>>*& dmk_c,
    bool gamma_only_local
);

#endif
