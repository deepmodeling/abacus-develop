#ifndef FORCE_STRESS_LCAO_H
#define FORCE_STRESS_LCAO_H

#include "FORCE_k.h"
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
#include "module_hamilt_lcao/module_gint/gint_gamma.h"
#include "module_hamilt_lcao/module_gint/gint_k.h"

template<typename T>
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
		Local_Orbital_Charge& loc,
		Parallel_Orbitals &pv,
		LCAO_Matrix &lm,
		const elecstate::ElecState* pelec,
        const psi::Psi<T>* psi,
		LCAO_Hamilt& uhm,
		Gint_Gamma &gint_gamma, // mohan add 2024-04-01
		Gint_k &gint_k, // mohan add 2024-04-01
        ModuleBase::matrix& fcs,
        ModuleBase::matrix& scs,
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
    Force_LCAO_k flk;
    //	Force_LCAO_gamma flg;
    Stress_Func<double> sc_pw;
    Forces<double> f_pw;

    void forceSymmetry(ModuleBase::matrix& fcs, ModuleSymmetry::Symmetry* symm);

    void calForcePwPart(ModuleBase::matrix& fvl_dvl,
                        ModuleBase::matrix& fewalds,
                        ModuleBase::matrix& fcc,
                        ModuleBase::matrix& fscc,
                        const double& etxc,
                        const ModuleBase::matrix& vnew,
                        const bool vnew_exist,
                        const Charge* const chr,
                        ModulePW::PW_Basis* rhopw,
                        const Structure_Factor& sf);

    void integral_part(
        const bool isGammaOnly,
        const bool isforce,
        const bool isstress,
        Local_Orbital_Charge& loc,
        const elecstate::ElecState* pelec,
        const psi::Psi<T>* psi,
        ModuleBase::matrix& foverlap,
        ModuleBase::matrix& ftvnl_dphi,
        ModuleBase::matrix& fvnl_dbeta,
        ModuleBase::matrix& fvl_dphi,
        ModuleBase::matrix& soverlap,
        ModuleBase::matrix& stvnl_dphi,
        ModuleBase::matrix& svnl_dbeta,
#if __DEEPKS
        ModuleBase::matrix& svl_dphi,
        ModuleBase::matrix& svnl_dalpha,
#else
        ModuleBase::matrix& svl_dphi,
#endif
		LCAO_Hamilt &uhm,
		Gint_Gamma &gint_gamma,
		Gint_k &gint_k,
	    Parallel_Orbitals &pv,
		LCAO_Matrix &lm,
		const K_Vectors& kv);

    void calStressPwPart(ModuleBase::matrix& sigmadvl,
                         ModuleBase::matrix& sigmahar,
                         ModuleBase::matrix& sigmaewa,
                         ModuleBase::matrix& sigmacc,
                         ModuleBase::matrix& sigmaxc,
                         const double& etxc,
                         const Charge* const chr,
                         ModulePW::PW_Basis* rhopw,
                         const Structure_Factor& sf);

    static double force_invalid_threshold_ev;
};

template<typename T>
double Force_Stress_LCAO<T>::force_invalid_threshold_ev = 0.00;

#endif
