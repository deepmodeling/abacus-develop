#ifndef FORCE_STRESS_LCAO_H
#define FORCE_STRESS_LCAO_H

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../module_base/matrix.h"
#include "FORCE_k.h"
//#include "./force_lcao_gamma.h"
#include "../src_pw/stress_func.h"
#include "../input_conv.h"
#include "../src_pw/forces.h"
#include "module_psi/psi.h"

class Force_Stress_LCAO  
{
	// mohan add 2021-02-09
	friend class md;
	friend class Run_MD_LCAO;
	friend void Input_Conv::Convert();
	friend class Update_input;
	friend class ions;
	friend class MD_func;

	public :
	
	Force_Stress_LCAO (Record_adj &ra);
    ~Force_Stress_LCAO();

	void getForceStress(
		const bool isforce, 
		const bool isstress, 
		const bool istestf, 
        const bool istests,
        Local_Orbital_Charge& loc,
        const psi::Psi<double>* psid,
		const psi::Psi<std::complex<double>>* psi,
        LCAO_Hamilt &uhm,
        ModuleBase::matrix& fcs,
		ModuleBase::matrix &scs);

private:
    
    Record_adj* RA;
	Force_LCAO_k flk;
//	Force_LCAO_gamma flg;
	Stress_Func sc_pw;
	Forces f_pw;
	
	void print_force(const std::string &name, ModuleBase::matrix& f, const bool screen, bool ry)const;
	void printforce_total (const bool ry, const bool istestf, ModuleBase::matrix& fcs);
	

	void forceSymmetry(ModuleBase::matrix &fcs);

	void calForcePwPart(
		ModuleBase::matrix &fvl_dvl, 
		ModuleBase::matrix &fewalds, 
		ModuleBase::matrix &fcc, 
		ModuleBase::matrix &fscc);

	void calForceStressIntegralPart(
		const bool isGammaOnly,
		const bool isforce,
        const bool isstress,
        Local_Orbital_Charge& loc,
        const psi::Psi<double>* psid,
		const psi::Psi<std::complex<double>>* psi,
        ModuleBase::matrix& foverlap,
		ModuleBase::matrix &ftvnl_dphi,
		ModuleBase::matrix &fvnl_dbeta,	
		ModuleBase::matrix &fvl_dphi,
		ModuleBase::matrix &soverlap,
		ModuleBase::matrix &stvnl_dphi,
		ModuleBase::matrix &svnl_dbeta,
#if __DEEPKS
		ModuleBase::matrix& svl_dphi,
	    ModuleBase::matrix& svnl_dalpha,
#else
	    ModuleBase::matrix& svl_dphi,
#endif
        LCAO_Hamilt &uhm);
    


	void calStressPwPart(
		ModuleBase::matrix &sigmadvl,
		ModuleBase::matrix &sigmahar,
		ModuleBase::matrix &sigmaewa,
		ModuleBase::matrix &sigmacc,
		ModuleBase::matrix &sigmaxc);
	
	static double force_invalid_threshold_ev;
	static double output_acc; // control the accuracy
};
#endif
