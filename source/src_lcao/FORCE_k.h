#ifndef FORCE_LCAO_K_H
#define FORCE_LCAO_K_H

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../module_base/matrix.h"
#include "LCAO_matrix.h" 
#include "FORCE_gamma.h"
#include "src_lcao/local_orbital_charge.h"
#include "src_lcao/LCAO_hamilt.h"

class Force_LCAO_k : public Force_LCAO_gamma
{
	public :
	
	friend class Force_Stress_LCAO;

	Force_LCAO_k ();
	~Force_LCAO_k ();

private:
    
    LCAO_Hamilt* UHM;
    const Parallel_Orbitals* ParaV;

	//orthonormal force + contribution from T and VNL
	void ftable_k (
		const bool isforce,
        const bool isstress,
        Record_adj &ra, 
        std::vector<ModuleBase::ComplexMatrix>& wfc_k,
        Local_Orbital_Charge &loc,
        ModuleBase::matrix& foverlap,
		ModuleBase::matrix& ftvnl_dphi,
		ModuleBase::matrix& fvnl_dbeta,	
		ModuleBase::matrix& fvl_dphi,
		ModuleBase::matrix& soverlap,
		ModuleBase::matrix& stvnl_dphi,
		ModuleBase::matrix& svnl_dbeta,
#ifdef __DEEPKS
		ModuleBase::matrix& svl_dphi,
		ModuleBase::matrix& svnl_dalpha,
#else
		ModuleBase::matrix& svl_dphi,
#endif
		LCAO_Hamilt &uhm);

	// get the ds, dt, dvnl.
	void allocate_k(const Parallel_Orbitals &pv);

	void finish_k(void);
	
	// calculate the force due to < dphi | beta > < beta | phi >
    void cal_ftvnl_dphi_k(double** dm2d, const bool isforce, const bool isstress, Record_adj& ra,
        ModuleBase::matrix& ftvnl_dphi, ModuleBase::matrix& stvnl_dphi);

	// calculate the overlap force
    void cal_foverlap_k(const bool isforce, const bool isstress, Record_adj &ra, std::vector<ModuleBase::ComplexMatrix>& wfc_k,
        Local_Orbital_Charge& loc, ModuleBase::matrix& foverlap, ModuleBase::matrix& soverlap);

	// calculate the force due to < phi | Vlocal | dphi >
	void cal_fvl_dphi_k(double** dm2d, const bool isforce, const bool isstress, ModuleBase::matrix& fvl_dphi, ModuleBase::matrix& svl_dphi);
	void cal_fvl_dphi_k_new(const bool isforce, const bool isstress, ModuleBase::matrix& fvl_dphi, ModuleBase::matrix& svl_dphi);

	// old method to calculate the force due to < phi | dbeta > < beta | phi >
	void cal_fvnl_dbeta_k(double** dm2d, const bool isforce, const bool isstress, ModuleBase::matrix& fvnl_dbeta, ModuleBase::matrix& svnl_dbeta);
	// new method to calculate the force due to < phi | dbeta > < beta | phi > , developed by wenfei-li
	void cal_fvnl_dbeta_k_new(double** dm2d, const bool isforce, const bool isstress, ModuleBase::matrix& fvnl_dbeta, ModuleBase::matrix& svnl_dbeta);

	void test(double* mm, const std::string &name);

	// calculate the force due to < phi | dbeta > < beta | phi >
	void calFvnlDbeta(
		double** dm2d, 
		const bool &isforce, 
		const bool &isstress, 
		ModuleBase::matrix& fvnl_dbeta, 
		ModuleBase::matrix& svnl_dbeta,
		const int &vnl_method);

};
#endif
