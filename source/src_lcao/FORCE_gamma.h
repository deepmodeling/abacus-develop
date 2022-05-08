#ifndef FORCE_LCAO_GAMMA_H
#define FORCE_LCAO_GAMMA_H

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../module_base/matrix.h"
#include "LCAO_matrix.h" 
#include "src_lcao/local_orbital_charge.h"
#include "src_lcao/LCAO_hamilt.h"

class Force_LCAO_gamma
{
	public :

	friend class Force_Stress_LCAO;

	Force_LCAO_gamma ();
	~Force_LCAO_gamma ();

private:

    LCAO_Hamilt *UHM;
    const Parallel_Orbitals* ParaV;

	//orthonormal force + contribution from T and VNL
	void ftable_gamma (
		const bool isforce,
        const bool isstress,
        vector<ModuleBase::matrix>& wfc_gamma,
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
	void allocate_gamma(const Parallel_Orbitals &pv);

	void finish_ftable_gamma(void);

	void average_force(double* fm);

	void test_gamma(double* mm, const std::string &name);


	//-------------------------------------------------------------
	// forces reated to overlap matrix
	// forces related to energy density matrix 
	//-------------------------------------------------------------

	void cal_foverlap(
		const bool isforce, 
        const bool isstress,
        vector<ModuleBase::matrix>& wfc_gamma,
        Local_Orbital_Charge &loc,
        ModuleBase::matrix& foverlap,
		ModuleBase::matrix& soverlap);	

	//-------------------------------------------------------------
	// forces related to kinetic and non-local pseudopotentials
	//--------------------------------------------------------------
	void cal_ftvnl_dphi(
		ModuleBase::matrix& dm2d, 
		const bool isforce, 
		const bool isstress, 
		ModuleBase::matrix& ftvnl_dphi, 
		ModuleBase::matrix& stvnl_dphi);

	void cal_ftvnl_dphi(
		const std::vector<ModuleBase::matrix> &dm2d, 
		const bool isforce, 
		const bool isstress, 
		ModuleBase::matrix& ftvnl_dphi, 
		ModuleBase::matrix& stvnl_dphi);

	void cal_fvnl_dbeta(
		ModuleBase::matrix& dm2d, 
		const bool isforce, 
		const bool isstress, 
		ModuleBase::matrix& fvnl_dbeta, 
		ModuleBase::matrix& svnl_dbeta);

	void cal_fvnl_dbeta(
		const std::vector<ModuleBase::matrix> &dm2d, 
		const bool isforce, 
		const bool isstress, 
		ModuleBase::matrix& fvnl_dbeta, 
		ModuleBase::matrix& svnl_dbeta);

	void cal_fvnl_dbeta_new(
		ModuleBase::matrix& dm2d, 
		const bool isforce, 
		const bool isstress, 
		ModuleBase::matrix& fvnl_dbeta, 
		ModuleBase::matrix& svnl_dbeta);

	void cal_fvnl_dbeta_new(
		const std::vector<ModuleBase::matrix> &dm2d, 
		const bool isforce, 
		const bool isstress, 
		ModuleBase::matrix& fvnl_dbeta, 
		ModuleBase::matrix& svnl_dbeta);

	//-------------------------------------------
	// forces related to local pseudopotentials
	//-------------------------------------------
	void cal_fvl_dphi(
		ModuleBase::matrix& dm2d, 
		const bool isforce, 
        const bool isstress,
        ModuleBase::matrix& fvl_dphi,
		ModuleBase::matrix& svl_dphi);

	void cal_fvl_dphi(
		const std::vector<ModuleBase::matrix> &dm2d, 
		const bool isforce, 
        const bool isstress,
        ModuleBase::matrix& fvl_dphi,
		ModuleBase::matrix& svl_dphi);

	void cal_fvl_dphi_new(
		double*** DM_in,
		const bool isforce, 
        const bool isstress,
        ModuleBase::matrix& fvl_dphi,
		ModuleBase::matrix& svl_dphi);

	void calFvnlDbeta(
		const std::vector<ModuleBase::matrix> &dm2d, 
		const bool &isforce, 
		const bool &isstress, 
		ModuleBase::matrix& fvnl_dbeta, 
		ModuleBase::matrix& svnl_dbeta,
		const int &vnl_method);

	public:
	// calculate dVnl=<phi|dVnl|dphi> in LCAO
    void NonlocalDphi(const int& nspin, const int& vnl_method, const bool& cal_deri,
        LCAO_gen_fixedH &genH);
};

// this namespace used to store global function for some stress operation
namespace StressTools{
	//set upper matrix to whole matrix
	void stress_fill( 
		const double& lat0_, 
		const double& omega_,
		ModuleBase::matrix& stress_matrix);
}
#endif
