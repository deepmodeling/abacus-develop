#ifndef FORCE_LCAO_GAMMA_H
#define FORCE_LCAO_GAMMA_H

#include "../src_pw/tools.h"
#include "LCAO_matrix.h" 

class Force_LCAO_gamma
{
	public :

	friend class Force_Stress_LCAO;

	Force_LCAO_gamma ();
	~Force_LCAO_gamma ();

	private:
	
	//orthonormal force + contribution from T and VNL
	void ftable_gamma (
		const bool isforce,
		const bool isstress,
		matrix& foverlap,
		matrix& ftvnl_dphi,
		matrix& fvnl_dbeta,	
		matrix& fvl_dphi,
		matrix& soverlap,
		matrix& stvnl_dphi,
		matrix& svnl_dbeta,
		matrix& svl_dphi
		);

	// get the ds, dt, dvnl.
	void allocate_gamma(void);

	void finish_ftable_gamma(void);

	void average_force(double* fm);

	void test_gamma(double* mm, const string &name);


	//-------------------------------------------------------------
	// forces reated to overlap matrix
	// forces related to energy density matrix 
	//-------------------------------------------------------------
	void set_EDM_gamma(matrix& dm, bool with_energy);

	double set_EDM_element(const int &ii, const int &jj, const bool with_energy, 
	 double*** coef1, double*** coef2, const int &is);

	void cal_foverlap(
		const bool isforce, 
		const bool isstress, 
		matrix& foverlap, 
		matrix& soverlap);	

	//-------------------------------------------------------------
	// forces related to kinetic and non-local pseudopotentials
	//--------------------------------------------------------------
	void cal_ftvnl_dphi(
		matrix& dm2d, 
		const bool isforce, 
		const bool isstress, 
		matrix& ftvnl_dphi, 
		matrix& stvnl_dphi);

	void cal_ftvnl_dphi(
		const std::vector<matrix> &dm2d, 
		const bool isforce, 
		const bool isstress, 
		matrix& ftvnl_dphi, 
		matrix& stvnl_dphi);

	void cal_fvnl_dbeta(
		matrix& dm2d, 
		const bool isforce, 
		const bool isstress, 
		matrix& fvnl_dbeta, 
		matrix& svnl_dbeta);

	void cal_fvnl_dbeta(
		const std::vector<matrix> &dm2d, 
		const bool isforce, 
		const bool isstress, 
		matrix& fvnl_dbeta, 
		matrix& svnl_dbeta);

	//-------------------------------------------
	// forces related to local pseudopotentials
	//-------------------------------------------
	void cal_fvl_dphi(
		matrix& dm2d, 
		const bool isforce, 
		const bool isstress, 
		matrix& fvl_dphi, 
		matrix& svl_dphi);

	void cal_fvl_dphi(
		const std::vector<matrix> &dm2d, 
		const bool isforce, 
		const bool isstress, 
		matrix& fvl_dphi, 
		matrix& svl_dphi);

};
#endif
