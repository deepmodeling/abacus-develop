#ifndef LCAO_MATRIX_H
#define LCAO_MATRIX_H

#include "../src_pw/tools.h"
#include "src_parallel/parallel_orbitals.h"

class LCAO_Matrix 
{
	friend class energy;
	friend class Mulliken_Charge;

	public:

	LCAO_Matrix();
	~LCAO_Matrix();

	void divide_HS_in_frag(const bool isGamma, Parallel_Orbitals &po);

	private:
	
	void allocate_HS_gamma(const long &nloc);

	void allocate_HS_k(const long &nloc);


	public:
	//------------------------------
	// H, S, Hfixed 
	// used in gamma only algorithm.
	// thse matrix are used to
	// diagonalize.
	//------------------------------
	double* Hloc;
	double* Sloc;
	double* Hloc_fixed;
	double* Sdiag; // used in pdiag_double.cpp
	
	//------------------------------
	// 1. Hamiltonian(vl),
	// 2. overlap matrix Sloc2
	// 3. fixed (vna+T+Vnl) matrix.
	// used in kpoint algorithm.
	// these matrix are used to
	// diagonalize.
	//------------------------------
	complex<double> *Hloc2;
	complex<double> *Sloc2;
	complex<double> *Hloc_fixed2;
	complex<double> *Sdiag2; // used in pdiag_double.cpp
	//with soc, zhengdy-soc
/*	ComplexMatrix Hloc2_soc;
	ComplexMatrix Sloc2_soc;
	ComplexMatrix Hloc_fixed2_soc;
	ComplexMatrix Sdiag2_soc;*/


	//------------------------------
	// Store H(mu,nu')
	// nu' : nu in near unitcell R.
	// used in kpoint algorithm.
	// these matrixed are used
	// for 'folding_matrix' in lcao_nnr,
	// HlocR -> Hloc2,
	// SlocR -> Sloc2, 
	//------------------------------
	double* HlocR;
	double* SlocR;
	double* Hloc_fixedR;
	//with soc, zhengdy-soc
	complex<double>* HlocR_soc;
	complex<double>* SlocR_soc;
	complex<double>* Hloc_fixedR_soc;
	//LiuXh add 2019-07-15
	double ****Hloc_fixedR_tr;
	double ****SlocR_tr;
	double ****HR_tr;
	complex<double> ****Hloc_fixedR_tr_soc;
	complex<double> ****SlocR_tr_soc;
	complex<double> ****HR_tr_soc;	


	//========================================
	// FORCE
	//========================================

	//-----------------------------------------
	// force in LCAO
	// used in gamma only algorithm.
	//-----------------------------------------
	double* DSloc_x;
	double* DSloc_y;
	double* DSloc_z;

	//-----------------------------------------
	// force in LCAO
	// used in k-points algorithm.
	//-----------------------------------------
	double* DSloc_Rx;
	double* DSloc_Ry;
	double* DSloc_Rz;

	//-----------------------------------------
	// dT + part of dVNL
	// used in gamma only algorithm.
	//-----------------------------------------
	double* DHloc_fixed_x; 
	double* DHloc_fixed_y;
	double* DHloc_fixed_z;

	//-----------------------------------------
	// dT + part of dVNL
	// used in kpoint algorithm.
	//-----------------------------------------
	double* DHloc_fixedR_x;
	double* DHloc_fixedR_y;
	double* DHloc_fixedR_z;

	//----------------------------------------
	// r_mu - r_nu
	//----------------------------------------
	double* DH_r;//zhengdy added 2017-07                        
	double* stvnl11;
	double* stvnl12;
	double* stvnl13;
	double* stvnl22;
	double* stvnl23;
	double* stvnl33;
	double* DSloc_11;
	double* DSloc_12;
	double* DSloc_13;
	double* DSloc_22;
	double* DSloc_23;
	double* DSloc_33;
	double* DHloc_fixed_11;
	double* DHloc_fixed_12;
	double* DHloc_fixed_13;
	double* DHloc_fixed_22;
	double* DHloc_fixed_23;
	double* DHloc_fixed_33;


	void set_HSgamma(const int &iw1_all, const int &iw2_all, const double &v, const char &dtype);
	void set_HSk(const int &iw1_all, const int &iw2_all, const complex<double> &v, const char &dtype, const int spin = 0);

	void set_force (const int& iw1_all, const int& iw2_all, const double& vx, const double& vy, 
		const double& vz, const char &dtype);
	void set_stress (const int& iw1_all, const int& iw2_all, const double& vx, const double& vy,
		const double& vz, const char &dtype, const Vector3<double> &dtau);

	void set_HR_tr(const int &Rx, const int &Ry, const int &Rz, const int &iw1_all, const int &iw2_all, const double &v);
	void set_HR_tr_soc(const int &Rx, const int &Ry, const int &Rz, 
		const int &iw1_all, const int &iw2_all, const complex<double> &v); //LiuXh add 2019-07-16

	void zeros_HSgamma(const char &mtype);
	void zeros_HSk(const char &mtype);
	void zeros_HSR(const char &mtype, const int &nnr);

	void print_HSgamma(const char &mtype, ostream &os=cout);
	void print_HSk(const char &mtype, const char &vtype = 'C', const double &accuracy = 1.0e-5, ostream &os=cout);
	void update_Hloc(void);
	void update_Hloc2(void);

	void allocate_HS_R(const int &nnr);

	void output_HSk(const char &mtype, string &fn);
	//LiuXh add 2019-07-15
	void allocate_Hloc_fixedR_tr(void);
	void allocate_HR_tr(void);
	void allocate_SlocR_tr(void);
	void destroy_Hloc_fixedR_tr(void);

};

#endif
