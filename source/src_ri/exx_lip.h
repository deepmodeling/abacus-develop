//==========================================================
// AUTHOR : Peize Lin
// DATE : 2015-03-10
//==========================================================
#ifndef EXX_LIP_H
#define EXX_LIP_H

#include "../module_base/complexmatrix.h"
#include "../module_base/vector3.h"
#include "../src_pw/hamilt.h"
#include "../module_xc/exx_global.h"
#include "../module_pw/pw_basis_k.h"

class K_Vectors;
class wavefunc;
class UnitCell_pseudo;

class Exx_Lip
{
public:
	Exx_Lip( const Exx_Global::Exx_Info &info_global );
	~Exx_Lip();

	struct Exx_Info
	{
		const Exx_Global::Hybrid_Type &hybrid_type;

		const double &hse_omega;

		double lambda;

		Exx_Info( const Exx_Global::Exx_Info &info_global );
	};
	Exx_Info info;

	void init(K_Vectors *kv_ptr_in, wavefunc *wf_ptr_in, ModulePW::PW_Basis_K *wfc_basis_in, ModulePW::PW_Basis *rho_basis_in, UnitCell_pseudo *ucell_ptr_in);
	void cal_exx();
	const std::complex<double> * const * const * get_exx_matrix() const { return exx_matrix; }
	double get_exx_energy() const { return exx_energy; }

	void write_q_pack() const;

private:

	bool init_finish;

	int gzero_rank_in_pool;

	struct k_package
	{
		K_Vectors *kv_ptr;
		wavefunc *wf_ptr;
		ModuleBase::matrix wf_wg;
		ModuleBase::ComplexMatrix *hvec_array;		
	} *k_pack, *q_pack;

	int iq_vecik;

	std::complex<double> **phi;
	std::complex<double> ***psi;
	double *recip_qkg2;
	double sum2_factor;
	std::complex<double> *b;
	std::complex<double> *b0;
	std::complex<double> *sum1;
	std::complex<double> **sum3;

	std::complex<double> ***exx_matrix;
	double exx_energy;

	void wf_wg_cal();
	void phi_cal(k_package *kq_pack, int ikq);
	void psi_cal();
	void judge_singularity( int ik);
	void qkg2_exp(int ik, int iq);
	void b_cal(int ik, int iq, int ib);
	void sum3_cal(int iq, int ib);
	void b_sum(int iq, int ib);
	void sum_all(int ik);
	void exx_energy_cal();
	void read_q_pack();

// mohan comment out 2021-02-24
	friend void Hamilt_PW::diagH_subspace(
		const int ik,
		const int nstart,
		const int n_band,
		const ModuleBase::ComplexMatrix &psi,
		ModuleBase::ComplexMatrix &evc,
		double *en);	
	
public:

	ModulePW::PW_Basis *rho_basis;
	ModulePW::PW_Basis_K *wfc_basis;

	UnitCell_pseudo *ucell_ptr;
};


#endif
