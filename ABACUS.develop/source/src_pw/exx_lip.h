//==========================================================
// AUTHOR : Peize Lin
// DATE : 2015-03-10
//==========================================================
#ifndef EXX_LIP_H
#define EXX_LIP_H

#include "module_base/complexmatrix.h"
#include "module_base/vector3.h"
#include "src_pw/hamilt_pw.h"
#include "src_pw/exx_global.h"

class K_Vectors;
class wavefunc;
class PW_Basis;
class Use_FFT;
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

	void init(K_Vectors *kv_ptr_in, wavefunc *wf_ptr_in, PW_Basis *pw_ptr_in, Use_FFT *UFFT_ptr_in, UnitCell_pseudo *ucell_ptr_in);
	void cal_exx();
	const complex<double> * const * const * get_exx_matrix() const { return exx_matrix; }
	double get_exx_energy() const { return exx_energy; }
	
	void write_q_pack() const;
	
private:

	bool init_finish;
	
	int gzero_rank_in_pool;
		
	struct k_package
	{
		K_Vectors *kv_ptr;
		wavefunc *wf_ptr;
		matrix wf_wg;
		ComplexMatrix *hvec_array;		
	} *k_pack, *q_pack;
	
	int iq_vecik;
	
	complex<double> **phi;
	complex<double> ***psi;
	double *recip_qkg2;
	double sum2_factor;
	complex<double> *b;
	complex<double> *b0;
	complex<double> *sum1;
	complex<double> **sum3;
	
	complex<double> ***exx_matrix;
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
		const ComplexMatrix &psi,
		ComplexMatrix &evc,
		double *en);	
	
public:

	PW_Basis *pw_ptr;
	Use_FFT *UFFT_ptr;
	UnitCell_pseudo *ucell_ptr;
};
	

#endif
