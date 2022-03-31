#ifndef CHARGE_PULAY_H
#define CHARGE_PULAY_H

//==================================
// (1) Plain Mixing
// (2) KerKer Mixing
// (3) Pulay Mixing
// (4) Modified Broden Mixing
//===================================
#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../module_base/matrix.h"
#include "charge_mixing.h"

#ifdef __EXX // mohan add 2021-05-23
#include "../src_ri/exx_abfs-dm.h"
#include "../src_ri/exx_abfs-parallel-communicate-hexx.h"
#endif

class Charge_Pulay: public Charge_Mixing
{
	public:
	Charge_Pulay();
	~Charge_Pulay();

	int irstep; //mohan add 2012-02-10
	int idstep;
	int totstep;
	int rstep; // the record step;
	int dstep; // Delta step " dstep = rstep-1 ".
	double* alpha; // - sum (Abar * dRR)
	
	// Peize Lin add 2018-11-01, and set new_e_iteration protected
	const bool &get_new_e_iteration(){ return new_e_iteration; }
	void set_new_e_iteration( const bool new_e_iteration_in );			

	// mohan add 2010-07-16
	bool new_e_iteration;

	protected:

	
	// Pulay mixing method.
	void Pulay_mixing();
	double*** Rrho;// Rrho(i) = rho(i) - rho_save(i), (GlobalV::NSPIN, rstep, pw.nrxx)
	double*** dRrho;// dRrho(i) = Rrho(i+1) - Rrho(i), (GlobalV::NSPIN, dstep, pw.nrxx)
	double*** drho;// drho(i)= rho_save(i+1) - rho_save2(i), (GlobalV::NSPIN, dstep, pw.nrxx)
	double** rho_save2;//rho_save: rho_in, rho_save2: rho_in(last step)
	bool initp; // p stands for pulay algorithms
	std::complex<double>*** dF; // dF(i) = rhog(i) - rhog_save(i), (GlobalV::NSPIN, rstep, pw.ngmc)
	std::complex<double>*** dn; // dn(i) = rhog(i+1) - rhog(i), (GlobalV::NSPIN, rstep, pw.ngmc)
	
	ModuleBase::matrix Abar; // <dR_j|dR_i>^{-1}
	double* dRR; // <dR_j|R_m>
	
	void allocate_pulay(const int &scheme);
	void generate_datas(const int &irstep, const int &idstep, const int &totstep);
	void generate_Abar(const int &scheme, ModuleBase::matrix &A)const;
	void inverse_preA(const int &dim, ModuleBase::matrix &preA)const;
	void inverse_real_symmetry_matrix(const int &scheme, ModuleBase::matrix &A)const; // indicate the spin.
	void generate_dRR(const int &m);
	void generate_alpha(const int &scheme);
	void generate_new_rho(const int &is,const int &m);

	void generate_residual_vector(double *residual, const double* rho_out, const double* rho_in)const;
	double calculate_residual_norm(double *residual1, double *residual2)const;

	// Sophisticated mixing method.
	void Modified_Broyden_mixing();

#ifdef __EXX
    friend class Exx_Abfs::DM;
    #ifdef __MPI
    friend class Exx_Abfs::Parallel::Communicate::Hexx;
    #endif
#endif

};

#endif
