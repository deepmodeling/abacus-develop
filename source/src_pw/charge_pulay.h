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

class Charge_Pulay: public Charge_Mixing
{
	public:
	Charge_Pulay();
	~Charge_Pulay();
	
	void reset(const bool final_scf); //if first electronic step, then reset charge mixing

	// mohan add 2010-07-16
	bool new_e_iteration;

	// normally these parameters will not be used
	// outside charge mixing, but Exx is using them
	int get_totstep() const {return totstep;}
	int get_rstep() const {return rstep;}
	int get_dstep() const {return dstep;}
	int get_idstep() const {return idstep;}
	double* get_alpha() const {return alpha;}

	protected:

	int irstep; //mohan add 2012-02-10
	int idstep;
	int totstep;
	int rstep; // the record step;
	int dstep; // Delta step " dstep = rstep-1 ".
	double* alpha; // - sum (Abar * dRR)

	// Pulay mixing method.
	void Pulay_mixing(double** rho, double**rho_save);
	double*** Rrho;// Rrho(i) = rho(i) - rho_save(i), (GlobalV::NSPIN, rstep, pw.nrxx)
	double*** dRrho;// dRrho(i) = Rrho(i+1) - Rrho(i), (GlobalV::NSPIN, dstep, pw.nrxx)
	double*** drho;// drho(i)= rho_save(i+1) - rho_save2(i), (GlobalV::NSPIN, dstep, pw.nrxx)
	double** rho_save2;//rho_save: rho_in, rho_save2: rho_in(last step)
	bool initp; // p stands for pulay algorithms
	std::complex<double>*** dF; // dF(i) = rhog(i) - rhog_save(i), (GlobalV::NSPIN, rstep, rhopw->npw)
	std::complex<double>*** dn; // dn(i) = rhog(i+1) - rhog(i), (GlobalV::NSPIN, rstep, rhopw->npw)
	
	ModuleBase::matrix Abar; // <dR_j|dR_i>^{-1}
	double* dRR; // <dR_j|R_m>
	
	void allocate_pulay(const int &scheme);
	void generate_datas(const int &irstep, const int &idstep, const int &totstep, double** rho, double** rho_save);
	void generate_Abar(const int &scheme, ModuleBase::matrix &A)const;
	void inverse_preA(const int &dim, ModuleBase::matrix &preA)const;
	void inverse_real_symmetry_matrix(const int &scheme, ModuleBase::matrix &A)const; // indicate the spin.
	void generate_dRR(const int &m);
	void generate_alpha(const int &scheme);
	void generate_new_rho(const int &is,const int &m, double**rho, double** rho_save);

	void generate_residual_vector(double *residual, const double* rho_out, const double* rho_in)const;
	double calculate_residual_norm(double *residual1, double *residual2)const;

};

#endif
