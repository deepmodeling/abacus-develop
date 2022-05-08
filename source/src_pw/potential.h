#ifndef POTENTIAL_H
#define POTENTIAL_H

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../module_base/matrix.h"
#include "../module_base/complexmatrix.h"
#include "module_esolver/esolver_ks_lcao.h"

class Potential
{
	public:

	friend class Hamilt_PW; 
	friend class Electrons;
    friend class ELEC_scf;
    friend class ModuleESolver::ESolver_KS_LCAO;

    // constructor and deconstructor
    Potential();
    ~Potential();

    //==========================================================
    // init_chg : "atomic" or "file"
	// chg_extrap : extrapolation methods for potential
    // vr(nspin,ncxyz) : Hartree + xc potentials in real space
    // vr_eff(nspin,ncxyz) : effective potential in real space 
    // vnew(nspin,ncxyz) : V_out - V_in, needed in scf
	// vltot: the local potential in real space
	// out_pot: options to print out potentials 
    //==========================================================

    std::string init_chg;
	std::string chg_extrap;
    ModuleBase::matrix vr;
    ModuleBase::matrix vr_eff;
    ModuleBase::matrix vnew;
	
	ModuleBase::matrix vofk; //kinetic energy density, for meta-GGA; wenfei 2021-07-28

    double *vr_eff1; 
#ifdef __CUDA
	double *d_vr_eff1;
#endif
    double *vltot;
	int out_pot; // mohan add 2011-02-28

    void allocate(const int nrxx);

	void init_pot(
		const int &istep, // ionic steps
		ModuleBase::ComplexMatrix &sf // structure factors
	);

    ModuleBase::matrix v_of_rho(
		const double*const*const rho_in,
		const double * const rho_core_in);

    void set_vr_eff(void);

    void newd(void);

	public:

	// mohan add 2011-02-28
	// here vh is std::complex because the array is got after std::complex FFT.
	void write_potential(const int &is, const int &iter, const std::string &fn,
		const ModuleBase::matrix &v, const int &precision, const int &hartree = 0)const;

    void write_elecstat_pot(const std::string &fn, const std::string &fn_ave);
	
	private:


	void set_local_pot(
		double* vl_pseudo, // store the local pseudopotential
		const int &ntype, // number of atom types
		const int &ngmc, // number of |g|, g is plane wave
		ModuleBase::matrix &vloc, // local pseduopotentials
		int* ig2ngg, // ig2ngg
		ModuleBase::ComplexMatrix &sf // structure factors	
	)const;

	// TDDFT related, fuxiang add
    double *vext;

    double *vextold;

    void set_vrs_tddft(const int istep);

};

#endif	// POTENTIAL_H
