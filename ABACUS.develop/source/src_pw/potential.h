#ifndef POTENTIAL_H
#define POTENTIAL_H

#include "tools.h"

class Potential
{
	public:

	friend class Hamilt_PW; 
	friend class Electrons;
	friend class ELEC_scf;

    // constructor and deconstructor
    Potential();
    ~Potential();

    //==========================================================
    // start_pot : "atomic" or "file"
	// extra_pot : extrapolation methods for potential
    // vr(nspin,ncxyz) : Hartree + xc potentials in real space
    // vr_eff(nspin,ncxyz) : effective potential in real space 
    // vnew(nspin,ncxyz) : V_out - V_in, needed in scf
	// vltot: the local potential in real space
	// out_potential: options to print out potentials 
    //==========================================================

    string start_pot;
	string extra_pot;
    matrix vr;
    matrix vr_eff;
    matrix vnew;
    double *vr_eff1; 
    double *vltot;
	int out_potential; // mohan add 2011-02-28

    void allocate(const int nrxx);

	void init_pot(
		const int &istep, // ionic steps
		ComplexMatrix &sf // structure factors
	);

    matrix v_of_rho(
		const double*const*const rho_in,
		const double * const rho_core_in);

    void set_vr_eff(void);

    void newd(void);

	public:

	// mohan add 2011-02-28
	// here vh is complex because the array is got after complex FFT.
	void write_potential(const int &is, const int &iter, const string &fn,
		const matrix &v, const int &precision, const int &hartree = 0)const;

    void write_elecstat_pot(const string &fn, const string &fn_ave);
	
	private:


	void set_local_pot(
		double* vl_pseudo, // store the local pseudopotential
		const int &ntype, // number of atom types
		const int &ngmc, // number of |g|, g is plane wave
		matrix &vloc, // local pseduopotentials
		int* ig2ngg, // ig2ngg
		ComplexMatrix &sf // structure factors	
	)const;

	// TDDFT related, fuxiang add
    double *vext;

    double *vextold;

    void set_vrs_tddft(const int istep);

};

#endif	// POTENTIAL_H
