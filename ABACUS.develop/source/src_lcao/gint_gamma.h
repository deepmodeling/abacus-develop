//=========================================================
//AUTHOR : mohan
//DATE : 2009-09-16
//=========================================================
#ifndef GINT_GAMMA_H
#define GINT_GAMMA_H

#include "../src_pw/tools.h"
#include "grid_base_beta.h"
#include "grid_technique.h"
#include "LCAO_matrix.h"
#include <omp.h>

//=========================================================
// Integral On 3D Grids, different from Grid_Integral
// Feature : Matrix Elements Of Local Potential For 
// Numerical Orbitals
//=========================================================

class Gint_Gamma : public Grid_Base_Beta
{
	public:

	Gint_Gamma();
	~Gint_Gamma();

	// (1) calculate the H matrix in terms of effective potentials
	void cal_vlocal( const double* vlocal_in);

	// (2) calculate charge density
	double cal_rho(void);

	// (3) calcualte the forces related to grid
	void cal_force( const double* vlocal_in);

	// (4) calcualte the envelope function
	void cal_env(const double* wfc, double* rho);

	// (5) calculate the Mulliken charge
	void cal_mulliken(double** mulliken);


	private:	

	double* transformer;
	double psiv1;
	double psiv2;
	double* ylm1;
	double* ylm2;

	int grid_index;			// may delete?
	int max_size;
	
	// these parameters are for interpolation.
	// we store these parameters at first to speed
	// up the calculation.
	double *x0;
	double *x1;
	double *x2;
	double *x3;
	double* x12;
	double* x03;
	int *iq;

	// peize add, I guess, mohan add 2021-01-31
	omp_lock_t lock;

	void save_atoms_on_grid(const Grid_Technique &gt);

	// for calculation of < phi_i | Vlocal | phi_j >
	void gamma_vlocal(void);  

	// for calculation of charege 
	double gamma_charge(void);

	// for calculation of Mulliken charge.
	void gamma_mulliken(double** mulliken);

	// for calculation of envelope functions.
	void gamma_envelope(const double* wfc, double* rho);// mohan add 2011-07-01


	// for calculatin of < dphi_i | Vlocal | phi_j > for foce calculation
	// on regular FFT real space grid.
	void gamma_force(void);

	void cal_meshball_vlocal(
		const int na_grid,
		const int LD_pool,
		const int*const block_iw,
		const int*const block_size,
		const int*const block_index,
		const bool*const*const cal_flag,
		const double*const vldr3,
		const double*const*const psir_ylm,
		const double*const*const psir_vlbr3,
		const int lgd_now,
		double*const*const GridVlocal);

	void cal_band_rho(
		const int na_grid, 
		const int LD_pool, 
		const int*const block_iw, 
		const int*const block_size, 
		const int*const block_index,
		const bool*const*const cal_flag, 
		const double*const*const psir_ylm,
		const int*const vindex);
	
	// extract the local potentials.
	// vldr3[pw.bxyz]
	double* get_vldr3( const int ncyz, const int ibx, const int jby, const int kbz) const;
};

#endif
