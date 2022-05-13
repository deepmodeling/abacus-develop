//=========================================================
//AUTHOR : mohan
//DATE : 2009-09-16
//REFACTOR : Peize Lin, 2021.06.28
//=========================================================
#ifndef GINT_GAMMA_H
#define GINT_GAMMA_H

#include "gint_tools.h"
#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "grid_technique.h"
#include "LCAO_matrix.h"
#include "../src_pw/charge.h"
#include <omp.h>

//=========================================================
// ModuleBase::Integral On 3D Grids, different from Grid_Integral
// Feature : Matrix Elements Of Local Potential For 
// Numerical Orbitals
//=========================================================

class Gint_Gamma
{
	public:

	Gint_Gamma();
	~Gint_Gamma();

	// (1) calculate the H matrix in terms of effective potentials
	void cal_vlocal( const double*const vlocal, LCAO_Matrix &lm);

	// (2) calculate charge density
	void cal_rho(double*** DM_in, Charge* chr);

	// (3) calcualte the forces related to grid
	void cal_force(double*** DM_in, const double*const vlocal,
			ModuleBase::matrix& force, ModuleBase::matrix& stress,
			const bool is_force, const bool is_stress);

	// (4) calcualte the envelope function
	void cal_env(const double* wfc, double* rho);

	// (5) calculate the Mulliken charge
	void cal_mulliken(double** mulliken);

private:

    double***  DM;   //pointer to LOC.DM

	int grid_index;			// may delete?
	int max_size;

    ///===============================
    /// Use MPI_Alltoallv to convert a grid distributed matrix
    /// to 2D - block cyclic distributed matrix.
    ///===============================
    int sender_index_size;
    int *sender_local_index;
    int sender_size;
    int *sender_size_process;
    int *sender_displacement_process;
    double* sender_buffer;

    int receiver_index_size;
    int *receiver_global_index;
    int receiver_size;
    int *receiver_size_process;
    int *receiver_displacement_process;
    double* receiver_buffer;

	// for calculation of < phi_i | Vlocal | phi_j >
	// Input:	vlocal[ir]
	// Output:	GridVlocal.ptr_2D[iw1_lo][iw2_lo]
	Gint_Tools::Array_Pool<double> gamma_vlocal(const double*const vlocal) const;  
	// for calculation of Mulliken charge.
	void gamma_mulliken(double** mulliken);

	// for calculation of envelope functions.
	void gamma_envelope(const double* wfc, double* rho);// mohan add 2011-07-01

	void cal_meshball_vlocal(
		const int na_grid,  						// how many atoms on this (i,j,k) grid
		const int LD_pool,
		const int*const block_iw,					// block_iw[na_grid],	index of wave functions for each block
		const int*const block_size, 				// block_size[na_grid],	number of columns of a band
		const int*const block_index,				// block_index[na_grid+1], count total number of atomis orbitals
		const bool*const*const cal_flag,			// cal_flag[GlobalC::pw.bxyz][na_grid],	whether the atom-grid distance is larger than cutoff
		const double*const vldr3,					// vldr3[GlobalC::pw.bxyz]
		const double*const*const psir_ylm,			// psir_ylm[GlobalC::pw.bxyz][LD_pool]
		const double*const*const psir_vlbr3,		// psir_vlbr3[GlobalC::pw.bxyz][LD_pool]
		const int lgd_now,
		double*const*const GridVlocal) const;		// GridVlocal[lgd_now][lgd_now]

	void cal_band_rho(
		const int na_grid,   							// how many atoms on this (i,j,k) grid
		const int LD_pool, 
		const int*const block_iw, 						// block_iw[na_grid],	index of wave functions for each block
		const int*const block_size, 					// block_size[na_grid],	band size: number of columns of a band
		const int*const block_index,					// block_index[na_grid+1], count total number of atomis orbitals
		const bool*const*const cal_flag, 				// cal_flag[GlobalC::pw.bxyz][na_grid],	whether the atom-grid distance is larger than cutoff
		const double*const*const psir_ylm,				// psir_ylm[GlobalC::pw.bxyz][LD_pool]
		const int*const vindex,							// vindex[GlobalC::pw.bxyz]
		const double*const*const*const DM,				// DM[GlobalV::NSPIN][lgd_now][lgd_now]
		Charge* chr) const;		// rho[GlobalV::NSPIN][GlobalC::pw.nrxx]

	void cal_meshball_force(
		const int grid_index,
		const int na_grid,  					    // how many atoms on this (i,j,k) grid
		const int*const block_size, 			    // block_size[na_grid],	number of columns of a band
		const int*const block_index,		    	// block_index[na_grid+1], count total number of atomis orbitals
		const double*const*const psir_vlbr3_DM,	    // psir_vlbr3[GlobalC::pw.bxyz][LD_pool]
		const double*const*const dpsir_x,	    // psir_vlbr3[GlobalC::pw.bxyz][LD_pool]
		const double*const*const dpsir_y,	    // psir_vlbr3[GlobalC::pw.bxyz][LD_pool]
		const double*const*const dpsir_z,	    // psir_vlbr3[GlobalC::pw.bxyz][LD_pool]
		ModuleBase::matrix &force
	);

	void cal_meshball_stress(
		const int na_grid,  					    // how many atoms on this (i,j,k) grid
		const int*const block_index,		    	// block_index[na_grid+1], count total number of atomis orbitals
		const double*const*const psir_vlbr3_DM,
		const double*const*const dpsir_xx,
		const double*const*const dpsir_xy,
		const double*const*const dpsir_xz,
		const double*const*const dpsir_yy,
		const double*const*const dpsir_yz,
		const double*const*const dpsir_zz,
		ModuleBase::matrix &stress
	);

    void vl_grid_to_2D(const Gint_Tools::Array_Pool<double>& GridVlocal, LCAO_Matrix& lm);
};

#endif
