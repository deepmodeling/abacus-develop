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

	template<typename T>
	class Array_Pool
	{
	public:
		Array_Pool(const int nr, const int nc);
		Array_Pool(Array_Pool<T> &&array);
		~Array_Pool();
		T** ptr_2D;
		T* ptr_1D;
		Array_Pool(const Array_Pool<T> &array) = delete;
		Array_Pool(Array_Pool<T> &array) = delete;
	};			

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

	// vindex[pw.bxyz]
	static int* get_vindex(
		const int ncyz,
		const int ibx,
		const int jby,
		const int kbz);
	
	// index of wave functions for each block
	// block_iw[na_grid]
	static int* get_block_iw(
		const int na_grid,  		// how many atoms on this (i,j,k) grid
		const int grid_index,		// 1d index of FFT index (i,j,k))
		const int max_size);
		
	// block_index[na_grid+1]
	static int* get_block_index(
		const int na_grid,  		// how many atoms on this (i,j,k) grid
		const int grid_index);		// 1d index of FFT index (i,j,k)
		
	// band size: number of columns of a band
	// block_size[na_grid]
	static int* get_block_size(
		const int na_grid,			// how many atoms on this (i,j,k) grid
		const int grid_index);		// 1d index of FFT index (i,j,k)

	// whether the atom-grid distance is larger than cutoff
	// cal_flag[pw.bxyz][na_grid]
	static bool** get_cal_flag(
		const int na_grid, 		// number of atoms on this grid 
		const int grid_index);		

	// psir_ylm[pw.bxyz][LD_pool]
	static Array_Pool<double> cal_psir_ylm(
		const int na_grid, // number of atoms on this grid 
		const int LD_pool,
		const int grid_index, // 1d index of FFT index (i,j,k) 
		const double delta_r, // delta_r of the uniform FFT grid
		const int*const block_index,  // count total number of atomis orbitals
		const int*const block_size, 
		const bool*const*const cal_flag); // whether the atom-grid distance is larger than cutoff

	// psir_vlbr3[pw.bxyz][LD_pool]
	static Array_Pool<double> get_psir_vlbr3(
		const int na_grid,
		const int LD_pool,
		const int*const block_index,
		const bool*const*const cal_flag,
		const double*const vldr3,
		const double*const*const psir_ylm);	
};


template<typename T>
Gint_Gamma::Array_Pool<T>::Array_Pool(const int nr, const int nc)	// Attention: uninitialized
{
	ptr_1D = (T*)malloc(nr*nc*sizeof(T));
	ptr_2D = (T**)malloc(nr*sizeof(T*));
	for (int ir=0; ir<nr; ++ir)
		ptr_2D[ir] = &ptr_1D[ir*nc];
}

template<typename T>
Gint_Gamma::Array_Pool<T>::Array_Pool(Array_Pool<T> &&array)
{
	ptr_1D = array.ptr_1D;
	ptr_2D = array.ptr_2D;
	free(array.ptr_2D);		array.ptr_2D=nullptr;
	free(array.ptr_1D);		array.ptr_1D=nullptr;
}

template<typename T>
Gint_Gamma::Array_Pool<T>::~Array_Pool()
{
	free(ptr_2D);
	free(ptr_1D);
}

#endif
