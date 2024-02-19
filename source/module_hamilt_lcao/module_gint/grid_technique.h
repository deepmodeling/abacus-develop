#ifndef GRID_TECHNIQUE_H
#define GRID_TECHNIQUE_H

#include "grid_meshball.h"
#include "grid_index.h"
#include "module_basis/module_ao/parallel_orbitals.h"
#if ((defined __CUDA) /* || (defined __ROCM) */)
#include <cuda_runtime.h>
#endif
#include "gint_gpu/vbatch_matrix_multiple/cuda_tools.cuh"
#include "gint_gpu/vbatch_matrix_multiple/vbatch_matrix_mul.cuh"

// Author: mohan
// Date: 2009-10-17
class Grid_Technique : public Grid_MeshBall
{
	// public variables.
	public:
	
	//------------------------------------
	// 1: Info about atom number on grid.
	//------------------------------------
	// record how many atoms on each grid.
	int* how_many_atoms;
	// max atom on grid
	int max_atom;
	// sum of how_many_atoms
	int total_atoms_on_grid;

	int* start_ind;

	//------------------------------------
	// 2: Info about which atom on grid.
	//------------------------------------
	// save the start position of each big cell's adjacent
	// atoms in 1D grid.
	int* bcell_start;
	// save the 'iat' atom. 
	// dim: total_atoms_on_grid.
	int* which_atom;

	//--------------------------------------
	// save the bigcell index in meshball.
	// dim: total_atoms_on_grid.
	//--------------------------------------
	int* which_bigcell;
	int* which_unitcell;
	
	//------------------------------------
	// 3: which atom on local grid.
	//------------------------------------
	bool* in_this_processor;
	std::vector<int> trace_iat;
	int lnat; // local nat.
	int lgd; // local grid dimension.  lgd * lgd symmetry matrix. 
	int* trace_lo; // trace local orbital.

    //---------------------------------------
	// nnrg: number of matrix elements on
	// each processor's real space grid.
	// use: GridT.in_this_processor
	//---------------------------------------
	int nnrg;
	int *nlocdimg;
	int *nlocstartg;
    
    int* nad; // number of adjacent atoms for each atom.
	int **find_R2;
	int **find_R2_sorted_index;
	int **find_R2st;
    bool allocate_find_R2;
	int binary_search_find_R2_offset(int val, int iat) const;

	//indexes for nnrg -> orbital index + R index
	std::vector<gridIntegral::gridIndex> nnrg_index;
    
    // public functions
	public:

	Grid_Technique();
	~Grid_Technique();

	void set_pbc_grid(
			const int &ncx_in,
			const int &ncy_in,
			const int &ncz_in,
			const int &bx_in,
			const int &by_in,
			const int &bz_in,
			const int &nbx_in,
			const int &nby_in,
			const int &nbz_in,
			const int &nbxx_in,
			const int &nbzp_start_in,
            const int& nbzp_in,
            const int& ny,
            const int& nplane,
            const int& startz_current);

    /// number of elements(basis-pairs) in this processon
    /// on all adjacent atoms-pairs(Grid division)
    void cal_nnrg(Parallel_Orbitals* pv);
    int cal_RindexAtom(const int& u1, const int& u2, const int& u3, const int& iat2) const;
    
private:

    void cal_max_box_index(void);
    
    int maxB1;
    int maxB2;
	int maxB3;

	int minB1;
	int minB2;
	int minB3;

	int nB1;
	int nB2;
	int nB3;

    int nbox;

	// atoms on meshball
	void init_atoms_on_grid(const int& ny, const int& nplane, const int& startz_current);
	void init_atoms_on_grid2(const int* index2normal);
	void cal_grid_integration_index(void);
	void cal_trace_lo(void);
	void check_bigcell(int* &ind_bigcell, bool* &bigcell_on_processor);
	void get_startind(const int& ny, const int& nplane, const int& startz_current);

#if ((defined __CUDA) /* || (defined __ROCM) */)
public:
	double* ylmcoef_g;

	int *atom_nw_g;
	int *ucell_atom_nwl_g;
	double *psi_u_g;
	bool *atom_iw2_new_g;
	int *atom_iw2_ylm_g;
	int *atom_iw2_l_g;
	double **GridVlocal_v2_g;
	int nr_max;
	int psi_size_max;
	int psi_size_max_per_z;
	int psir_size;
	int atom_pair_size_of_meshcell;
	int atom_pair_size_over_nbz;

	const int nstreams = 4;
	cudaStream_t streams[4];
	// streams[nstreams]
	// TODO it needs to be implemented through configuration files

	double *psir_ylm_left_global_g;
	double *dpsir_ylm_left_x_global_g;
	double *dpsir_ylm_left_y_global_g;
	double *dpsir_ylm_left_z_global_g;

	double *ddpsir_ylm_left_xx_global_g;
	double *ddpsir_ylm_left_xy_global_g;
	double *ddpsir_ylm_left_xz_global_g;
	double *ddpsir_ylm_left_yy_global_g;
	double *ddpsir_ylm_left_yz_global_g;
	double *ddpsir_ylm_left_zz_global_g;
	double *psir_ylm_right_global_g;
	double *psir_ylm_dm_global_g;

	double *atom_pair_alpha_global;
	double *atom_pair_alpha_global_g;
	int *atom_pair_left_info_global;
	int *atom_pair_left_info_global_g;
	int *atom_pair_right_info_global;
	int *atom_pair_right_info_global_g;
	int *atom_pair_k_info_global;
	int *atom_pair_k_info_global_g;

	int *atom_pair_lda_global;
	int *atom_pair_lda_global_g;
	int *atom_pair_ldb_global;
	int *atom_pair_ldb_global_g;
	int *atom_pair_ldc_global;
	int *atom_pair_ldc_global_g;

	double **atom_pair_left_global;
	double **atom_pair_right_global;
	double **atom_pair_output_global;

	double **atom_pair_left_global_g;
	double **atom_pair_right_global_g;
	double **atom_pair_output_global_g;

	double *psi_input_double_global;
	double *psi_input_double_global_g;

	int *psi_input_int_global;
	int *psi_input_int_global_g;

	int *num_psir_global;
	int *num_psir_global_g;

	// additional variables for rho calculating
	int num_mcell;
	double *rho_g;
	int *vec_len;
	int *vec_len_g;
	double **vec_l;
	double **vec_l_g;
	double **vec_r;
	double **vec_r_g;
	double **dot_product;
	double **dot_product_g;



	matrix_multiple_func_type fastest_matrix_mul;
private:
	void init_gpu_gint_variables();
	
#endif

};
#endif
