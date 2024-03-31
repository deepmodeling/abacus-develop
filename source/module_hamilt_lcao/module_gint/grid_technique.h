#ifndef GRID_TECHNIQUE_H
#define GRID_TECHNIQUE_H

#include "grid_index.h"
#include "grid_meshball.h"
#include "module_basis/module_ao/parallel_orbitals.h"
#if ((defined __CUDA) /* || (defined __ROCM) */)
#include <cuda_runtime.h>

#include "kernels/cuda/cuda_tools.cuh"
#include "kernels/cuda/vbatch_matrix_mul.cuh"
#endif

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
    int lnat;      // local nat.
    int lgd;       // local grid dimension.  lgd * lgd symmetry matrix.
    int* trace_lo; // trace local orbital.

    //---------------------------------------
    // nnrg: number of matrix elements on
    // each processor's real space grid.
    // use: GridT.in_this_processor
    //---------------------------------------
    int nnrg;
    int* nlocdimg;
    int* nlocstartg;

    int* nad; // number of adjacent atoms for each atom.
    int** find_R2;
    int** find_R2_sorted_index;
    int** find_R2st;
    bool allocate_find_R2;
    int binary_search_find_R2_offset(int val, int iat) const;

    // indexes for nnrg -> orbital index + R index
    std::vector<gridIntegral::gridIndex> nnrg_index;

    // public functions
  public:
    Grid_Technique();
    ~Grid_Technique();

    void set_pbc_grid(const int& ncx_in,
                      const int& ncy_in,
                      const int& ncz_in,
                      const int& bx_in,
                      const int& by_in,
                      const int& bz_in,
                      const int& nbx_in,
                      const int& nby_in,
                      const int& nbz_in,
                      const int& nbxx_in,
                      const int& nbzp_start_in,
                      const int& nbzp_in,
                      const int& ny,
                      const int& nplane,
                      const int& startz_current);

    /// number of elements(basis-pairs) in this processon
    /// on all adjacent atoms-pairs(Grid division)
    void cal_nnrg(Parallel_Orbitals* pv);
    int cal_RindexAtom(const int& u1,
                       const int& u2,
                       const int& u3,
                       const int& iat2) const;

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
    void init_atoms_on_grid(const int& ny,
                            const int& nplane,
                            const int& startz_current);
    void init_atoms_on_grid2(const int* index2normal);
    void cal_grid_integration_index(void);
    void cal_trace_lo(void);
    void check_bigcell(int*& ind_bigcell, bool*& bigcell_on_processor);
    void get_startind(const int& ny,
                      const int& nplane,
                      const int& startz_current);

#if ((defined __CUDA) /* || (defined __ROCM) */)
  public:
    double* ylmcoef_g;
    bool is_malloced;

    int* atom_nw_g;
    int* ucell_atom_nwl_g;
    double* psi_u_g;
    bool* atom_iw2_new_g;
    int* atom_iw2_ylm_g;
    int* atom_iw2_l_g;
    double** grid_vlocal_g;
    int nr_max;
    int psi_size_max;
    int psi_size_max_z;
    int psir_size;
    int atom_pair_mesh;
    int atom_pair_nbz;

    const int nstreams = 4;
    cudaStream_t streams[4];
    // streams[nstreams]
    // TODO it needs to be implemented through configuration files

    double* left_global_g;
    double* d_left_x_g;
    double* d_left_y_g;
    double* d_left_z_g;

    double* dd_left_xx_g;
    double* dd_left_xy_g;
    double* dd_left_xz_g;
    double* dd_left_yy_g;
    double* dd_left_yz_g;
    double* dd_left_zz_g;
    double* right_global_g;
    double* dm_global_g;

    double* alpha_global;
    double* alpha_global_g;
    int* l_info_global;
    int* l_info_global_g;
    int* r_info_global;
    int* r_info_global_g;
    int* k_info_global;
    int* k_info_global_g;

    int* lda_info_global;
    int* lda_info_global_g;
    int* ldb_info_global;
    int* ldb_info_global_g;
    int* ldc_info_global;
    int* ldc_info_global_g;

    double** ap_left_glo;
    double** ap_right_glo;
    double** ap_output_glo;

    double** ap_left_glo_g;
    double** ap_right_glo_g;
    double** ap_output_glo_g;

    double* psi_dou_glo;
    double* psi_dou_glo_g;

    int* psi_int_glo;
    int* psi_int_glo_g;

    int* num_psir_glo;
    int* num_psir_glo_g;

    // additional variables for rho calculating
    int num_mcell;
    double* rho_g;
    int* vec_len;
    int* vec_len_g;
    double** vec_l;
    double** vec_l_g;
    double** vec_r;
    double** vec_r_g;
    double** dot_product;
    double** dot_product_g;

    matrix_multiple_func_type fastest_matrix_mul;

  private:
    void init_gpu_gint_variables();
    void free_gpu_gint_variables();

#endif
};
#endif
