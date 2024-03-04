//=========================================================
// REFACTOR : Peize Lin, 2021.06.28
//=========================================================
#ifndef GINT_TOOLS_H
#define GINT_TOOLS_H
#include "grid_technique.h"
#include "module_elecstate/module_charge/charge.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"
#include <cstdlib>

namespace Gint_Tools {
enum class job_type {
  vlocal,
  rho,
  force,
  tau,
  vlocal_meta,
  force_meta,
  dvlocal
};
// Hamiltonian, electron density, force, kinetic energy density, Hamiltonian for
// mGGA
} // namespace Gint_Tools

// the class is used to pass input/output variables
// into the unified interface gint
// not sure if this is the best practice though ..
class Gint_inout {
public:
  // input
  double **DM_R;
  double ***DM;
  const double *vl;
  const double *vofk;
  bool isforce;
  bool isstress;
  int ispin;
  bool if_symm =
      false; // if true, use dsymv in gint_kernel_rho; if false, use dgemv.

  // output
  double **rho;
  ModuleBase::matrix *fvl_dphi;
  ModuleBase::matrix *svl_dphi;

  Gint_Tools::job_type job;

  // electron density and kin_r, multi-k
  Gint_inout(double **DM_R_in, double **rho_in, Gint_Tools::job_type job_in,
             bool if_symm_in = true) {
    DM_R = DM_R_in;
    rho = rho_in;
    job = job_in;
    if_symm = if_symm_in;
  }

  // force, multi-k
  Gint_inout(double **DM_R_in, const int ispin_in, const double *vl_in,
             bool isforce_in, bool isstress_in, ModuleBase::matrix *fvl_dphi_in,
             ModuleBase::matrix *svl_dphi_in, Gint_Tools::job_type job_in) {
    DM_R = DM_R_in;
    vl = vl_in;
    isforce = isforce_in;
    isstress = isstress_in;
    fvl_dphi = fvl_dphi_in;
    svl_dphi = svl_dphi_in;
    job = job_in;
    ispin = ispin_in;
  }

  // force (mGGA), multi-k
  Gint_inout(double **DM_R_in, const int ispin_in, const double *vl_in,
             const double *vofk_in, const bool isforce_in,
             const bool isstress_in, ModuleBase::matrix *fvl_dphi_in,
             ModuleBase::matrix *svl_dphi_in, Gint_Tools::job_type job_in) {
    DM_R = DM_R_in;
    vl = vl_in;
    vofk = vofk_in;
    isforce = isforce_in;
    isstress = isstress_in;
    fvl_dphi = fvl_dphi_in;
    svl_dphi = svl_dphi_in;
    job = job_in;
    ispin = ispin_in;
  }

  // vlocal, multi-k
  Gint_inout(const double *vl_in, int ispin_in, Gint_Tools::job_type job_in) {
    vl = vl_in;
    ispin = ispin_in;
    job = job_in;
  }

  // mGGA vlocal, multi-k
  Gint_inout(const double *vl_in, const double *vofk_in, int ispin_in,
             Gint_Tools::job_type job_in) {
    vl = vl_in;
    vofk = vofk_in;
    ispin = ispin_in;
    job = job_in;
  }

  // electron density and kin_r, gamma point
  Gint_inout(double ***DM_in, double **rho_in, Gint_Tools::job_type job_in,
             bool if_symm_in = true) {
    DM = DM_in;
    rho = rho_in;
    job = job_in;
    if_symm = if_symm_in;
  }

  // force, gamma point
  Gint_inout(double ***DM_in, const int ispin_in, const double *vl_in,
             const bool isforce_in, const bool isstress_in,
             ModuleBase::matrix *fvl_dphi_in, ModuleBase::matrix *svl_dphi_in,
             Gint_Tools::job_type job_in) {
    DM = DM_in;
    vl = vl_in;
    isforce = isforce_in;
    isstress = isstress_in;
    fvl_dphi = fvl_dphi_in;
    svl_dphi = svl_dphi_in;
    job = job_in;
    ispin = ispin_in;
  }

  // force (mGGA), gamma point
  Gint_inout(double ***DM_in, const int ispin_in, const double *vl_in,
             const double *vofk_in, const bool isforce_in,
             const bool isstress_in, ModuleBase::matrix *fvl_dphi_in,
             ModuleBase::matrix *svl_dphi_in, Gint_Tools::job_type job_in) {
    DM = DM_in;
    vl = vl_in;
    vofk = vofk_in;
    isforce = isforce_in;
    isstress = isstress_in;
    fvl_dphi = fvl_dphi_in;
    svl_dphi = svl_dphi_in;
    job = job_in;
    ispin = ispin_in;
  }

  // vlocal, gamma point
  Gint_inout(const double *vl_in, Gint_Tools::job_type job_in) {
    vl = vl_in;
    job = job_in;
  }

  // mGGA vlocal, gamma point
  Gint_inout(const double *vl_in, const double *vofk_in,
             Gint_Tools::job_type job_in) {
    vl = vl_in;
    vofk = vofk_in;
    job = job_in;
  }
};

namespace Gint_Tools {
template <typename T> class Array_Pool {
public:
  T **ptr_2D;
  T *ptr_1D;
  Array_Pool(const int nr, const int nc);
  Array_Pool(Array_Pool<T> &&array);
  ~Array_Pool();
  Array_Pool(const Array_Pool<T> &array) = delete;
  Array_Pool(Array_Pool<T> &array) = delete;
};

// vindex[pw.bxyz]
int *get_vindex(const int bxyz, const int bx, const int by, const int bz,
                const int nplane, const int ncyz, const int ibx, const int jby,
                const int kbz);

int *get_vindex(const int bxyz, const int bx, const int by, const int bz,
                const int nplane, const int start_ind, const int ncyz);

// extract the local potentials.
// vldr3[bxyz]
double *get_vldr3(const double *const vlocal, const int bxyz, const int bx,
                  const int by, const int bz, const int nplane, const int ncyz,
                  const int ibx, const int jby, const int kbz, const double dv);

double *get_vldr3(const double *const vlocal, const int bxyz, const int bx,
                  const int by, const int bz, const int nplane,
                  const int start_ind, const int ncyz, const double dv);

//------------------------------------------------------
// na_grid : #. atoms for this group of grids
// block_iw : size na_grid, index of the first orbital on this atom
// block_size : size na_grid, number of orbitals on this atom
// block_index : size na_grid+1, start from 0, accumulates block_size
// cal_flag : whether the atom-grid distance is larger than cutoff
//------------------------------------------------------
void get_block_info(const Grid_Technique &gt, const int bxyz, const int na_grid,
                    const int grid_index, int *&block_iw, int *&block_index,
                    int *&block_size, bool **&cal_flag);

// psir_ylm[pw.bxyz][LD_pool]
void cal_psir_ylm(
    const Grid_Technique &gt, const int bxyz,
    const int na_grid,            // number of atoms on this grid
    const int grid_index,         // 1d index of FFT index (i,j,k)
    const double delta_r,         // delta_r of the uniform FFT grid
    const int *const block_index, // count total number of atomis orbitals
    const int *const block_size, const bool *const *const cal_flag,
    double *const *const
        psir_ylm); // whether the atom-grid distance is larger than cutoff

// psir_ylm and dpsir_ylm, both[pw.bxyz][LD_pool]
void cal_dpsir_ylm(
    const Grid_Technique &gt, const int bxyz,
    const int na_grid,            // number of atoms on this grid
    const int grid_index,         // 1d index of FFT index (i,j,k)
    const double delta_r,         // delta_r of the uniform FFT grid
    const int *const block_index, // block_index[na_grid+1], count total number
                                  // of atomis orbitals
    const int *const
        block_size, // block_size[na_grid],	number of columns of a band
    const bool *const
        *const cal_flag, // cal_flag[bxyz][na_grid],	whether the atom-grid
                         // distance is larger than cutoff
    double *const *const psir_ylm, double *const *const dpsir_ylm_x,
    double *const *const dpsir_ylm_y, double *const *const dpsir_ylm_z);

// dpsir_ylm * (r-R), R is the atomic position
void cal_dpsirr_ylm(
    const Grid_Technique &gt, const int bxyz,
    const int na_grid,            // number of atoms on this grid
    const int grid_index,         // 1d index of FFT index (i,j,k)
    const int *const block_index, // block_index[na_grid+1], count total number
                                  // of atomis orbitals
    const int *const
        block_size, // block_size[na_grid],	number of columns of a band
    const bool *const
        *const cal_flag, // cal_flag[bxyz][na_grid],	whether the atom-grid
                         // distance is larger than cutoff
    double *const *const dpsir_ylm_x, double *const *const dpsir_ylm_y,
    double *const *const dpsir_ylm_z, double *const *const dpsir_ylm_xx,
    double *const *const dpsir_ylm_xy, double *const *const dpsir_ylm_xz,
    double *const *const dpsir_ylm_yy, double *const *const dpsir_ylm_yz,
    double *const *const dpsir_ylm_zz);

void cal_ddpsir_ylm(
    const Grid_Technique &gt, const int bxyz,
    const int na_grid,            // number of atoms on this grid
    const int grid_index,         // 1d index of FFT index (i,j,k)
    const double delta_r,         // delta_r of the uniform FFT grid
    const int *const block_index, // block_index[na_grid+1], count total number
                                  // of atomis orbitals
    const int *const
        block_size, // block_size[na_grid],	number of columns of a band
    const bool *const
        *const cal_flag, // cal_flag[bxyz][na_grid],	whether the atom-grid
                         // distance is larger than cutoff
    double *const *const ddpsir_ylm_xx, double *const *const ddpsir_ylm_xy,
    double *const *const ddpsir_ylm_xz, double *const *const ddpsir_ylm_yy,
    double *const *const ddpsir_ylm_yz, double *const *const ddpsir_ylm_zz);

// psir_ylm * vldr3
Gint_Tools::Array_Pool<double>
get_psir_vlbr3(const int bxyz,
               const int na_grid, // how many atoms on this (i,j,k) grid
               const int LD_pool,
               const int *const block_index, // block_index[na_grid+1], count
                                             // total number of atomis orbitals
               const bool *const
                   *const cal_flag, // cal_flag[bxyz][na_grid],	whether the
                                    // atom-grid distance is larger than cutoff
               const double *const vldr3,            // vldr3[bxyz]
               const double *const *const psir_ylm); // psir_ylm[bxyz][LD_pool]

// sum_nu rho_mu,nu psi_nu, for gamma point
void mult_psi_DM(
    const Grid_Technique &gt, const int bxyz,
    const int na_grid, // how many atoms on this (i,j,k) grid
    const int LD_pool,
    const int *const
        block_iw, // block_iw[na_grid],	index of wave functions for each block
    const int *const
        block_size, // block_size[na_grid],	number of columns of a band
    const int *const block_index, // block_index[na_grid+1], count total number
                                  // of atomis orbitals
    const bool *const
        *const cal_flag, // cal_flag[bxyz][na_grid],	whether the atom-grid
                         // distance is larger than cutoff
    const double *const *const psi, // psir_vlbr3[bxyz][LD_pool]
    double **psi_DM, const double *const *const DM, const bool if_symm);

// sum_nu,R rho_mu,nu(R) psi_nu, for multi-k
void mult_psi_DMR(const Grid_Technique &gt, const int bxyz,
                  const int &grid_index, const int &na_grid,
                  const int *const block_index, const int *const block_size,
                  bool **cal_flag, double **psi, double **psi_DMR, double *DMR,
                  const hamilt::HContainer<double> *DM, const bool if_symm);

// sum_nu rho_mu,nu psi_nu, for gamma point
void mult_psi_DM_new(
    const Grid_Technique &gt, const int bxyz, const int &grid_index,
    const int na_grid, // how many atoms on this (i,j,k) grid
    const int LD_pool,
    const int *const
        block_iw, // block_iw[na_grid],	index of wave functions for each block
    const int *const
        block_size, // block_size[na_grid],	number of columns of a band
    const int *const block_index, // block_index[na_grid+1], count total number
                                  // of atomis orbitals
    const bool *const
        *const cal_flag, // cal_flag[bxyz][na_grid],	whether the atom-grid
                         // distance is larger than cutoff
    const double *const *const psi, // psir_vlbr3[bxyz][LD_pool]
    double **psi_DM, const hamilt::HContainer<double> *DM, const bool if_symm);

} // namespace Gint_Tools

namespace Gint_Tools {
template <typename T>
Array_Pool<T>::Array_Pool(const int nr,
                          const int nc) // Attention: uninitialized
{
  ptr_1D = new T[nr * nc];
  ptr_2D = new T *[nr];
  for (int ir = 0; ir < nr; ++ir)
    ptr_2D[ir] = &ptr_1D[ir * nc];
}

template <typename T> Array_Pool<T>::Array_Pool(Array_Pool<T> &&array) {
  ptr_1D = array.ptr_1D;
  ptr_2D = array.ptr_2D;
  delete[] array.ptr_2D;
  delete[] array.ptr_1D;
}

template <typename T> Array_Pool<T>::~Array_Pool() {
  delete[] ptr_2D;
  delete[] ptr_1D;
}
} // namespace Gint_Tools

#endif