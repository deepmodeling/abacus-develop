//=========================================================
// REFACTOR : Peize Lin, 2021.06.28
//=========================================================
#ifndef GINT_TOOLS_H
#define GINT_TOOLS_H
#include "grid_technique.h"
#include "module_elecstate/module_charge/charge.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"
#include "module_base/array_pool.h"

#include <cstdlib>
#include <utility> // for std::pair

namespace Gint_Tools
{
enum class job_type
{
    vlocal,
    rho,
    force,
    tau,
    vlocal_meta,
    force_meta,
    dvlocal
};
// Hamiltonian, electron density, force, kinetic energy density, Hamiltonian for mGGA
} // namespace Gint_Tools

// the class is used to pass input/output variables
// into the unified interface gint
// not sure if this is the best practice though ..
class Gint_inout
{
  public:
    // input
    double*** DM;
    const double* vl;
    const double* vofk;
    bool isforce;
    bool isstress;
    int ispin;
    int nspin_rho;  // usually, but not always, equal to global nspin
    bool if_symm = false; // if true, use dsymv in gint_kernel_rho; if false, use dgemv.

    // output
    double** rho;
    ModuleBase::matrix* fvl_dphi;
    ModuleBase::matrix* svl_dphi;
    Gint_Tools::job_type job;

    // electron density and kin_r, multi-k
    Gint_inout(double** rho_in, Gint_Tools::job_type job_in, const int& nspin_rho_in, bool if_symm_in = true)
    {
        rho = rho_in;
        job = job_in;
        nspin_rho = nspin_rho_in;
        if_symm = if_symm_in;
    }

    // force
    Gint_inout(const int ispin_in,
               const double* vl_in,
               bool isforce_in,
               bool isstress_in,
               ModuleBase::matrix* fvl_dphi_in,
               ModuleBase::matrix* svl_dphi_in,
               Gint_Tools::job_type job_in)
    {
        vl = vl_in;
        isforce = isforce_in;
        isstress = isstress_in;
        fvl_dphi = fvl_dphi_in;
        svl_dphi = svl_dphi_in;
        job = job_in;
        ispin = ispin_in;
    }

    // force (mGGA)
    Gint_inout(const int ispin_in,
               const double* vl_in,
               const double* vofk_in,
               const bool isforce_in,
               const bool isstress_in,
               ModuleBase::matrix* fvl_dphi_in,
               ModuleBase::matrix* svl_dphi_in,
               Gint_Tools::job_type job_in)
    {
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
    Gint_inout(const double* vl_in, int ispin_in, Gint_Tools::job_type job_in)
    {
        vl = vl_in;
        ispin = ispin_in;
        job = job_in;
    }

    // mGGA vlocal, multi-k
    Gint_inout(const double* vl_in, const double* vofk_in, int ispin_in, Gint_Tools::job_type job_in)
    {
        vl = vl_in;
        vofk = vofk_in;
        ispin = ispin_in;
        job = job_in;
    }

    // vlocal, gamma point
    Gint_inout(const double* vl_in, Gint_Tools::job_type job_in)
    {
        vl = vl_in;
        job = job_in;
    }

    // mGGA vlocal, gamma point
    Gint_inout(const double* vl_in, const double* vofk_in, Gint_Tools::job_type job_in)
    {
        vl = vl_in;
        vofk = vofk_in;
        job = job_in;
    }
};

namespace Gint_Tools
{
// if exponent is an integer between 0 and 5 (the most common cases in gint),
// pow_int is much faster than std::pow
inline double pow_int(const double base, const int exp)
{
    switch (exp)
    {
    case 0:
        return 1.0;
    case 1:
        return base;
    case 2:
        return base * base;
    case 3:
        return base * base * base;
    case 4:
        return base * base * base * base;
    case 5:
        return base * base * base * base * base;
    default:
        double result = std::pow(base, exp);
        return result;
    }
}
// vindex[pw.bxyz]

/**
 * @brief Get the vindex form the grid index
 * @param bxyz number of big grids
 * @param bx number of big grids in x direction
 * @param by number of big grids in y direction
 * @param bz number of big grids in z direction
 * @param nplane Currently using Z-axis 1D division, 
 * recording the number of the Z-axis process
 * (nbz in the current process).
 * @param start_ind start index of the grid in the 1D FFT grid
 * @param ncyz number of grids in yz plane
 * @param vindex the index of the grid 
*/
void get_vindex(const int bxyz, const int bx, const int by,
                    const int bz, const int nplane, 
                    const int start_ind,const int ncyz,int* vindex);

/**
 * @brief Get the vldr3 form the grid index
 * @param vldr3 the local potential multiplied by the grid volume
 * @param vlocal the local potential
 * @param bxyz number of grids
 * @param bx number of grids in x direction
 * @param by number of grids in y direction
 * @param bz number of grids in z direction
 * @param nplane Currently using Z-axis 1D division, 
 * recording the number of the Z-axis process
 * (nbz in the current process).
 * @param start_ind start index of the grid in the 1D FFT grid
 * @param ncyz number of grids in yz plane
 * @param dv the volume of the grid
*/
void get_gint_vldr3(double* vldr3,
                    const double* const vlocal,
                    const int bxyz,
                    const int bx,
                    const int by,
                    const int bz,
                    const int nplane,
                    const int start_ind,
                    const int ncyz,
                    const double dv);

/**
 * @brief Get the information of a big grid index
 * @param gt the grid technique, which contains the tools of the grid intergration
 * @param bxyz number of grids
 * @param na_grid number of atoms on this grid
 * @param grid_index 1d index of FFT index (i,j,k)
 * @param block_iw track the atom orbitals in all atoms
 * @param block_index count total number of atomis orbitals
 * @param block_size count the number of atomis orbitals in each atom
 * @param cal_flag whether the atom-grid distance is larger than cutoff
*/                    
void get_block_info(const Grid_Technique& gt, const int bxyz, const int na_grid, const int grid_index,
                    int* block_iw, int* block_index, int* block_size, bool** cal_flag);

void init_orb(double& dr_uniform,
              std::vector<double>& rcuts,
              UnitCell& ucell,
              const LCAO_Orbitals& orb,
              std::vector<std::vector<double>>& psi_u,
              std::vector<std::vector<double>>& dpsi_u,
              std::vector<std::vector<double>>& d2psi_u);

// psir_ylm[pw.bxyz][LD_pool]
void cal_psir_ylm(const Grid_Technique& gt,
                  const int bxyz,
                  const int na_grid,            // number of atoms on this grid
                  const int grid_index,         // 1d index of FFT index (i,j,k)
                  const double delta_r,         // delta_r of the uniform FFT grid
                  const int* const block_index, // count total number of atomis orbitals
                  const int* const block_size,
                  const bool* const* const cal_flag,
                  double* const* const psir_ylm); // whether the atom-grid distance is larger than cutoff

// psir_ylm and dpsir_ylm, both[pw.bxyz][LD_pool]
void cal_dpsir_ylm(
    const Grid_Technique& gt,
    const int bxyz,
    const int na_grid,                 // number of atoms on this grid
    const int grid_index,              // 1d index of FFT index (i,j,k)
    const double delta_r,              // delta_r of the uniform FFT grid
    const int* const block_index,      // block_index[na_grid+1], count total number of atomis orbitals
    const int* const block_size,       // block_size[na_grid],	number of columns of a band
    const bool* const* const cal_flag, // cal_flag[bxyz][na_grid],	whether the atom-grid distance is larger than cutoff
    double* const* const psir_ylm,
    double* const* const dpsir_ylm_x,
    double* const* const dpsir_ylm_y,
    double* const* const dpsir_ylm_z);

// dpsir_ylm * (r-R), R is the atomic position
void cal_dpsirr_ylm(
    const Grid_Technique& gt, const int bxyz,
    const int na_grid,                 // number of atoms on this grid
    const int grid_index,              // 1d index of FFT index (i,j,k)
    const int* const block_index,      // block_index[na_grid+1], count total number of atomis orbitals
    const int* const block_size,       // block_size[na_grid],	number of columns of a band
    const bool* const* const cal_flag, // cal_flag[bxyz][na_grid],	whether the atom-grid distance is larger than cutoff
    double* const* const dpsir_ylm_x, double* const* const dpsir_ylm_y, double* const* const dpsir_ylm_z,
    double* const* const dpsir_ylm);

void cal_ddpsir_ylm(
    const Grid_Technique& gt,
    const int bxyz,
    const int na_grid,                 // number of atoms on this grid
    const int grid_index,              // 1d index of FFT index (i,j,k)
    const double delta_r,              // delta_r of the uniform FFT grid
    const int* const block_index,      // block_index[na_grid+1], count total number of atomis orbitals
    const int* const block_size,       // block_size[na_grid],	number of columns of a band
    const bool* const* const cal_flag, // cal_flag[bxyz][na_grid],	whether the atom-grid distance is larger than cutoff
    double* const* const ddpsir_ylm_xx,
    double* const* const ddpsir_ylm_xy,
    double* const* const ddpsir_ylm_xz,
    double* const* const ddpsir_ylm_yy,
    double* const* const ddpsir_ylm_yz,
    double* const* const ddpsir_ylm_zz);

// psir_ylm * vldr3
ModuleBase::Array_Pool<double> get_psir_vlbr3(
    const int bxyz,
    const int na_grid, // how many atoms on this (i,j,k) grid
    const int LD_pool,
    const int* const block_index,      // block_index[na_grid+1], count total number of atomis orbitals
    const bool* const* const cal_flag, // cal_flag[bxyz][na_grid],	whether the atom-grid distance is larger than cutoff
    const double* const vldr3,         // vldr3[bxyz]
    const double* const* const psir_ylm); // psir_ylm[bxyz][LD_pool]

// sum_nu,R rho_mu,nu(R) psi_nu, for multi-k and gamma point
void mult_psi_DMR(
    const Grid_Technique& gt,
    const int bxyz,
    const int LD_pool,
    const int &grid_index,
    const int &na_grid,
    const int*const block_index,
    const int*const block_size,
    const bool*const*const cal_flag,
    const double*const*const psi,
    double*const*const psi_DMR,
    const hamilt::HContainer<double>*const DM,
    const bool if_symm);


// pair.first is the first index of the meshcell which is inside atoms ia1 and ia2.
// pair.second is the number of meshcells which should be calculated in the following gemm.
// If no meshcell is inside both ia1 and ia2, return [bxyz, 0].
std::pair<int, int> cal_info(const int bxyz, 
			                 const int ia1,
			                 const int ia2,
			                 const bool* const* const cal_flag);
            
} // namespace Gint_Tools
#endif
