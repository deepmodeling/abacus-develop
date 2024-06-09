#ifndef GINT_FORCE_GPU_H
#define GINT_FORCE_GPU_H

#include "module_hamilt_lcao/module_gint/gint.h"
#include "module_hamilt_lcao/module_gint/grid_technique.h"
namespace GintKernel
{

typedef struct
{
    int stream_num;
    double* input_dou;
    int* input_int;
    int* num_psir;
    int* atom_pair_A_m;
    int* atom_pair_B_n;
    int* atom_pair_K;
    int* atom_pair_lda;
    int* atom_pair_ldb;
    int* atom_pair_ldc;
    double* input_double_g;
    int* input_int_g;
    int* num_psir_g;
    double* psir_dm_device;
    double* psir_r_device;
    double* psir_lx_device;
    double* psir_ly_device;
    double* psir_lz_device;
    double* psir_lxx_device;
    double* psir_lxy_device;
    double* psir_lxz_device;
    double* psir_lyy_device;
    double* psir_lyz_device;
    double* psir_lzz_device;
    int* A_m_device;
    int* B_n_device;
    int* K_device;
    int* lda_device;
    int* ldb_device;
    int* ldc_device;
    double** matrix_A;
    double** matrix_B;
    double** matrix_C;
    double** matrix_A_device;
    double** matrix_B_device;
    double** matrix_C_device;
} grid_para;

typedef struct
{
    double* stress_device;
    double* force_device;
    int* iat_device;
    std::vector<int> iat_host;

} frc_strs_iat;

typedef struct
{
    double* stress_global;
    double* force_global;
    int* iat_global;
} frc_strs_iat_gbl;

typedef struct
{
    std::vector<double> density_mat_h;
    double* density_mat_d;
} DensityMat;

/**
 * @brief Calculate forces using GPU.
 *
 * This function calculates forces and stress for a given set of parameters.
 *
 * @param dm A pointer to hamilt::HContainer<double>.
 * @param vfactor Scaling factor for forces.
 * @param vlocal Local potential values.
 * @param force Output array for forces.
 * @param stress Output array for stress.
 * @param nczp Size parameter.
 * @param ylmcoef_now Coefficients for spherical harmonics.
 * @param gridt Reference to Grid_Technique object.
 */
void gint_fvl_gamma_gpu(hamilt::HContainer<double>* dm,
                          const double vfactor,
                          const double* vlocal,
                          std::vector<double>& force,
                          std::vector<double>& stress,
                          const int nczp,
                          double dr,
                          double* rcut,
                          const int isforce,
                          const int isstress,
                          const Grid_Technique& gridt,
                          const UnitCell& ucell);

/**
 * @brief GPU task generator for forces.
 *
 * This function generates GPU tasks for force calculations.
 *
 * @param gridt Reference to Grid_Technique .
 * @param ucell Reference to UnitCell .
 * @param grid_index_ij Index of the grid.
 * @param psi_size_max Maximum size of psi.
 * @param max_size Maximum size of atoms on a grid.
 * @param nczp Size parameter,stand for the current z-axis grids.
 * @param vfactor Scaling factor,stand for the Local potential.
 * @param rcut distance for each atom orbits
 * @param vlocal_global_value Global values of local potential.
 * @param iat_per_nbz save the number of the iat on per nbz grids.
 * @param atom_pair_num Number of atom pairs,stand for the max number of mat_n.
 * @param gpu_mat_cal_flag Establish whether to perform calculations between
 * atoms and grid points.
 * @param para Grid parameter in task generator,
 */

void gpu_task_generator_force(const Grid_Technique& gridt,
                              const UnitCell& ucell,
                              const int grid_index_ij,
                              const int psi_size_max,
                              const int max_size,
                              const int nczp,
                              const double vfactor,
                              const double* ruct,
                              const double* vlocal_global_value,
                              std::vector<int>& iat_per_nbz,
                              int& atom_pair_num,
                              std::vector<bool>& gpu_mat_cal_flag,
                              grid_para& para);
/**
 * @brief Calculate atom pair parameters.
 *
 * This function calculates the parameters for atom pairs.
 * @param gridt Reference to Grid_Technique.
 * @param ucell Reference to UnitCell.
 * @param grid_index_ij Index of the grid.
 * @param max_size Maximum size of atoms on a grid.
 * @param lgd Value of local grid dimension.
 * @param dm_matrix_g GPU array for dm_matrix.
 * @param max_m The maximum length of matrix A.
 * @param max_n The maximum length of matrix B.
 * @param gpu_mat_cal_flag Establish whether to perform calculations between
 *  atoms and grid points
 * @param para Grid parameter in multi matrix multiplication.
 */
void alloc_mult_force(const Grid_Technique& gridt,
                                    const UnitCell& ucell,
                                    const int grid_index_ij,
                                    const int max_size,
                                    const int lgd,
                                    double* dm_matrix_g,
                                    int& max_m,
                                    int& max_n,
                                    int& atom_pair_num,
                                    std::vector<bool>& gpu_mat_cal_flag,
                                    grid_para& para);
/**
 * @brief Density Matrix,force Stress Iat Init
 *
 * Using structure to init the parameter
 *
 * @param denstiy_mat DensityMat,contained the density_mat_dice and
 * destiyMatHost
 * @param f_s_iat_dev frc_strs_iat_gbl,contined the Force Stress and
 * Iat Number
 * @param dm hamilt::HContainer,denstiy stored in the Hcontainer
 * @param gridt Grid_Technique,stored the major method in the the gint.
 * @param UnitCell ucell,stored the cell tools
 * @param lgd Value of lgd,stand for the local grid dimension.
 * @param cuda_block in stress compute,used for Block nums
 * @param atom_num_grid in force calculate,used for Block nums
 */
void calculateInit(DensityMat& denstiy_mat,
                   frc_strs_iat_gbl& f_s_iat_dev,
                   hamilt::HContainer<double>* dm,
                   const Grid_Technique& gridt,
                   const UnitCell& ucell,
                   const int lgd,
                   const int cuda_block,
                   const int nat,
                   const int atom_num_grid);

/**
 * @brief Density Matrix,from Hcontainer to structure
 *
 * Using structure to init the parameter
 *
 * @param matrix_host double *,contained the destiyMatHost
 * @param dm hamilt::HContainer,denstiy stored in the Hcontainer
 * @param gridt Grid_Technique,stored the major method in the the gint.
 * @param lgd Value of lgd,stand for the local grid dimension.
 */
void allocateDm(std::vector<double>* matrix_host,
                hamilt::HContainer<double>* dm,
                const Grid_Technique& gridt,
                const UnitCell& ucell);

/**
 * @brief grid parameter Init
 *
 * GridParameter init
 *
 * @param para double *,contained the destiyMatHost
 * @param iter_num int , used for calcute the stream
 * @param nbz int,stand for the number of Z-axis
 * @param gridt Grid_Technique,stored the major method in the the gint.
 */
void para_init(grid_para& para,
                       const int iter_num,
                       const int nbz,
                       const int pipeline_index,
                       const Grid_Technique& gridt);
/**
 * @brief frc_strs_iat on host and device Init
 *
 * GridParameter init
 *
 * @param frc_strs_iat frc_strs_iat,contains the Force Stree Iat on Host
 * @param stream_num int , record the stream in GPU
 * @param cuda_block in stress compute,used for Block nums
 * @param atom_num_grid in force calculate,used for Block nums
 * @param max_size Maximum size of atoms on a grid.
 * @param frc_strs_iat_gbl frc_strs_iat_gbl,contains the Force Stree Iat on Host
 */
void cal_init(frc_strs_iat& f_s_iat,
                        const int stream_num,
                        const int cuda_block,
                        const int atom_num_grid,
                        const int nat,
                        const int max_size,
                        frc_strs_iat_gbl& f_s_iatg);
/**
 * @brief GridParameter memCpy,from Host to Device
 *
 * parameter init,which contains the gpu task and multi matrix multiplication
 *
 * @param para Grid parameter in task generator.
 * @param gridt Grid_Technique,stored the major method in the the gint.
 * @param nbz int,stand for the number of Z-axis.
 * @param atom_num_grid in force calculate,used for Block nums.
 */
void mem_copy(grid_para& para,
                    frc_strs_iat& f_s_iat,
                    const Grid_Technique& gridt,
                    const int nbz,
                    const int pipeline_index,
                    const int atom_num_grid);

} // namespace GintKernel
#endif
