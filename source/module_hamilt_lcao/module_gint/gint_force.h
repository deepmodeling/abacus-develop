#ifndef GINT_FORCE_H
#define GINT_FORCE_H

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
    double* psi_input_double_g;
    int* input_int_g;
    int* num_psirDevice;
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
} SGridParameter;

typedef struct
{
    double* stress_device;
    double* stress_host;
    double* force_device;
    double* force_host;
    int* iat_device;
    int* iat_host;

} ForceStressIat;

typedef struct
{
    double* stress_global;
    double* force_global;
    int* iat_global;
} ForceStressIatGlobal;

typedef struct
{
    double* density_mat_h;
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
void gint_gamma_force_gpu(hamilt::HContainer<double>* dm,
                          const double vfactor,
                          const double* vlocal,
                          double* force,
                          double* stress,
                          const int nczp,
                          double dr,
                          double* rcut,
                          const Grid_Technique& gridt,
                          const UnitCell& ucell);

/**
 * @brief GPU task generator for forces.
 *
 * This function generates GPU tasks for force calculations.
 *
 * @param gridt Reference to Grid_Technique object.
 * @param i Value of i,stand for the x-axis gird.
 * @param j Value of j.stand for the y-axis grid.
 * @param psi_size_max Maximum size of psi.
 * @param max_size Maximum size of atoms on a grid.
 * @param nczp Size parameter,stand for the current z-axis grids.
 * @param vfactor Scaling factor,stand for the Local potential.
 * @param rcut distance for each atom orbits
 * @param vlocal_global_value Global values of local potential.
 * @param iat_per_nbz save the number of the iat on per nbz grids.
 * @param lgd Value of lgd,stand for the local grid dimension.
 * @param num_psir Array for num_psir values,contained the each number of the
 * atom psir on a grid.
 * @param dm_matrix_g GPU array for dm_matrix,send as the denstiy matrix.
 * @param max_m Maximum value of m,stand for the max number of mat_m.
 * @param max_n Maximum value of n,stand for the max number of mat_n.
 * @param atom_pair_num Number of atom pairs,stand for the max number of mat_n.
 * @param para Grid parameter in task generator,
 */

void gpu_task_generator_force(const Grid_Technique& gridt,
                              const UnitCell& ucell,
                              const int i,
                              const int j,
                              const int psi_size_max,
                              const int max_size,
                              const int nczp,
                              const double vfactor,
                              double* ruct,
                              const double* vlocal_global_value,
                              int* iat_per_nbz,
                              const int lgd,
                              double* dm_matrix_g,
                              int& max_m,
                              int& max_n,
                              int& atom_pair_num,
                              SGridParameter& para);
/**
 * @brief Density Matrix,force Stress Iat Init
 *
 * Using structure to init the parameter
 *
 * @param denstiy_mat DensityMat,contained the density_mat_dice and
 * destiyMatHost
 * @param calcualteG ForceStressIatGlobal,contined the Force Stress and
 * Iat Number
 * @param dm hamilt::HContainer,denstiy stored in the Hcontainer
 * @param gridt Grid_Technique,stored the major method in the the gint.
 * @param UnitCell ucell,stored the cell tools
 * @param lgd Value of lgd,stand for the local grid dimension.
 * @param cuda_block in stress compute,used for Block nums
 * @param atom_num_grid in force calculate,used for Block nums
 */
void calculateInit(DensityMat& denstiy_mat,
                   ForceStressIatGlobal& calcualteG,
                   hamilt::HContainer<double>* dm,
                   const Grid_Technique& gridt,
                   const UnitCell& ucell,
                   int lgd,
                   int cuda_block,
                   int atom_num_grid);

/**
 * @brief Density Matrix,from Hcontainer to structure
 *
 * Using structure to init the parameter
 *
 * @param MatrixHost double *,contained the destiyMatHost
 * @param dm hamilt::HContainer,denstiy stored in the Hcontainer
 * @param gridt Grid_Technique,stored the major method in the the gint.
 * @param lgd Value of lgd,stand for the local grid dimension.
 */
void allocateDm(double* MatrixHost,
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
void para_init(SGridParameter& para,
                       int iter_num,
                       int nbz,
                       const Grid_Technique& gridt);
/**
 * @brief ForceStressIat on host and device Init
 *
 * GridParameter init
 *
 * @param ForceStressIat ForceStressIat,contains the Force Stree Iat on Device
 * and Host
 * @param stream_num int , record the stream in GPU
 * @param cuda_block in stress compute,used for Block nums
 * @param atom_num_grid in force calculate,used for Block nums
 * @param max_size Maximum size of atoms on a grid.
 * @param stress_global sorted in the ForceStressIatG,the array of global stress
 * @param force_global sorted in the ForceStressIatG,the array of global force
 * @param iat_global sorted in the ForceStressIatG,the array of global iat
 */
void cal_init(ForceStressIat& calcualte,
                        int stream_num,
                        int cuda_block,
                        int atom_num_grid,
                        int max_size,
                        ForceStressIatGlobal& calcualteg);
/**
 * @brief GridParameter memCpy,from Host to Device
 *
 * parameter init,which contains the gpu task and multi matrix multiplication
 *
 * @param para Grid parameter in task generator,
 * @param gridt Grid_Technique,stored the major method in the the gint.
 * @param nbz int,stand for the number of Z-axis
 * @param atom_num_grid in force calculate,used for Block nums
 */
void para_mem_copy(SGridParameter& para,
                         const Grid_Technique& gridt,
                         int nbz,
                         int atom_num_grid);
/**
 * @brief Force Stress Force Iat memCpy,from Host to Device
 *
 *  @param ForceStressIat ForceStressIat,contains the Force Stree Iat on Device
 * and Host
 *  @param gridt Grid_Technique,stored the major method in the the gint.
 *  @param atom_num_grid in force calculate,used for Block nums
 *  @param cuda_block in stress compute,used for Block nums
 *  @param stream_num int , record the stream in GPU
 */
void cal_mem_cpy(ForceStressIat& calcualte,
                          const Grid_Technique& gridt,
                          int atom_num_grid,
                          int cuda_block,
                          int stream_num);
/**
 * @brief Force Calculate on Host
 *
 * @param ForceStressIat ForceStressIat,contains the Force Stree Iat on Device
 * and Host
 * @param force stored the force for each atom on each directions
 * @param atom_num_grid in force calculate,used for Block nums
 */
void cal_calculate_cpu(ForceStressIat& calcualte,
                    double* force,
                    int atom_num_grid);
/**
 * @brief Stress Calculate on Host
 *
 * @param ForceStressIat ForceStressIat,contains the Force Stree Iat on Device
 * and Host
 * @param stress stored the stress for each directions
 * @param cuda_block in stress compute,used for Block nums
 */
void cal_stress_cpu(ForceStressIat& calcualte,
                     double* stress,
                     int cuda_block);
} // namespace GintKernel
#endif
