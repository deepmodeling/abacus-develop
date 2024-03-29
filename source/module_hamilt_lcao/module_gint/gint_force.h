#ifndef GINT_FORCE_H
#define GINT_FORCE_H

#include "module_hamilt_lcao/module_gint/gint.h"
#include "module_hamilt_lcao/module_gint/grid_technique.h"
namespace GintKernel{

typedef struct 
{
    int      streamNum;
    double  *input_dou;
    int     *input_int;
    int     *num_psir;
    int     *atom_pair_A_m;
    int     *atom_pair_B_n;
    int     *atom_pair_K;
    int     *atom_pair_lda;
    int     *atom_pair_ldb;
    int     *atom_pair_ldc;
    double  *psi_input_double_g;
    int     *input_int_g;
    int     *num_psirDevice ;
    double  *psir_dm_device;
    double  *psir_r_device ;
    double  *psir_lx_device ;
    double  *psir_ly_device ;
    double  *psir_lz_device ;
    double  *psir_lxx_device;
    double  *psir_lxy_device;
    double  *psir_lxz_device;
    double  *psir_lyy_device;
    double  *psir_lyz_device;
    double  *psir_lzz_device;
    int     *A_m_device;
    int     *B_n_device;
    int     *K_device;
    int     *lda_device;
    int     *ldb_device;
    int     *ldc_device;
    double  **matrix_A;
    double  **matrix_B;
    double  **matrix_C;
    double  **matrix_ADev;
    double  **matrix_BDev;
    double  **matrix_CDev;
}SGridParameter;

typedef struct 
{
    double  *stressDev;
    double  *stressHost;
    double  *forceDev;
    double  *forceHost;
    int     *iatDev;
    int     *iatHost;

}ForceStressIat;

typedef struct 
{
    double  *stressGlobal;
    double  *forceGlobal;
    int     *iatGlobal;
}ForceStressIatGlobal;

typedef struct 
{
    double  *densityMatHost;
    double  *densityMatDev;
}DensityMat;

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
void gint_gamma_force_gpu(hamilt::HContainer<double> *dm,
                          const double vfactor,
                          const double *vlocal,
                          double *force, 
                          double *stress,
                          const int nczp, 
                          const double *ylmcoef_now,
                          const Grid_Technique &gridt,
                          const LCAO_Orbitals &ORB,
                          const UnitCell &ucell);

/**
 * @brief GPU task generator for forces.
 *
 * This function generates GPU tasks for force calculations.
 *
 * @param gridt Reference to Grid_Technique object.
 * @param i Value of i.
 * @param j Value of j.
 * @param psi_size_max Maximum size of psi.
 * @param max_size Maximum size.
 * @param nczp Size parameter.
 * @param vfactor Scaling factor.
 * @param vlocal_global_value Global values of local potential.
 * @param iat_per_nbz Array of iat values per nbz.
 * @param psi_input_double Double array for input psi values.
 * @param input_int Integer array for input psi values.
 * @param num_psir Array for num_psir values.
 * @param lgd Value of lgd.
 * @param psir_ylm_g GPU array for psir_ylm.
 * @param psir_zeros_g GPU array for psir_zeros.
 * @param dm_matrix_g GPU array for dm_matrix.
 * @param mat_m Array for mat_m values.
 * @param mat_n Array for mat_n values.
 * @param mat_k Array for mat_k values.
 * @param mat_lda Array for mat_lda values.
 * @param mat_ldb Array for mat_ldb values.
 * @param mat_ldc Array for mat_ldc values.
 * @param mat_A Double pointer for mat_A.
 * @param mat_B Double pointer for mat_B.
 * @param mat_C Double pointer for mat_C.
 * @param max_m Maximum value of m.
 * @param max_n Maximum value of n.
 * @param atom_pair_num Number of atom pairs.
 * @param rho_g GPU array for rho.
 * @param vec_l Double pointer for vec_l.
 * @param vec_r Double pointer for vec_r.
 * @param dot_product Double pointer for dot_product.
 * @param vec_len Array for vec_len values.
 * @param dot_count Reference to dot_count.
 */
void gtask_force(
    const Grid_Technique &gridt, const double *rcut,
    const UnitCell &ucell, const int i, const int j,
    const int psi_size_max, const int max_size, const int nczp,
    const double vfactor, const double *vlocal_global_value, int *iat_per_nbz,
    double *psi_input_double, int *input_int, int *num_psir, const int lgd,
    double *psir_ylm_g, double *psir_zeros_g, double *dm_matrix_g, int *mat_m,
    int *mat_n, int *mat_k, int *mat_lda, int *mat_ldb, int *mat_ldc,
    double **mat_A, double **mat_B, double **mat_C, int &max_m, int &max_n,
    int &atom_pair_num);
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
 * @param vlocal_global_value Global values of local potential.
 * @param iat_per_nbz save the number of the iat on per nbz grids.
 * @param input_dou Double array for input psi values,contains the x,y,z,distance of the grids.
 * @param input_int Integer array for input psi values,contains the index of the girds.
 * @param num_psir Array for num_psir values,contained the each number of the atom psir on a grid.
 * @param lgd Value of lgd,stand for the local grid dimension.
 * @param psir_ylm_g GPU array for psir_ylm,send as the right matrix.
 * @param psir_zeros_g GPU array for psir_zeros,send as the zero matirx.
 * @param dm_matrix_g GPU array for dm_matrix,send as the denstiy matrix.
 * @param mat_m Array for mat_m values,max number to choose on the psir_ylm_g.
 * @param mat_n Array for mat_n values.
 * @param mat_k Array for mat_k values.
 * @param mat_lda Array for mat_lda values,contained each leading dimension for psir_ylm_g.
 * @param mat_ldb Array for mat_ldb values,contained each leading dimension for dm_matrix_g.
 * @param mat_ldc Array for mat_ldc values,contained each leading dimension for psir_zeros_g.
 * @param mat_A Double pointer for mat_A,batch matrix to compute,according to the mat_m.
 * @param mat_B Double pointer for mat_B,batch matrix to compute,according to the mat_m.
 * @param mat_C Double pointer for mat_C.
 * @param max_m Maximum value of m,stand for the max number of mat_m.
 * @param max_n Maximum value of n,stand for the max number of mat_n.
 * @param atom_pair_num Number of atom pairs,stand for the max number of mat_n.
 */

void gpu_task_generator_force(const Grid_Technique &gridt, 
                        const LCAO_Orbitals &ORB,
                        const UnitCell &ucell, 
                        const int i, 
                        const int j,
                        const int psi_size_max, 
                        const int max_size, 
                        const int nczp,
                        const double vfactor, 
                        const double *vlocal_global_value, 
                        int *iat_per_nbz,
                        const int lgd, 
                        double *dm_matrix_g, 
                        int &max_m, 
                        int &max_n,
                        int &atom_pair_num,
                        SGridParameter &para);
/**
 * @brief Density Matrix,force Stress Iat Init
 *
 * Using structure to init the parameter
 *
 * @param densityMat DensityMat,contained the densityMatDevice and destiyMatHost
 * @param forceStressIatG ForceStressIatGlobal,contined the Force Stress and Iat
 * Number
 * @param dm hamilt::HContainer,denstiy stored in the Hcontainer
 * @param gridt Grid_Technique,stored the major method in the the gint.
 * @param UnitCell ucell,stored the cell tools
 * @param lgd Value of lgd,stand for the local grid dimension.
 * @param cudaBlocks in stress compute,used for Block nums
 * @param atomNumOnGrids in force calculate,used for Block nums
 */
void CalculateInit(DensityMat &densityMat,
                ForceStressIatGlobal &forceStressIatG,
                hamilt::HContainer<double> *dm,
                const Grid_Technique &gridt, const UnitCell &ucell,
                int lgd,int cudaBlocks,int atomNumOnGrids);

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
void AllocateDm(double *MatrixHost,
                hamilt::HContainer<double> *dm,
                const Grid_Technique &gridt,
                const UnitCell &ucell);

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
void CalculateGridInit(SGridParameter &para,
                int iter_num,
                int nbz,
                const Grid_Technique &gridt);
/**
 * @brief ForceStressIat on host and device Init
 *
 * GridParameter init
 *
 * @param ForceStressIat ForceStressIat,contains the Force Stree Iat on Device
 * and Host
 * @param streamNum int , record the stream in GPU
 * @param cudaBlocks in stress compute,used for Block nums
 * @param atomNumOnGrids in force calculate,used for Block nums
 * @param max_size Maximum size of atoms on a grid.
 * @param stressGlobal sorted in the ForceStressIatG,the array of global stress
 * @param forceGlobal sorted in the ForceStressIatG,the array of global force
 * @param iatGlobal sorted in the ForceStressIatG,the array of global iat
 */
void ForceStressIatInit(ForceStressIat &forceStressIat,int streamNum,int cudaBlocks,int atomNumOnGrids,
                        int max_size,ForceStressIatGlobal &forceStressIatg );
/**
 * @brief GridParameter memCpy,from Host to Device
 *
 * parameter init,which contains the gpu task and multi matrix multiplication
 *
 * @param para Grid parameter in task generator,
 * @param gridt Grid_Technique,stored the major method in the the gint.
 * @param nbz int,stand for the number of Z-axis
 * @param atomNumOnGrids in force calculate,used for Block nums
 */
void CalculateGridMemCpy(SGridParameter &para,
                        const Grid_Technique &gridt,
                        int nbz,
                        int atomNumOnGrids);
/**
 * @brief Force Stress Force Iat memCpy,from Host to Device
 *
 *  @param ForceStressIat ForceStressIat,contains the Force Stree Iat on Device
 * and Host
 *  @param gridt Grid_Technique,stored the major method in the the gint.
 *  @param atomNumOnGrids in force calculate,used for Block nums
 *  @param cudaBlocks in stress compute,used for Block nums
 *  @param streamNum int , record the stream in GPU
 */
void ForceStressIatMemCpy(ForceStressIat &forceStressIat,
                        const Grid_Technique &gridt,
                        int atomNumOnGrids,int cudaBlocks,int streamNum);
/**
 * @brief Force Calculate on Host
 *
 * @param ForceStressIat ForceStressIat,contains the Force Stree Iat on Device
 * and Host
 * @param force stored the force for each atom on each directions
 * @param atomNumOnGrids in force calculate,used for Block nums
 */
void ForceCalculate(ForceStressIat &forceStressIat,
                    double *force,int atomNumOnGrids);
/**
 * @brief Stress Calculate on Host
 *
 * @param ForceStressIat ForceStressIat,contains the Force Stree Iat on Device
 * and Host
 * @param stress stored the stress for each directions
 * @param cudaBlocks in stress compute,used for Block nums
 */
void StressCalculate(ForceStressIat &forceStressIat,
                    double *stress,int cudaBlocks);
} // namespace GintKernel
#endif
