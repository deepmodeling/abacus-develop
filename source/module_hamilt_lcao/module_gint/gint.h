#ifndef GINT_INTERFACE
#define GINT_INTERFACE

#include "gint_tools.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_hamilt_lcao/module_gint/grid_technique.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"
#include <functional>

//----------------------------------------------------------
//！This class provides a unified interface to the
//！grid intergration operation used to calculate
//！electron density, and the contribution of local
//！potential to Hamiltonian and force/stress.
//！There are two derived classes of this class
//! namely Gint_Gamma and Gint_k, which contain
//! specific operations for gamma point/multi-k calculations
//----------------------------------------------------------

class Gint {
  public:
    ~Gint();

    //! move operator for the next ESolver to directly use its infomation
    Gint& operator=(Gint&& rhs);

    hamilt::HContainer<double>* get_hRGint() const { return hRGint; }

    std::vector<hamilt::HContainer<double>*> get_DMRGint() const { return DMRGint; }

    int get_ncxyz() const { return ncxyz; }

    //! the unified interface to grid integration
    void cal_gint(Gint_inout* inout);

    //! preparing FFT grid
    void prep_grid(const Grid_Technique& gt,
                   const int& nbx_in,
                   const int& nby_in,
                   const int& nbz_in,
                   const int& nbz_start_in,
                   const int& ncxyz_in,
                   const int& bx_in,
                   const int& by_in,
                   const int& bz_in,
                   const int& bxyz_in,
                   const int& nbxx_in,
                   const int& ny_in,
                   const int& nplane_in,
                   const int& startz_current_in,
                   const UnitCell* ucell_in,
                   const LCAO_Orbitals* orb_in);

    /**
     * @brief calculate the neighbor atoms of each atom in this processor
     * size of BaseMatrix with be the non-parallel version
     */
    void initialize_pvpR(const UnitCell& unitcell, Grid_Driver* gd, const int& nspin);

    /**
     * @brief resize DMRGint to nspin and reallocate the memory
     */
    void reset_DMRGint(const int& nspin);

    /**
     * @brief transfer DMR (2D para) to DMR (Grid para) in elecstate_lcao.cpp
     */
    void transfer_DM2DtoGrid(std::vector<hamilt::HContainer<double>*> DM2D);

    const Grid_Technique* gridt = nullptr;
    const UnitCell* ucell;

    // psir_ylm_new = psir_func(psir_ylm)
    // psir_func==nullptr means psir_ylm_new=psir_ylm
    using T_psir_func = std::function<
        const ModuleBase::Array_Pool<double>&(
            const ModuleBase::Array_Pool<double> &psir_ylm,
            const Grid_Technique &gt,
            const int grid_index,
            const int is,
            const std::vector<int> &block_iw,
            const std::vector<int> &block_size,
            const std::vector<int> &block_index,
            const ModuleBase::Array_Pool<bool> &cal_flag)>;

    T_psir_func psir_func_1 = nullptr;
    T_psir_func psir_func_2 = nullptr;

  protected:

    //! variables related to FFT grid
    int nbx;
    int nby;
    int nbz;
    int ncxyz;
    int nbz_start;
    int bx;
    int by;
    int bz;
    int bxyz;
    int nbxx;
    int ny;
    int nplane;
    int startz_current; // from rhopw

    //! in cal_gint_gpu.cpp
    void gpu_vlocal_interface(Gint_inout* inout);

    void gpu_rho_interface(Gint_inout* inout);

    void gpu_force_interface(Gint_inout* inout);

    //! in cal_gint_cpu.cpp
    void gint_kernel_vlocal(Gint_inout* inout);

    //! calculate H_mu_nu(local)=<phi_0|vlocal|dphi_R>
    void gint_kernel_dvlocal(Gint_inout* inout);

    //! calculate vlocal in meta-GGA functionals
    void gint_kernel_vlocal_meta(Gint_inout* inout);

    //! calculate charge density rho(r)=\int D_munu \phi_mu \phi_nu
    void gint_kernel_rho(Gint_inout* inout);

    //! used in meta-GGA functional
    void gint_kernel_tau(Gint_inout* inout);

    //! compute forces
    void gint_kernel_force(Gint_inout* inout);

    //! compute forces related to meta-GGA functionals
    void gint_kernel_force_meta(Gint_inout* inout);

    //! calculate local potential contribution to the Hamiltonian
    //! na_grid: how many atoms on this (i,j,k) grid
    //! block_iw: dim is [na_grid], index of wave function for each block
    //! block_size: dim is [block_size], number of columns of a band
    //! block_index: dim is [na_grid+1], total number of atomic orbitals
    //! grid_index: index of grid group, for tracing iat
    //! cal_flag: dim is [bxyz][na_grid], whether the atom-grid distance is larger than cutoff
    //! psir_ylm: dim is [bxyz][LD_pool]
    //! psir_vlbr3: dim is [bxyz][LD_pool]
    //! hR: HContainer for storing the <phi_0|V|phi_R> matrix elements

    void cal_meshball_vlocal(
        const int na_grid,
        const int LD_pool,
        const int* const block_iw,
        const int* const block_size,
        const int* const block_index,
        const int grid_index,
        const bool* const* const cal_flag,
        const double* const* const psir_ylm,
        const double* const* const psir_vlbr3,
        hamilt::HContainer<double>* hR);


    //! in gint_fvl.cpp
    //! calculate vl contributuion to force & stress via grid integrals
    void gint_kernel_force(const int na_grid,
                           const int grid_index,
                           const double delta_r,
                           double* vldr3,
                           const int is,
                           const bool isforce,
                           const bool isstress,
                           ModuleBase::matrix* fvl_dphi,
                           ModuleBase::matrix* svl_dphi,
                           const UnitCell& ucell);

    //! in gint_fvl.cpp
    //! calculate vl contributuion to force & stress via grid integrals
    //! used in meta-GGA calculations
    void gint_kernel_force_meta(const int na_grid,
                                const int grid_index,
                                const double delta_r,
                                double* vldr3,
                                double* vkdr3,
                                const int is,
                                const bool isforce,
                                const bool isstress,
                                ModuleBase::matrix* fvl_dphi,
                                ModuleBase::matrix* svl_dphi,
                                const UnitCell& ucell);

    //! Use grid integrals to compute the atomic force contributions
    //! na_grid: how many atoms on this (i,j,k) grid
    //! block_size: dim is [na_grid], number of columns of a band
    //! block_index: dim is [na_grid+1], total number of atomis orbitals
    //! psir_vlbr3_DMR: dim is [bxyz][LD_pool]
    //! dpsir_x: dim is [bxyz][LD_pool]
    //! dpsir_y: dim is [bxyz][LD_pool]
    //! dpsir_z: dim is [bxyz][LD_pool]
    void cal_meshball_force(
        const int grid_index,
        const int na_grid,
        const int* const block_size,
        const int* const block_index,
        const double* const* const psir_vlbr3_DMR,
        const double* const* const dpsir_x,        // psir_vlbr3[bxyz][LD_pool]
        const double* const* const dpsir_y,        // psir_vlbr3[bxyz][LD_pool]
        const double* const* const dpsir_z,        // psir_vlbr3[bxyz][LD_pool]
        ModuleBase::matrix* force);

    //! Use grid integrals to compute the stress contributions
    //! na_grid: how many atoms on this (i,j,k) grid
    //! block_index: dim is [na_grid+1], total number of atomis orbitals
    void cal_meshball_stress(
        const int na_grid,
        const int*const block_index,
        const double*const psir_vlbr3_DMR,
        const double*const dpsirr,
        ModuleBase::matrix *stress);
    
    //! Use grid integrals to compute charge density
    //! in gint_k_rho.cpp
    //! calculate the charge density & kinetic energy density (tau) via grid integrals
    void gint_kernel_rho(const int na_grid,
                         const int grid_index,
                         const double delta_r,
                         int* vindex,
                         const int LD_pool,
                         const UnitCell& ucell,
                         Gint_inout* inout);

    //! Use grid integrals to compute charge density in a meshball
    void cal_meshball_rho(const int na_grid,
                          const int*const block_index,
                          const int*const vindex,
                          const double*const*const psir_ylm,
                          const double*const*const psir_DMR,
                          double*const rho);

    //! Use grid integrals to compute kinetic energy density tau 
    //！in meta-GGA functional 
    void gint_kernel_tau(const int na_grid,
                         const int grid_index,
                         const double delta_r,
                         int* vindex,
                         const int LD_pool,
                         Gint_inout* inout,
                         const UnitCell& ucell);

    //! Use grid integrals to compute kinetic energy density tau
    //！in a meshball, used in meta-GGA functional calculations
    void cal_meshball_tau(const int na_grid,
                          int* block_index,
                          int* vindex,
                          double** dpsix,
                          double** dpsiy,
                          double** dpsiz,
                          double** dpsix_dm,
                          double** dpsiy_dm,
                          double** dpsiz_dm,
                          double* rho);

    //! save the < phi_0i | V | phi_Rj > in sparse H matrix.
    //! stores Hamiltonian in sparse format
    hamilt::HContainer<double>* hRGint = nullptr; 

    //! size of vec is 4, only used when nspin = 4
    std::vector<hamilt::HContainer<double>*> hRGint_tmp; 

    //! stores Hamiltonian in sparse format
    hamilt::HContainer<std::complex<double>>* hRGintCd = nullptr; 

    //! stores DMR in sparse format
    std::vector<hamilt::HContainer<double>*> DMRGint; 

    //! tmp tools used in transfer_DM2DtoGrid 
    hamilt::HContainer<double>* DMRGint_full = nullptr;

    std::vector<hamilt::HContainer<double>> pvdpRx_reduced;
    std::vector<hamilt::HContainer<double>> pvdpRy_reduced;
    std::vector<hamilt::HContainer<double>> pvdpRz_reduced;
};

#endif
