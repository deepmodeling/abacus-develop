#ifndef GINT_INTERFACE
#define GINT_INTERFACE

//This class provides a unified interface to the
// grid intergration operation used to calculate
// electron density, and the contribution of local potential
// to Hamiltonian and force/stress
// There are two derived classes of this class
// namely Gint_k and Gint_Gamma, which contains some
// specific operations for gamma point/multi-k calculations

#include "gint_tools.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"

class Gint
{
    public:

    ~Gint();

    // the unified interface to grid integration
	void cal_gint(Gint_inout *inout);

    // preparing FFT grid
    void prep_grid(
        const Grid_Technique& gt,
        const int& nbx_in,
        const int &nby_in,
        const int &nbz_in,
        const int &nbz_start_in,
        const int& ncxyz_in,
        const int& bx_in,
        const int& by_in,
        const int& bz_in,
        const int& bxyz_in,
        const int& nbxx_in,
        const int& ny_in,
        const int& nplane_in,
        const int& startz_current_in);
    
    /**
     * @brief calculate the neighbor atoms of each atom in this processor
     * size of BaseMatrix with be the non-parallel version
    */
    void initialize_pvpR(
        const UnitCell& unitcell,
        Grid_Driver* gd
    );
    
    /**
     * @brief transfer DMR (2D para) to DMR (Grid para) in elecstate_lcao.cpp
    */
    void transfer_DM2DtoGrid(std::vector<hamilt::HContainer<double>*> DM2D);

    const Grid_Technique* gridt = nullptr;
    protected:
    // variables related to FFT grid
 	int nbx;
	int nby;
	int nbz;
	int ncxyz;
    int nbz_start;
    int bx, by, bz, bxyz;
    int nbxx;
    int ny, nplane, startz_current; // from rhopw

    //------------------------------------------------------
    // in gint_vl.cpp 
    //------------------------------------------------------
    // calculate the matrix elements of Hamiltonian matrix,
    // < phi_0 | Vl + Vh + Vxc | phi_R> or if the Vna is used,
    // < phi_0 | delta_Vh + Vxc | phi_R>.
    void gint_kernel_vlocal(
        const int na_grid,
        const int grid_index,
        const double delta_r,
        double* vldr3,
        const int LD_pool,
        double* pvpR_reduced,
        hamilt::HContainer<double>* hR = nullptr);

    // calculate < phi_0 | vlocal | dphi_R >
    void gint_kernel_dvlocal(
        const int na_grid,
        const int grid_index,
        const double delta_r,
        double* vldr3,
        const int LD_pool,
        double* pvdpRx_reduced,
        double* pvdpRy_reduced,
        double* pvdpRz_reduced);

    void gint_kernel_vlocal_meta(
        const int na_grid,
        const int grid_index,
        const double delta_r,
        double* vldr3,
        double* vkdr3,
        const int LD_pool,
        double* pvpR_reduced,
        hamilt::HContainer<double>* hR = nullptr);

	void cal_meshball_vlocal_gamma(
		const int na_grid,  						// how many atoms on this (i,j,k) grid
		const int LD_pool,
		const int*const block_iw,					// block_iw[na_grid],	index of wave functions for each block
		const int*const block_size, 				// block_size[na_grid],	number of columns of a band
		const int*const block_index,				// block_index[na_grid+1], count total number of atomis orbitals
		const int grid_index,                       // index of grid group, for tracing iat
        const bool*const*const cal_flag,			// cal_flag[bxyz][na_grid],	whether the atom-grid distance is larger than cutoff
		const double*const*const psir_ylm,			// psir_ylm[bxyz][LD_pool]
		const double*const*const psir_vlbr3,		// psir_vlbr3[bxyz][LD_pool]
		hamilt::HContainer<double>* hR);		// HContainer for storing the <phi_0 | V | phi_R> matrix element.

    void cal_meshball_vlocal_k(
        int na_grid,
        int LD_pool,
        int grid_index, 
        int* block_size,
        int* block_index,
        int* block_iw,
        bool** cal_flag, 
        double** psir_ylm,
        double** psir_vlbr3,
        double* pvpR);

    //------------------------------------------------------
    // in gint_fvl.cpp 
    //------------------------------------------------------
    // calculate vl contributuion to force & stress via grid integrals
    void gint_kernel_force(
        const int na_grid,
        const int grid_index,
        const double delta_r,
        double* vldr3,
        const int LD_pool,
        double** DM_R,
        const int is,
        const bool isforce,
        const bool isstress,
        ModuleBase::matrix* fvl_dphi,
        ModuleBase::matrix* svl_dphi);

    void gint_kernel_force_meta(
        const int na_grid,
        const int grid_index,
        const double delta_r,
        double* vldr3,
        double* vkdr3,
        const int LD_pool,
        double** DM_in,
        const int is,
        const bool isforce,
        const bool isstress,
        ModuleBase::matrix* fvl_dphi,
        ModuleBase::matrix* svl_dphi);

    void cal_meshball_force(
        const int grid_index,
        const int na_grid,  					    // how many atoms on this (i,j,k) grid
        const int*const block_size, 			    // block_size[na_grid],	number of columns of a band
        const int*const block_index,		    	// block_index[na_grid+1], count total number of atomis orbitals
        const double*const*const psir_vlbr3_DMR,	    // psir_vlbr3[bxyz][LD_pool]
        const double*const*const dpsir_x,	    // psir_vlbr3[bxyz][LD_pool]
        const double*const*const dpsir_y,	    // psir_vlbr3[bxyz][LD_pool]
        const double*const*const dpsir_z,	    // psir_vlbr3[bxyz][LD_pool]
        ModuleBase::matrix *force);

    void cal_meshball_stress(
        const int na_grid,  					    // how many atoms on this (i,j,k) grid
        const int*const block_index,		    	// block_index[na_grid+1], count total number of atomis orbitals
        const double*const*const psir_vlbr3_DMR,
        const double*const*const dpsir_xx,
        const double*const*const dpsir_xy,
        const double*const*const dpsir_xz,
        const double*const*const dpsir_yy,
        const double*const*const dpsir_yz,
        const double*const*const dpsir_zz,
        ModuleBase::matrix *stress);

    //------------------------------------------------------
    // in gint_k_rho.cpp 
    //------------------------------------------------------
    // calculate the charge density & kinetic energy density (tau) via grid integrals
    void gint_kernel_rho(
        const int na_grid,
        const int grid_index,
        const double delta_r,
        int* vindex,
        const int LD_pool,
        Gint_inout *inout);

    void cal_meshball_rho(
        const int na_grid,
        int* block_index,
        int* vindex,
        double** psir_ylm,
        double** psir_DMR,
        double* rho);

    void gint_kernel_tau(
        const int na_grid,
        const int grid_index,
        const double delta_r,
        int* vindex,
        const int LD_pool,
        Gint_inout *inout);

    void cal_meshball_tau(
        const int na_grid,
        int* block_index,
        int* vindex,
        double** dpsix,
        double** dpsiy,
        double** dpsiz,
        double** dpsix_dm,
        double** dpsiy_dm,
        double** dpsiz_dm,
        double* rho);

    // dimension: [GlobalC::LNNR.nnrg] 
    // save the < phi_0i | V | phi_Rj > in sparse H matrix.
    bool pvpR_alloc_flag = false;
    double** pvpR_reduced = nullptr; //stores Hamiltonian in reduced format, for multi-l
    hamilt::HContainer<double>* hRGint = nullptr; //stores Hamiltonian in sparse format
    hamilt::HContainer<std::complex<double>>* hRGintCd = nullptr; //stores Hamiltonian in sparse format
    std::vector<hamilt::HContainer<double>*> DMRGint; //stores DMR in sparse format
    hamilt::HContainer<double>* DMRGint_full = nullptr; //tmp tools used in transfer_DM2DtoGrid
    double** pvdpRx_reduced = nullptr;
    double** pvdpRy_reduced = nullptr;
    double** pvdpRz_reduced = nullptr;

	double* pvpR_grid = nullptr; //stores Hamiltonian in grid format, for gamma-point
};

#endif