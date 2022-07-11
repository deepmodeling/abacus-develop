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

class Gint
{
    public:

    // the unified interface to grid integration
	void cal_gint(Gint_inout *inout);

    // preparing FFT grid
    void prep_grid(
        const int &nbx_in,
        const int &nby_in,
        const int &nbz_in,
        const int &nbz_start_in,
        const int& ncxyz_in);

    protected:
    // variables related to FFT grid
 	int nbx;
	int nby;
	int nbz;
	int ncxyz;
	int nbz_start;

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
        double* pvpR_reduced);

    void gint_kernel_vlocal_meta(
        const int na_grid,
        const int grid_index,
        const double delta_r,
        double* vldr3,
        double* vkdr3,
        const int LD_pool,
        double* pvpR_reduced);

	void cal_meshball_vlocal_gamma(
		const int na_grid,  						// how many atoms on this (i,j,k) grid
		const int LD_pool,
		const int*const block_iw,					// block_iw[na_grid],	index of wave functions for each block
		const int*const block_size, 				// block_size[na_grid],	number of columns of a band
		const int*const block_index,				// block_index[na_grid+1], count total number of atomis orbitals
		const bool*const*const cal_flag,			// cal_flag[GlobalC::bigpw->bxyz][na_grid],	whether the atom-grid distance is larger than cutoff
		const double*const*const psir_ylm,			// psir_ylm[GlobalC::bigpw->bxyz][LD_pool]
		const double*const*const psir_vlbr3,		// psir_vlbr3[GlobalC::bigpw->bxyz][LD_pool]
		double* GridVlocal);		// GridVlocal[lgd_now][lgd_now]

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
        double**DM_R,
        const bool isforce,
        const bool isstress,
        ModuleBase::matrix* fvl_dphi,
        ModuleBase::matrix* svl_dphi);

    void cal_meshball_force(
        const int grid_index,
        const int na_grid,  					    // how many atoms on this (i,j,k) grid
        const int*const block_size, 			    // block_size[na_grid],	number of columns of a band
        const int*const block_index,		    	// block_index[na_grid+1], count total number of atomis orbitals
        const double*const*const psir_vlbr3_DMR,	    // psir_vlbr3[GlobalC::bigpw->bxyz][LD_pool]
        const double*const*const dpsir_x,	    // psir_vlbr3[GlobalC::bigpw->bxyz][LD_pool]
        const double*const*const dpsir_y,	    // psir_vlbr3[GlobalC::bigpw->bxyz][LD_pool]
        const double*const*const dpsir_z,	    // psir_vlbr3[GlobalC::bigpw->bxyz][LD_pool]
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
    // calculate the charge density via grid integrals
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

    //------------------------------------------------------
    // in gint_k_rho.cpp 
    //------------------------------------------------------
    // calculate the charge density via grid integrals
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

    void cal_meshball_crosstaus(
        const int na_grid,
        int* block_index,
        int* vindex,
        double** dpsix,
        double** dpsiy,
        double** dpsiz,
        double** dpsix_dm,
        double** dpsiy_dm,
        double** dpsiz_dm,
        double** crosstaus);

    // dimension: [GlobalC::LNNR.nnrg] 
    // save the < phi_0i | V | phi_Rj > in sparse H matrix.
    bool pvpR_alloc_flag = false;
    double** pvpR_reduced = nullptr; //stores Hamiltonian in reduced format, for multi-l
    
	double* pvpR_grid = nullptr; //stores Hamiltonian in grid format, for gamma-point
};

#endif