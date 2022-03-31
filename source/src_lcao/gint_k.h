#ifndef GINT_K_H
#define GINT_K_H

#include "gint_k_init.h"
#include "../module_orbital/ORB_atomic_lm.h"
#include "grid_technique.h"
#include "LCAO_matrix.h"

// add by jingan for map<> in 2021-12-2, will be deleted in the future
#include "../src_ri/abfs-vector3_order.h"

class Gint_k : public Gint_k_init
{
    public:

    Gint_k();
    ~Gint_k();

    // allocate the <phi_0 | V | phi_R> matrix element.
    void allocate_pvpR(void);
    void allocate_pvpR_tr(void); //LiuXh add 2019-07-15


    // destroy the temporary <phi_0 | V | phi_R> matrix element.
    void destroy_pvpR(void);
    //LiuXh add 2019-07-15
    void destroy_pvpR_tr(void);
    void distribute_pvpR_tr(void);

    // jingan add 2021-6-4, modify 2021-12-2
    void distribute_pvpR_sparseMatrix(
        const int current_spin, 
        const double &sparse_threshold, 
        const std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, double>>> &pvpR_sparseMatrix
    );
    void distribute_pvpR_soc_sparseMatrix(
        const double &sparse_threshold, 
        const std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, std::complex<double>>>> &pvpR_soc_sparseMatrix
    );
    void cal_vlocal_R_sparseMatrix(const int &current_spin, const double &sparse_threshold);

    // reset the spin.
    void reset_spin(const int &spin_now);

    // get the spin.
    int get_spin(void)const{return spin_now;}


    //------------------------------------------------------
    // in gint_k_vl.cpp 
    //------------------------------------------------------
    // calculate the matrix elements of Hamiltonian matrix,
    // < phi_0 | Vl + Vh + Vxc | phi_R> or if the Vna is used,
    // < phi_0 | delta_Vh + Vxc | phi_R>.
    void cal_vlocal_k(const double* vrs1, const Grid_Technique &gt, const int spin=0);

    //------------------------------------------------------
    // in gint_k.cpp 
    //------------------------------------------------------
    void cal_vlocal_R(const int current_spin); //LiuXh add 2019-07-15

    // folding the < phi_0 | V | phi_R> matrix to 
    // <phi_0i | V | phi_0j>
    // V is (Vl + Vh + Vxc) if no Vna is used,
    // and is (Vna + delta_Vh + Vxc) if Vna is used.
    void folding_vl_k(const int &ik);

    void folding_vl_k_nc(const int &ik);//zhengdy-soc

    // folding the < dphi_0 | V | phi_R> matrix to 
    // < dphi_0i | V | phi_0j>
    void folding_force(
        ModuleBase::matrix& fvl_dphi,
        double* pvdpx, 
        double* pvdpy, 
        double* pvdpz);//mohan add 2012-1-6

    // folding the < dphi_0 | V * R_beta | phi_R> matrix
    // < dphi_0i | V | phi_0j>
    void folding_stress(
        const bool isforce,
        const bool isstress,
        ModuleBase::matrix& fvl_dphi, 
        ModuleBase::matrix& svl_dphi,
        double* pvdpx, 
        double* pvdpy, 
        double* pvdpz,
        double* pvdp11, 
        double* pvdp22, 
        double* pvdp33, 
        double* pvdp12, 
        double* pvdp13, 
        double* pvdp23);//zhengdy add 2016-10-18

    //------------------------------------------------------
    // in gint_k_rho.cpp 
    //------------------------------------------------------
    // calculate the charge density via grid integrals
    void cal_rho_k(double** DM_R_in);

    //------------------------------------------------------
    // in gint_k_fvl.cpp 
    //------------------------------------------------------
    // calculate the force (many k-points).
    void fvl_k_RealSpace(
        ModuleBase::matrix& fvl_dphi, 
        const double* vl);//mohan add 2011-06-19

    void svl_k_RealSpace(
        const bool isforce,
        const bool isstress,
        ModuleBase::matrix& fvl_dphi, 
        ModuleBase::matrix& svl_dphi, 
        const double* vl);//zhengdy add 2016-10-18

    private:
    
    //------------------------------------------------------
    // in gint_k.cpp 
    //------------------------------------------------------
    // set the orbital info 
    // set the orbital/Ylm information on each real space grid.
    void set_ijk_atom(
        const int &grid_index, 
        const int &size,
        double*** psir_ylm, 
        double*** dr, 
        bool** cal_flag, 
        double** distance, 
        const double &delta_r);

    //------------------------------------------------------
    // in gint_k_vl.cpp 
    //------------------------------------------------------
    // evaluate the matrix element < phi0 | V | phiR> and store them in
    // a full H matrix.
    void evaluate_pvpR_full(
        const int &grid_index, 
        const int &size, 
        double*** psir_ylm,
        bool** cal_flag, 
        double* vldr3);

    // reduced means the H storage take the advance of adjacent atoms.
    void evaluate_pvpR_reduced(
        double* pvpR, 
        const int &grid_index, 
        const int &size, 
        const int &i, 
        const int &j, 
        const int &k,
        double*** psir_ylm, 
        bool** cal_flag, 
        double* vldr3, 
        double** distance, 
        const Grid_Technique &gt);

    //------------------------------------------------------
    // in gint_k_rho.cpp 
    //------------------------------------------------------
    // evaluate the <phi0 | Density Matrix | phiR> to get the charge density.
    void evaluate_pDMp(
        const int &grid_index, 
        const int &size,
        bool** cal_flag, 
        double*** psir_ylm, 
        int* vindex);

    //------------------------------------------------------
    // in gint_k_fvl.cpp 
    //------------------------------------------------------
    // set the orbital info 
    // set the derivative/Ylm information on each real space grid.
    void set_ijk_atom_force(
        const int &grid_index, 
        const int &size,
        double*** psir_ylm, 
        double*** dr, 
        bool** cal_flag, 
        double** distance, 
        double* ylma, 
        const double &delta_r,
        double*** dphi_x, 
        double ***dphi_y, 
        double*** dphi_z);

    // evaluate the force due to local potential.
    void evaluate_vl_force(
        const int &grid_index, 
        const int &size, 
        const int &i, 
        const int &j, 
        const int &k,
        double*** psir_ylm, 
        bool** cal_flag, 
        double* vldr3, 
        double** distance,
        double*** dphi_x, // gradient of orbital phi along x direction
        double*** dphi_y, // gradient of orbital phi along y direction
        double*** dphi_z, // gradient of orbital phi along z direction
        double* pvdpx, 
        double* pvdpy, 
        double* pvdpz,
        const Grid_Technique &gt);

    // evaluate the stresses due to local potential
    void evaluate_vl_stress(
        const int &grid_index, 
        const int &size, 
        const int &i, 
        const int &j, 
        const int &k,
        double*** psir_ylm, 
        bool** cal_flag, 
        double* vldr3, 
        double** distance,
        double*** dphi_x, 
        double*** dphi_y, 
        double*** dphi_z,
        double* pvdp11, 
        double* pvdp22, 
        double* pvdp33, 
        double* pvdp12, 
        double* pvdp13, 
        double* pvdp23,
        double*** dr,
        const Grid_Technique &gt);

    private:

    //----------------------------
    // key variable 
    //----------------------------
    // dimension: [GridT.lgd, GridT.nutot]
    // used only in vlocal with full H matrix.
    double* pvpR_pool;
    double** pvpR;

    double***** pvpR_tr; //LiuXh add 2019-07-15
    std::complex<double>***** pvpR_tr_soc; //LiuXh add 2019-07-15

    //----------------------------
    // key variable 
    //----------------------------
    // dimension: [GlobalC::LNNR.nnrg] 
    // save the < phi_0i | V | phi_Rj > in sparse H matrix.
    double** pvpR_reduced;

    //----------------------------
    // key variable 
    //----------------------------
    // dimension: [GridT.lgd, GridT.lgd]	
    // used only when folding the H matrix.
    std::complex<double>** pvp;
    std::complex<double>** pvp_nc[4];

    // used only in vlocal.
    int ik_now;
    int spin_now;

    // just pointer.
    bool pvpR_alloc_flag;
    bool reduced;

    double** DM_R; //pointer to Local_Orbital_Charge::DM_R
};

#endif
