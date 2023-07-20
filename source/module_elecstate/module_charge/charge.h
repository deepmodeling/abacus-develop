#ifndef CHARGE_H
#define CHARGE_H

#include "module_base/complexmatrix.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/parallel_global.h"
#include "module_basis/module_pw/pw_basis.h"
#include "module_elecstate/fp_energy.h"

//a forward declaration of UnitCell
class UnitCell;

//==========================================================
// Electron Charge Density
//==========================================================
class Charge
{

  public:
    Charge();
    ~Charge();

    //==========================================================
    // MEMBER VARIABLES :
    // init_chg : "atomic" or "file"
    // NAME : total number of electrons
    // NAME : rho (nspin,ncxyz), the charge density in real space
    // NAME : rho_save (nspin,ncxyz), for charge mixing
    // NAME : rhog, charge density in G space
    // NAME : rhog_save, chage density in G space
    // NAME : rho_core [nrxx], the core charge in real space
    // NAME : rhog_core [ngm], the core charge in reciprocal space
    //==========================================================

    double **rho = nullptr;
    double **rho_save = nullptr;

    std::complex<double> **rhog = nullptr;
    std::complex<double> **rhog_save = nullptr;

    double **kin_r = nullptr; // kinetic energy density in real space, for meta-GGA
    double **kin_r_save = nullptr; // kinetic energy density in real space, for meta-GGA
                                   // wenfei 2021-07-28

    double *rho_core = nullptr;
    std::complex<double> *rhog_core = nullptr;

    int prenspin = 1;

    void set_rhopw(ModulePW::PW_Basis* rhopw_in);

    /**
     * @brief Init charge density from file or atomic pseudo-wave-functions
     * 
     * @param eferm_iout fermi energy to be initialized
     * @param strucFac [in] structure factor 
     * @param nbz [in] number of big grids in z direction
     * @param bz [in] number of small grids in big grids for z dirction
     */
    void init_rho(elecstate::efermi& eferm_iout, const ModuleBase::ComplexMatrix& strucFac, const int& nbz, const int& bz);
    
    void allocate(const int &nspin_in);

    void atomic_rho(const int spin_number_need,
                    const double& omega,
                    double** rho_in,
                    const ModuleBase::ComplexMatrix& strucFac,
                    const UnitCell& ucell) const;

    void set_rho_core(const ModuleBase::ComplexMatrix &structure_factor);

    void renormalize_rho(void);

    void save_rho_before_sum_band(void);

	// for non-linear core correction
    void non_linear_core_correction
    (
        const bool &numeric,
        const int mesh,
        const double *r,
        const double *rab,
        const double *rhoc,
        double *rhocg
    ) const;

	double check_ne(const double *rho_in) const;

    void init_final_scf(); //LiuXh add 20180619

	public:
    /**
     * @brief init some arrays for mpi_inter_pools, rho_mpi
     * 
     * @param nbz number of bigz in big grids
     * @param bz  number of z for each bigz
     */
    void init_chgmpi(const int& nbz, const int& bz);

    //sum rho of all pools
    void rho_mpi();

	  /**
	   * @brief 	Reduce among different pools 
     *          If NPROC_IN_POOLs are all the same, use GlobalV::INTER_POOL
     *          else, gather rho in a POOL, and then reduce among different POOLs
	   * 
	   * @param array_rho f(rho): an array [nrxx]
	   */
	  void reduce_diff_pools(double* array_rho) const;

    // mohan add 2021-02-20
    int nrxx; // number of r vectors in this processor
    int nxyz; // total nuber of r vectors
    int ngmc; // number of g vectors in this processor
    int nspin; // number of spins
    ModulePW::PW_Basis* rhopw = nullptr;
  private:
    double sum_rho(void) const;

    void destroy();    // free arrays  liuyu 2023-03-12

    bool allocate_rho;

    bool allocate_rho_final_scf; // LiuXh add 20180606
#ifdef __MPI
  private:
    bool use_intel_pool = false; //use INTER_POOL when NPROC_IN_POOLs are all the same
    int *rec = nullptr; //The number of elements each process should receive into the receive buffer.
    int *dis = nullptr; //The displacement (relative to recvbuf) for each process in the receive buffer.
#endif
    
};

#endif // charge
