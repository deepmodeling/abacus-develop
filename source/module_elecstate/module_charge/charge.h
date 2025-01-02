#ifndef CHARGE_H
#define CHARGE_H

#include "module_base/complexmatrix.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/parallel_global.h"
#include "module_basis/module_pw/pw_basis.h"
#include "module_cell/module_symmetry/symmetry.h"
#include "module_elecstate/fp_energy.h"
#include "module_hamilt_pw/hamilt_pwdft/parallel_grid.h"

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
    const Parallel_Grid* pgrid = nullptr;
  private:
    //temporary
    double *_space_rho = nullptr, *_space_rho_save = nullptr;
    std::complex<double> *_space_rhog = nullptr, *_space_rhog_save = nullptr;
    double *_space_kin_r = nullptr, *_space_kin_r_save = nullptr;

  public:
    double **nhat = nullptr; //compensation charge for PAW
    double **nhat_save = nullptr; //compensation charge for PAW
                                 // wenfei 2023-09-05

    double *rho_core = nullptr;
    std::complex<double> *rhog_core = nullptr;

    int prenspin = 1;

    void set_rhopw(ModulePW::PW_Basis* rhopw_in);

    /**
     * @brief Init charge density from file or atomic pseudo-wave-functions
     *
     * @param eferm_iout [out] fermi energy to be initialized
     * @param ucell [in] unit cell
     * @param strucFac [in] structure factor
     * @param symm [in] symmetry
     * @param klist [in] k points list if needed
     * @param wfcpw [in] PW basis for wave function if needed
     */
    void init_rho(elecstate::efermi& eferm_iout,
                  const UnitCell& ucell,
                  const Parallel_Grid& pgrid,
                  const ModuleBase::ComplexMatrix& strucFac,
                  ModuleSymmetry::Symmetry& symm,
                  const void* klist = nullptr,
                  const void* wfcpw = nullptr);

    void allocate(const int &nspin_in);

    void atomic_rho(const int spin_number_need,
                    const double& omega,
                    double** rho_in,
                    const ModuleBase::ComplexMatrix& strucFac,
                    const UnitCell& ucell) const;

    void set_rho_core(const UnitCell& ucell,
                      const ModuleBase::ComplexMatrix& structure_factor, 
                      const bool* numeric);
    void set_rho_core_paw();

    void renormalize_rho();

    double sum_rho() const;

    void save_rho_before_sum_band();

	// for non-linear core correction
    void non_linear_core_correction
    (
        const bool &numeric,
        const double omega,
        const double tpiba2,
        const int mesh,
        const double *r,
        const double *rab,
        const double *rhoc,
        double *rhocg
    ) const;

	double cal_rho2ne(const double *rho_in) const;

  void check_rho(); // to check whether the charge density is normal

  void init_final_scf(); //LiuXh add 20180619

	public:
    /**
     * @brief init some arrays for mpi_inter_pools, rho_mpi
     * 
     */
    void init_chgmpi();

    /**
     * @brief Sum rho at different pools (k-point parallelism).
     *        Only used when GlobalV::KPAR > 1
     */
    void rho_mpi();

	  /**
	   * @brief 	Reduce among different pools 
     *          If NPROC_IN_POOLs are all the same, use GlobalV::INTER_POOL
     *          else, gather rho in a POOL, and then reduce among different POOLs
	   * 
	   * @param array_rho f(rho): an array [nrxx]
	   */
	  void reduce_diff_pools(double* array_rho) const;

    void set_omega(double* omega_in){this->omega_ = omega_in;};

    // mohan add 2021-02-20
    int nrxx=0; // number of r vectors in this processor
    int nxyz = 0; // total number of r vectors
    int ngmc=0; // number of g vectors in this processor
    int nspin=0; // number of spins
    ModulePW::PW_Basis* rhopw = nullptr;// When double_grid is used, rhopw = rhodpw (dense grid)
    bool cal_elf = false; // whether to calculate electron localization function (ELF)
  private:

    void destroy();    // free arrays  liuyu 2023-03-12

    double* omega_ = nullptr; // omega for non-linear core correction

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
