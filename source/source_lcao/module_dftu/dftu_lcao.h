#ifndef DFTU_H
#define DFTU_H

#include "source_cell/klist.h"
#include "source_cell/unitcell.h"
#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_estate/module_charge/charge_mixing.h"
#include "source_pw/module_pwdft/dftu_base.h"
#ifdef __LCAO
#include "source_basis/module_ao/orb_read.h"
#include "source_hamilt/hamilt.h"
#include "source_hamilt/module_hcontainer/hcontainer.h"
#include "source_estate/module_dm/density_matrix.h"
#include "source_lcao/force_stress_arrays.h" // mohan add 2024-06-15
#endif

#include <string>
#include <vector>


class Plus_U : public Plus_U_Base
{

  public:
    Plus_U();
    ~Plus_U();

  public:
    // allocate relevant data strcutures
    void init(UnitCell& cell,
                const Parallel_Orbitals* pv,
                const int npol,
                const int nspin,
                const std::vector<int>& orbital_corr,
                const bool yukawa_potential,
                const double yukawa_lambda,
                const std::string& global_readin_dir,
                const std::string& global_out_dir,
                const std::string& init_chg,
                const int nlocal,
                const bool gamma_only_local,
                const std::string& ks_solver,
                const bool cal_force,
                const bool cal_stress,
                const std::string& device,
                const int kpar,
                const std::vector<double>& hubbard_u,
                const double uramping,
                const int occ_mat_ctrl,
                const int mixing_dftu
#ifdef __LCAO
                , const LCAO_Orbitals* orb = nullptr
#endif
                );

    // calculate the energy correction
    void cal_energy_correction(const UnitCell& ucell, const int istep);

  private:

    const Parallel_Orbitals* paraV = nullptr;

    double yukawa_lambda = 0.0;
    int npol = 1;
    int nlocal = 0;
    bool gamma_only_local = false;
    std::string ks_solver;
    bool cal_force = false;
    bool cal_stress = false;

#ifdef __LCAO
    const LCAO_Orbitals* ptr_orb_ = nullptr;
    std::vector<double> orb_cutoff_;
#endif

#ifdef __LCAO
    //=============================================================
    // In dftu_hamilt.cpp
    // For calculating contribution to Hamiltonian matrices
    //=============================================================
  public:
	void cal_eff_pot_mat_complex(const int ik,
			std::complex<double>* eff_pot,
			const std::vector<int>& isk,
			const std::complex<double>* sk,
			const int npol);

	void cal_eff_pot_mat_real(const int ik,
			double* eff_pot,
			const std::vector<int>& isk,
			const double* sk,
			const int npol);

    void cal_eff_pot_mat_R_double(const int ispin, double* SR, double* HR, const int npol);

	void cal_eff_pot_mat_R_complex_double(const int ispin,
			std::complex<double>* SR,
			std::complex<double>* HR,
			const int npol);
#endif

#ifdef __LCAO
    // calculate the local occupation number matrix
    void cal_occup_m_k(const int iter,
                       const UnitCell& ucell,
                       const std::vector<std::vector<std::complex<double>>>& dm_k,
                       const K_Vectors& kv,
                       const double& mixing_beta,
                       hamilt::Hamilt<std::complex<double>>* p_ham);

    void cal_occup_m_gamma(const int iter,
                           const UnitCell& ucell,
                           const std::vector<std::vector<double>>& dm_gamma,
                           const double& mixing_beta,
                           hamilt::Hamilt<double>* p_ham);
#endif

#ifdef __LCAO
private:
    //=============================================================
    // In dftu_tools.cpp
    // For calculating onsite potential, which is used
    // for both Hamiltonian and force/stress
    //=============================================================

    void cal_VU_pot_mat_complex(const int spin, const bool newlocale, std::complex<double>* VU, const int npol);
    void cal_VU_pot_mat_real(const int spin, const bool newlocale, double* VU, const int npol);

    double get_onebody_eff_pot(const int T,
                               const int iat,
                               const int L,
                               const int N,
                               const int spin,
                               const int m0,
                               const int m1,
                               const bool newlocale);

    //=============================================================
    // In dftu_folding.cpp
    // Subroutines for folding S and dS matrix
    //=============================================================

    void fold_dSR_gamma(const UnitCell& ucell,
                        const Parallel_Orbitals& pv,
                        const Grid_Driver* gd,
                        double* dsloc_x,
                        double* dsloc_y,
                        double* dsloc_z,
                        double* dh_r,
                        const int dim1,
                        const int dim2,
                        double* dSR_gamma);

    // dim = 0 : S, for Hamiltonian
    // dim = 1-3 : dS, for force
    // dim = 4-6 : dS * dR, for stress

    void folding_matrix_k(const UnitCell& ucell,
                          const Grid_Driver& gd,
                          ForceStressArrays& fsr,
                          const Parallel_Orbitals& pv,
                          const int ik,
                          const int dim1,
                          const int dim2,
                          std::complex<double>* mat_k,
                          const ModuleBase::Vector3<double>& kvec_d);

    /**
     * @brief new function of folding_S_matrix
     * only for Hamiltonian now, for force and stress will be developed later
     * use HContainer as input and output in mat_k
    */
	void folding_matrix_k_new(const int ik,
			hamilt::Hamilt<std::complex<double>>* p_ham);

    //=============================================================
    // In dftu_force.cpp
    // For calculating force and stress fomr DFT+U
    //=============================================================
 public:
   void force_stress(const UnitCell& ucell,
                     const Grid_Driver& gd,
					 std::vector<std::vector<double>>* dmk_d,
					 std::vector<std::vector<std::complex<double>>>* dmk_c,
					 const Parallel_Orbitals& pv,
                     ForceStressArrays& fsr,
                     ModuleBase::matrix& force_dftu,
                     ModuleBase::matrix& stress_dftu,
                     const K_Vectors& kv,
                     const int npol);

 private:
   void cal_force_k(const UnitCell& ucell,
                    const Grid_Driver& gd,
                    ForceStressArrays& fsr,
                    const Parallel_Orbitals& pv,
                    const int ik,
                    const std::complex<double>* rho_VU,
                    ModuleBase::matrix& force_dftu,
                    const ModuleBase::Vector3<double>& kvec_d);

   void cal_stress_k(const UnitCell& ucell,
                     const Grid_Driver& gd,
                     ForceStressArrays& fsr,
                     const Parallel_Orbitals& pv,
                     const int ik,
                     const std::complex<double>* rho_VU,
                     ModuleBase::matrix& stress_dftu,
                     const ModuleBase::Vector3<double>& kvec_d);

   void cal_force_gamma(const UnitCell& ucell,
                        const double* rho_VU,
                        const Parallel_Orbitals& pv,
                        double* dsloc_x,
                        double* dsloc_y,
                        double* dsloc_z,
                        ModuleBase::matrix& force_dftu);

   void cal_stress_gamma(const UnitCell& ucell,
                         const Parallel_Orbitals& pv,
                         const Grid_Driver* gd,
                         double* dsloc_x,
                         double* dsloc_y,
                         double* dsloc_z,
                         double* dh_r,
                         const double* rho_VU,
                         ModuleBase::matrix& stress_dftu);
#endif

    //=============================================================
    // In dftu_yukawa.cpp
    // Relevant for calculating U using Yukawa potential
    //=============================================================

  public:
    void cal_slater_UJ(const UnitCell& ucell, double** rho, const int& nrxx);

  private:
    void cal_slater_Fk(const UnitCell& ucell,const int L, const int T); // L:angular momnet, T:atom type
    void cal_yukawa_lambda(double** rho, const int& nrxx);

    double spherical_Bessel(const int k, const double r, const double lambda);
    double spherical_Hankel(const int k, const double r, const double lambda);

#ifdef __LCAO
  public:
    /**
     * @brief get the density matrix of target spin
     * nspin = 1 and 4 : ispin should be 0
     * nspin = 2 : ispin should be 0/1
    */
    const hamilt::HContainer<double>* get_dmr(int ispin) const;
    /**
     * @brief set the density matrix for DFT+U calculation
     * if the density matrix is not set or set to nullptr, the DFT+U calculation will not be performed
    */
    void set_dmr(const elecstate::DensityMatrix<double, double>* dm_in_dftu_d);
    void set_dmr(const elecstate::DensityMatrix<std::complex<double>, double>* dm_in_dftu_cd);

  private:
    const UnitCell* ucell = nullptr;
    const elecstate::DensityMatrix<double, double>* dm_in_dftu_d = nullptr;
    const elecstate::DensityMatrix<std::complex<double>, double>* dm_in_dftu_cd = nullptr;
#endif
};


#ifdef __LCAO
template <typename T>
void dftu_cal_occup_m(const int iter,
                      const UnitCell& ucell,
                      const std::vector<std::vector<T>>& dm,
                      const K_Vectors& kv,
                      const double& mixing_beta,
                      hamilt::Hamilt<T>* p_ham,
                      Plus_U &dftu);
#endif


#endif
