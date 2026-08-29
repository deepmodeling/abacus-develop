#ifndef DFTU_LCAO_H
#define DFTU_LCAO_H

#include "source_cell/klist.h"
#include "source_cell/unitcell.h"
#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_pw/module_pwdft/dftu_base.h"
#ifdef __LCAO
#include "source_basis/module_ao/orb_read.h"
#include "source_hamilt/hamilt.h"
#include "source_hamilt/module_hcontainer/hcontainer.h"
#include "source_estate/module_dm/density_matrix.h"
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
    void cal_eff_pot_mat_R_double(const int ispin, double* SR, double* HR, const int npol);

	void cal_eff_pot_mat_R_complex_double(const int ispin,
			std::complex<double>* SR,
			std::complex<double>* HR,
			const int npol);
#endif

#ifdef __LCAO
    // calculate the local occupation number matrix
    void cal_occ_mat_k(const int iter,
                       const UnitCell& ucell,
                       const std::vector<std::vector<std::complex<double>>>& dm_k,
                       const K_Vectors& kv,
                       const double& mixing_beta,
                       hamilt::Hamilt<std::complex<double>>* p_ham);

    void cal_occ_mat_gamma(const int iter,
                           const UnitCell& ucell,
                           const std::vector<std::vector<double>>& dm_gamma,
                           const double& mixing_beta,
                           hamilt::Hamilt<double>* p_ham);
#endif

#ifdef __LCAO
    //=============================================================
    // In dftu_tools.cpp
    // For calculating onsite potential, which is used
    // for both Hamiltonian and force/stress
    //=============================================================
  public:
    void pot_onsite_complex(const int spin, const bool newlocale, std::complex<double>* pot_onsite, const int npol);
    void pot_onsite_real(const int spin, const bool newlocale, double* pot_onsite, const int npol);

  private:
    double get_onebody_eff_pot(const int T,
                               const int iat,
                               const int L,
                               const int N,
                               const int spin,
                               const int m0,
                               const int m1,
                               const bool newlocale);

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

    /// read-only accessors for state needed by DFTU_LCAO free functions
    const Parallel_Orbitals* get_paraV() const { return paraV; }
    int get_npol() const { return npol; }
    int get_nlocal() const { return nlocal; }
    const std::string& get_ks_solver() const { return ks_solver; }
    const std::vector<double>& get_orb_cutoff() const { return orb_cutoff_; }
    bool is_gamma_only_local() const { return gamma_only_local; }
    bool is_cal_force() const { return cal_force; }
    bool is_cal_stress() const { return cal_stress; }

  private:
    const UnitCell* ucell = nullptr;
    const elecstate::DensityMatrix<double, double>* dm_in_dftu_d = nullptr;
    const elecstate::DensityMatrix<std::complex<double>, double>* dm_in_dftu_cd = nullptr;
#endif
};


#endif
