#ifndef ISTATE_CHARGE_H
#define ISTATE_CHARGE_H
#include "module_basis/module_pw/pw_basis.h"
#include "module_cell/klist.h"
#include "module_elecstate/module_dm/density_matrix.h"
#include "module_hamilt_lcao/module_gint/gint.h"
#include "module_hamilt_lcao/module_gint/gint_gamma.h"
#include "module_hamilt_lcao/module_gint/gint_k.h"
#include "module_hamilt_pw/hamilt_pwdft/parallel_grid.h"
#include "module_psi/psi.h"

#include <module_base/complexmatrix.h>
#include <module_base/matrix.h>
#include <stdexcept>
#include <vector>

/**
 * @brief Manages the computation of the charge densities for different bands (band-decomposed charge densities).
 *
 * This class is responsible for initializing and managing the
 * charge state computation process, offering functionality to
 * calculate and plot the decomposed charge density for specified bands.
 */
class IState_Charge
{
  public:
    IState_Charge(psi::Psi<double>* psi_gamma_in, const Parallel_Orbitals* ParaV_in);
    IState_Charge(psi::Psi<std::complex<double>>* psi_k_in, const Parallel_Orbitals* ParaV_in);

    ~IState_Charge();

    // For gamma_only
    void begin(Gint_Gamma& gg,
               double** rho,
               const ModuleBase::matrix& wg,
               const std::vector<double>& ef_all_spin,
               const int rhopw_nrxx,
               const int rhopw_nplane,
               const int rhopw_startz_current,
               const int rhopw_nx,
               const int rhopw_ny,
               const int rhopw_nz,
               const int bigpw_bz,
               const int bigpw_nbz,
               const bool gamma_only_local,
               const int nbands_istate,
               const std::vector<int>& out_pchg,
               const int nbands,
               const double nelec,
               const int nspin,
               const int nlocal,
               const std::string& global_out_dir,
               std::ofstream& ofs_warning,
               const UnitCell* ucell_in,
               Grid_Driver* GridD_in,
               const K_Vectors& kv);

    // For multi-k
    void begin(Gint_k& gk,
               double** rho,
               std::complex<double>** rhog,
               const ModuleBase::matrix& wg,
               const std::vector<double>& ef_all_spin,
               const ModulePW::PW_Basis* rho_pw,
               const int rhopw_nrxx,
               const int rhopw_nplane,
               const int rhopw_startz_current,
               const int rhopw_nx,
               const int rhopw_ny,
               const int rhopw_nz,
               const int bigpw_bz,
               const int bigpw_nbz,
               const bool gamma_only_local,
               const int nbands_istate,
               const std::vector<int>& out_pchg,
               const int nbands,
               const double nelec,
               const int nspin,
               const int nlocal,
               const std::string& global_out_dir,
               std::ofstream& ofs_warning,
               UnitCell* ucell_in,
               Grid_Driver* GridD_in,
               const K_Vectors& kv,
               const bool if_separate_k,
               Parallel_Grid* Pgrid,
               const int ngmc);

  private:
    /**
     * @brief Set this->bands_picked_ according to the mode, and process an error if the mode is not recognized.
     *
     * @param nbands_istate INPUT parameter nbands_istate.
     * @param out_pchg INPUT parameter out_pchg, vector.
     * @param nbands INPUT parameter nbands.
     * @param nelec Total number of electrons.
     * @param mode Selected mode.
     * @param fermi_band Calculated Fermi band.
     */
    void select_bands(const int nbands_istate,
                      const std::vector<int>& out_pchg,
                      const int nbands,
                      const double nelec,
                      const int mode,
                      const int fermi_band);

#ifdef __MPI
    /**
     * @brief Calculates the density matrix for a given band.
     *
     * This method calculates the density matrix for a given band using the wave function coefficients.
     * It performs a matrix multiplication to produce the density matrix.
     *
     * @param ib Band index.
     * @param nspin Number of spin channels.
     * @param nelec Total number of electrons.
     * @param nlocal Number of local orbitals.
     * @param wg Weight matrix for bands and spins (k-points).
     * @param DM Density matrix to be calculated.
     * @param kv K-vectors.
     */
    void idmatrix(const int& ib,
                  const int nspin,
                  const double& nelec,
                  const int nlocal,
                  const ModuleBase::matrix& wg,
                  elecstate::DensityMatrix<double, double>& DM,
                  const K_Vectors& kv);

    // For multi-k
    void idmatrix(const int& ib,
                  const int nspin,
                  const double& nelec,
                  const int nlocal,
                  const ModuleBase::matrix& wg,
                  elecstate::DensityMatrix<std::complex<double>, double>& DM,
                  const K_Vectors& kv,
                  const bool if_separate_k);

#endif
    std::vector<int> bands_picked_;
    psi::Psi<double>* psi_gamma;
    psi::Psi<std::complex<double>>* psi_k;
    const Parallel_Orbitals* ParaV;
};
#endif
