#ifndef GET_PCHG_LCAO_H
#define GET_PCHG_LCAO_H

#include "source_cell/klist.h"
#include "source_estate/module_dm/density_matrix.h"

/**
 * @brief Manages the computation of the charge densities for different bands (band-decomposed charge densities).
 *
 * This class is responsible for initializing and managing the
 * charge state computation process, offering functionality to
 * calculate and plot the decomposed charge density for specified bands.
 */
class Get_pchg_lcao
{
  public:
    Get_pchg_lcao(psi::Psi<double>* psi_gamma_in, const Parallel_Orbitals* ParaV_in);
    Get_pchg_lcao(psi::Psi<std::complex<double>>* psi_k_in, const Parallel_Orbitals* ParaV_in);

    ~Get_pchg_lcao();

    // For gamma_only
    void begin(double** rho,
               const ModuleBase::matrix& wg,
               const std::vector<double>& ef_all_spin,
               const int rhopw_nrxx,
               const std::vector<int>& out_pchg,
               const int nbands,
               const double nelec,
               const int nspin,
               const UnitCell* ucell_in,
               const Parallel_Grid& pgrid,
               const Grid_Driver* GridD_in,
               const K_Vectors& kv,
               const std::string& global_out_dir,
               std::ofstream& ofs_running);

    // For multi-k
    void begin(double** rho,
               std::complex<double>** rhog,
               const ModuleBase::matrix& wg,
               const std::vector<double>& ef_all_spin,
               const ModulePW::PW_Basis* rho_pw,
               const int rhopw_nrxx,
               const std::vector<int>& out_pchg,
               const int nbands,
               const double nelec,
               const int nspin,
               UnitCell* ucell_in,
               const Parallel_Grid& pgrid,
               const Grid_Driver* GridD_in,
               const K_Vectors& kv,
               const std::string& global_out_dir,
               std::ofstream& ofs_running,
               const bool if_separate_k,
               const int chr_ngmc);

  private:
    void prepare_get_pchg(std::ofstream& ofs_running);

    /**
     * @brief Set this->bands_picked_ according to the mode, and process an error if the mode is not recognized.
     *
     * @param out_pchg INPUT parameter out_pchg, vector.
     * @param nbands INPUT parameter nbands.
     * @param fermi_band Calculated Fermi band.
     */
    void select_bands(const std::vector<int>& out_pchg, const int nbands, const int fermi_band);

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
     * @param wg Weight matrix for bands and spins (k-points).
     * @param DM Density matrix to be calculated.
     * @param kv K-vectors.
     */
    void idmatrix(const int& ib,
                  const int nspin,
                  const double& nelec,
                  const ModuleBase::matrix& wg,
                  elecstate::DensityMatrix<double, double>& DM,
                  const K_Vectors& kv);

    // For multi-k
    void idmatrix(const int& ib,
                  const int nspin,
                  const double& nelec,
                  const ModuleBase::matrix& wg,
                  elecstate::DensityMatrix<std::complex<double>, double>& DM,
                  const K_Vectors& kv,
                  const bool if_separate_k);

#endif
    std::vector<int> bands_picked_;
    psi::Psi<double>* psi_gamma = nullptr;
    psi::Psi<std::complex<double>>* psi_k = nullptr;
    const Parallel_Orbitals* ParaV;
};
#endif // GET_PCHG_LCAO_H
