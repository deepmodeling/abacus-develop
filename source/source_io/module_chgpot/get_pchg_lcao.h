#ifndef GET_PCHG_LCAO_H
#define GET_PCHG_LCAO_H

#include "source_base/parallel_grid.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_cell/klist.h"
#include "source_estate/module_dm/density_matrix.h"
#include "source_psi/psi.h"

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
    Get_pchg_lcao(const psi::Psi<double>& psi, const Parallel_Orbitals& para_orb, int nspin, double nelec);
    Get_pchg_lcao(const psi::Psi<std::complex<double>>& psi, const Parallel_Orbitals& para_orb, int nspin, double nelec);

    // For gamma_only
    void begin_gamma(double* const* rho,
                     const ModuleBase::matrix& wg,
                     const UnitCell& ucell,
                     const Parallel_Grid& pgrid,
                     const Grid_Driver& grid_driver,
                     const std::vector<int>& out_pchg,
                     const std::string& global_out_dir,
                     std::ofstream& ofs_running);

    // For multi-k
    void begin_k(double* const* rho,
                 std::complex<double>* const* rhog,
                 const ModuleBase::matrix& wg,
                 const ModulePW::PW_Basis& rho_pw,
                 UnitCell& ucell,
                 const Parallel_Grid& pgrid,
                 const Grid_Driver& grid_driver,
                 const K_Vectors& kv,
                 const std::vector<int>& out_pchg,
                 bool if_separate_k,
                 const std::string& global_out_dir,
                 std::ofstream& ofs_running);

  private:
    const psi::Psi<double>* const psi_gamma_;
    const psi::Psi<std::complex<double>>* const psi_k_;
    const Parallel_Orbitals& para_orb_;
    const int nspin_;
    const int nbands_;
    const int fermi_band_;

    void prepare_get_pchg(std::ofstream& ofs_running);

    /**
     * @brief Build the selected-band mask and reject invalid selectors.
     *
     * @param out_pchg INPUT parameter out_pchg, vector.
     */
    std::vector<int> select_bands(const std::vector<int>& out_pchg) const;

#ifdef __MPI
    /**
     * @brief Calculates the density matrix for a given band.
     *
     * This method calculates the density matrix for a given band using the wave function coefficients.
     * It performs a matrix multiplication to produce the density matrix.
     *
     * @param ib Band index.
     * @param wg Weight matrix for bands and spins (k-points).
     * @param DM Density matrix to be calculated.
     */
    void idmatrix(const int& ib, const ModuleBase::matrix& wg, elecstate::DensityMatrix<double, double>& DM);

    // For multi-k
    void idmatrix(const int& ib,
                  const ModuleBase::matrix& wg,
                  elecstate::DensityMatrix<std::complex<double>, double>& DM,
                  const K_Vectors& kv,
                  bool if_separate_k);

#endif
};
#endif // GET_PCHG_LCAO_H
