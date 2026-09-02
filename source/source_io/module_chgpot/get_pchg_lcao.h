#ifndef GET_PCHG_LCAO_H
#define GET_PCHG_LCAO_H

#include "source_base/parallel_grid.h"
#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_cell/klist.h"
#include "source_cell/module_neighbor/sltk_grid_driver.h"
#include "source_cell/unitcell.h"
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
    Get_pchg_lcao(const psi::Psi<double>& psi, const Parallel_Orbitals& para_orb, int nspin);
    Get_pchg_lcao(const psi::Psi<std::complex<double>>& psi, const Parallel_Orbitals& para_orb, int nspin);

    // For gamma_only
    void begin_gamma(const UnitCell& ucell,
                     const Parallel_Grid& pgrid,
                     const Grid_Driver& grid_driver,
                     const std::vector<int>& out_pchg,
                     const std::string& global_out_dir,
                     std::ofstream& ofs_running);

    // For multi-k
    void begin_k(const ModulePW::PW_Basis& rho_pw,
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

    void prepare_get_pchg(std::ofstream& ofs_running);

    /**
     * @brief Build the selected-band mask and reject invalid selectors.
     *
     * @param out_pchg INPUT parameter out_pchg, vector.
     */
    std::vector<int> select_bands(const std::vector<int>& out_pchg) const;
};
#endif // GET_PCHG_LCAO_H
