#ifndef GET_PCHG_PW_H
#define GET_PCHG_PW_H

#include "source_base/parallel_grid.h"
#include "source_basis/module_pw/pw_basis_k.h"
#include "source_cell/klist.h"
#include "source_cell/unitcell.h"
#include "source_io/module_output/band_parallel_output.h"
#include "source_psi/psi.h"

#include <complex>
#include <string>
#include <vector>

namespace ModuleIO
{
/**
 * @brief Write band-resolved PW partial charges on the dense real-space grid.
 *
 * The caller owns the wavefunction and bases, which must outlive this object.
 * Scratch buffers are local to each begin() call.
 */
template <typename Device>
class Get_pchg_pw
{
  public:
    /** @brief Bind the wavefunction, grids, and global band/spin configuration. */
    Get_pchg_pw(const psi::Psi<std::complex<double>, Device>& psi,
                const ModulePW::PW_Basis_K& pw_wfc,
                const ModulePW::PW_Basis& pw_rho,
                const ModulePW::PW_Basis& pw_rhod,
                const int nspin,
                const int global_nbands);

    /**
     * @brief Write selected bands, either per k-point or after k summation and symmetry.
     * @param ucell Cell whose symmetry workspace may be updated during symmetrization.
     * @param noncolin Whether to retain transverse spinor magnetization components.
     */
    void begin(UnitCell* ucell,
               const Parallel_Grid& pgrid,
               const K_Vectors& kv,
               const std::vector<int>& out_pchg,
               const std::string& global_out_dir,
               const bool if_separate_k,
               const bool noncolin) const;

  private:
    const psi::Psi<std::complex<double>, Device>& psi_;
    const ModulePW::PW_Basis_K& pw_wfc_;
    const ModulePW::PW_Basis& pw_rho_;
    const ModulePW::PW_Basis& pw_rhod_;
    const int nspin_;
    const int global_nbands_;

    // Defined in the implementation to keep device buffers out of this header.
    class Workspace;

    // Validate the binary selector and return a zero-padded global band mask.
    std::vector<int> select_bands(const std::vector<int>& selection, const std::string& parameter_name) const;

    // Transform the owner's band and broadcast each local real-space slab.
    void transform_band(const int global_band, const int ik, const BandParallelLayout& layout, Workspace* work) const;
    // The returned data belongs to the selected workspace component and is valid
    // until that component is transformed again or the workspace is destroyed.
    const std::complex<double>* transform_wfc(const std::complex<double>* coefficients,
                                              const int ik,
                                              const int component,
                                              Workspace* work) const;

    void write_separate(const int band,
                        const UnitCell& ucell,
                        const Parallel_Grid& pgrid,
                        const K_Vectors& kv,
                        const std::string& out_dir,
                        const bool noncolin,
                        const BandParallelLayout& layout,
                        Workspace* work) const;
    void write_summed(const int band,
                      UnitCell* ucell,
                      const Parallel_Grid& pgrid,
                      const K_Vectors& kv,
                      const std::string& out_dir,
                      const bool noncolin,
                      const BandParallelLayout& layout,
                      Workspace* work) const;

    void calc_density(const int spin_index, const double weight, const bool noncolin, const bool accumulate, Workspace* work) const;
    void sum_pools(const Parallel_Grid& pgrid, const K_Vectors& kv, Workspace* work) const;
    void symmetrize(UnitCell* ucell, Workspace* work) const;

    void write_cube(const int band,
                    const int component,
                    const int k_number,
                    const UnitCell& ucell,
                    const Parallel_Grid& pgrid,
                    const std::string& out_dir,
                    const bool separate_k,
                    const std::vector<double>& values) const;
};
} // namespace ModuleIO

#endif // GET_PCHG_PW_H
