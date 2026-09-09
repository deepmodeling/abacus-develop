#ifndef GET_WF_PW_H
#define GET_WF_PW_H

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
 * @brief Write real-space norms and complex components of selected PW states.
 *
 * The caller owns the wavefunction and bases, which must outlive this object.
 * T follows the solver wavefunction precision; host grid processing uses double.
 * Scratch buffers are local to each begin() call.
 */
template <typename T, typename Device>
class Get_wf_pw
{
  public:
    /** @brief Bind the wavefunction, grids, and global band/spin configuration. */
    Get_wf_pw(const psi::Psi<T, Device>& psi,
              const ModulePW::PW_Basis_K& pw_wfc,
              const ModulePW::PW_Basis& pw_rho,
              const ModulePW::PW_Basis& pw_rhod,
              const int nspin,
              const int global_nbands);

    /**
     * @brief Write independently selected norm and Re/Im cubes for every k-point.
     *
     * A spinor has one combined norm and separate up/down complex components.
     * Re/Im output includes the Bloch phase; all fields scale as omega^(-1/2).
     */
    void begin(const UnitCell& ucell,
               const Parallel_Grid& pgrid,
               const K_Vectors& kv,
               const std::vector<int>& out_wfc_norm,
               const std::vector<int>& out_wfc_re_im,
               const std::string& global_out_dir) const;

  private:
    const psi::Psi<T, Device>& psi_;
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
    const std::complex<double>* transform_wfc(const T* coefficients, const int ik, const int component, Workspace* work) const;

    void write_norm(const int band,
                    const UnitCell& ucell,
                    const Parallel_Grid& pgrid,
                    const K_Vectors& kv,
                    const std::string& out_dir,
                    const BandParallelLayout& layout,
                    Workspace* work) const;
    void write_complex(const int band,
                       const UnitCell& ucell,
                       const Parallel_Grid& pgrid,
                       const K_Vectors& kv,
                       const std::string& out_dir,
                       const BandParallelLayout& layout,
                       Workspace* work) const;

    void calc_norm(const int spin_index, const double scale, Workspace* work) const;
    void calc_phase(const int ik, const K_Vectors& kv, std::vector<std::complex<double>>* phase) const;
    void calc_component(const std::vector<std::complex<double>>& component,
                        const std::vector<std::complex<double>>& phase,
                        const double scale,
                        std::vector<double>* real,
                        std::vector<double>* imag) const;

    void write_cube(const int band,
                    const int component,
                    const int k_number,
                    const std::string& part,
                    const UnitCell& ucell,
                    const Parallel_Grid& pgrid,
                    const std::string& out_dir,
                    const std::vector<double>& values) const;
};
} // namespace ModuleIO

#endif // GET_WF_PW_H
