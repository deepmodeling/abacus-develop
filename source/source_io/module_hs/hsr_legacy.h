#ifndef HSR_LEGACY_H
#define HSR_LEGACY_H

#include "source_base/matrix.h"
#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_basis/module_nao/two_center_bundle.h"
#include "source_cell/klist.h"
#include "source_cell/module_neighbor/sltk_grid_driver.h"
#include "source_hamilt/hamilt.h"
#include "source_lcao/lcao_hs_arrays.h"

#include <string>

// Legacy LCAO_HS_Arrays-based sparse matrix output path (dH/dR, dS/dR, T(R), S(R)).
// Kept as-is from the former write_hs_r.h; new code should prefer the
// HContainer-based writers in hsr_writer.h.

namespace ModuleIO
{
// Groups the shared format / path / runtime flags threaded through the
// LCAO_HS_Arrays-based writers. Replaces ten positional parameters that were
// identical across output_dHR/output_dSR/output_TR.
struct MatROutputOptions
{
    bool binary = false;
    double sparse_threshold = 0.0;
    int precision = 8;
    std::string global_out_dir;
    std::string global_matrix_dir;
    std::string calculation;
    bool out_app_flag = false;
    int nspin = 1;
};

void output_dHR(const int& istep,
                const ModuleBase::matrix& v_eff,
                const UnitCell& ucell,
                const Parallel_Orbitals& pv,
                LCAO_HS_Arrays& HS_Arrays,
                const Grid_Driver& grid, // mohan add 2024-04-06
                const TwoCenterBundle& two_center_bundle,
                const LCAO_Orbitals& orb,
                const MatROutputOptions& options,
                const bool gamma_only_local,
                const int npol,
                const int nlocal);

void output_dSR(const int& istep,
                const UnitCell& ucell,
                const Parallel_Orbitals& pv,
                LCAO_HS_Arrays& HS_Arrays,
                const Grid_Driver& grid, // mohan add 2024-04-06
                const TwoCenterBundle& two_center_bundle,
                const LCAO_Orbitals& orb,
                const MatROutputOptions& options,
                const bool gamma_only_local,
                const int npol,
                const int nlocal);

void output_TR(const int istep,
               const UnitCell& ucell,
               const Parallel_Orbitals& pv,
               LCAO_HS_Arrays& HS_Arrays,
               const Grid_Driver& grid,
               const TwoCenterBundle& two_center_bundle,
               const LCAO_Orbitals& orb,
               const std::string& TR_filename,
               const MatROutputOptions& options);

template <typename TK>
void output_SR(Parallel_Orbitals& pv,
               const Grid_Driver& grid,
               hamilt::Hamilt<TK>* p_ham,
               const std::string& SR_filename,
               const bool& binary,
               const double& sparse_threshold,
               const int precision,
               const std::string& global_out_dir,
               const std::string& global_matrix_dir,
               const std::string& calculation,
               const bool out_app_flag,
               const int nspin);

} // namespace ModuleIO

#endif
