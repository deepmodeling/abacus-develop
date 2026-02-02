#ifndef SNAP_PSIBETA_GPU_H
#define SNAP_PSIBETA_GPU_H

#include "source_base/vector3.h"
#include "source_basis/module_ao/ORB_read.h"
#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_cell/module_neighbor/sltk_grid_driver.h"
#include "source_cell/setup_nonlocal.h"
#include "source_cell/unitcell.h"

#include <complex>
#include <unordered_map>
#include <vector>

#ifdef __CUDA
#include <cuda_runtime.h>
#endif

namespace module_rt
{
namespace gpu
{

/**
 * @brief Initialize GPU resources for snap_psibeta module (copy grids to constant memory)
 *        Should be called at the start of each calculate_HR
 */
void init_snap_psibeta_gpu();

/**
 * @brief Atom-level GPU batch processing interface
 *
 * Processes ALL neighbors for a center atom in a SINGLE kernel launch.
 * This significantly reduces kernel launch overhead compared to neighbor-level batching.
 *
 * @param orb Orbital information
 * @param infoNL_ Non-local pseudopotential information
 * @param T0 Center atom type (projector location)
 * @param R0 Center atom position (already multiplied by lat0)
 * @param A Vector potential
 * @param adjs Adjacent atom information for this center atom
 * @param ucell Unit cell pointer
 * @param paraV Parallel orbitals information
 * @param npol Polarization number
 * @param nlm_dim 1 for no current, 4 for current calculation
 * @param nlm_tot Output: nlm_tot[ad][dir][iw_index] = nlm_vector
 */
void snap_psibeta_atom_batch_gpu(
    const LCAO_Orbitals& orb,
    const InfoNonlocal& infoNL_,
    const int T0,
    const ModuleBase::Vector3<double>& R0,
    const ModuleBase::Vector3<double>& A,
    const AdjacentAtomInfo& adjs,
    const UnitCell* ucell,
    const Parallel_Orbitals* paraV,
    const int npol,
    const int nlm_dim,
    std::vector<std::vector<std::unordered_map<int, std::vector<std::complex<double>>>>>& nlm_tot);

} // namespace gpu
} // namespace module_rt

#endif // SNAP_PSIBETA_GPU_H
