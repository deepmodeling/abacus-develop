#ifndef SNAP_PSIBETA_GPU_H
#define SNAP_PSIBETA_GPU_H

#include "source_base/vector3.h"
#include "source_basis/module_ao/ORB_read.h"
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
 * @brief Initialize GPU resources (copy grids to constant memory)
 *        Should be called at the start of each calculate_HR
 */
void initialize_gpu_resources();

/**
 * @brief Release GPU resources (clear any error states)
 *        Should be called at the end of each calculate_HR
 */
void finalize_gpu_resources();

/**
 * @brief Neighbor-level GPU batch processing interface
 *
 * Computes matrix elements for all orbitals on a neighbor atom with projectors
 * Replaces the original orbital loop:
 *   for (int iw1l = 0; iw1l < all_indexes.size(); iw1l += npol) {
 *       snap_psibeta_half_tddft(...);
 *   }
 *
 * @param orb Orbital information
 * @param infoNL_ Non-local pseudopotential information
 * @param T1 Neighbor atom type
 * @param R1 Neighbor atom position (already multiplied by lat0)
 * @param atom1 Neighbor atom pointer (for iw2l, iw2m, iw2n)
 * @param all_indexes All orbital indices on this neighbor
 * @param npol Polarization number
 * @param T0 Projector atom type
 * @param R0 Projector atom position (already multiplied by lat0)
 * @param A Vector potential
 * @param nlm_all Output: nlm_all[dir][iw_index] = nlm_vector
 * @param calc_r Whether to calculate position operator matrix elements
 * @return true if GPU processing succeeded, false if caller should use CPU path
 */
bool snap_psibeta_neighbor_batch_gpu(const LCAO_Orbitals& orb,
                                     const InfoNonlocal& infoNL_,
                                     const int T1,
                                     const ModuleBase::Vector3<double>& R1,
                                     const Atom* atom1,
                                     const std::vector<int>& all_indexes,
                                     const int npol,
                                     const int T0,
                                     const ModuleBase::Vector3<double>& R0,
                                     const ModuleBase::Vector3<double>& A,
                                     std::vector<std::unordered_map<int, std::vector<std::complex<double>>>>& nlm_all,
                                     const bool calc_r);

} // namespace gpu
} // namespace module_rt

#endif // SNAP_PSIBETA_GPU_H
