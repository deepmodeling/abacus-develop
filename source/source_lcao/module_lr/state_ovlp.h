#ifndef ABACUS_LR_STATE_OVERLAP_H
#define ABACUS_LR_STATE_OVERLAP_H
#include "source_base/vector3.h"
#include <vector>

// Cross-geometry AO/MO overlaps and reference projection, without tracking state.
// Distributed buffers are local column-major blocks; no matrices are gathered here.
namespace ModuleBase
{
class Matrix3;
}
class UnitCell;
class Parallel_Orbitals;
class Parallel_2D;
class TwoCenterIntegrator;
namespace LR
{
// Construct ordered old-bra/new-ket S(R), then fold at Gamma. Return this rank's
// column-major AO block; old_positions are Cartesian in cell.lat0 units, cutoffs in Bohr.
// Uses independent cross-geometry candidates and preserves nonsymmetry.
std::vector<double> cross_ao_overlap(const UnitCell& cell,
    const std::vector<ModuleBase::Vector3<double>>& old_positions,
    const std::vector<double>& cutoff, const TwoCenterIntegrator& integrator,
    const Parallel_Orbitals& pmat);

// Construct one occupied/virtual O = C_old^dagger S C_current block.
// first_band is the zero-based column offset into the full AO-by-band pc layout;
// po describes the selected square block. Coefficient columns are redistributed
// to origin (1,1) before multiplication; all layouts share grid and block size.
template <typename T>
std::vector<T> mo_overlap_dist(const std::vector<T>& old_coefficients,
    const std::vector<T>& ao, const std::vector<T>& current_coefficients,
    const Parallel_2D& pc, const Parallel_2D& ps, const Parallel_2D& po, int first_band);
// A = X^T (virtual rows, occupied columns): return O_vv^dagger A_previous O_oo.
// pi/pv describe already aligned occupied/virtual overlaps, reused for all roots.
// This function performs no redistribution or normalization.
template <typename T>
std::vector<T> project_reference_dist(const std::vector<T>& previous,
    const std::vector<T>& occupied_overlap, const std::vector<T>& virtual_overlap,
    const Parallel_2D& pi, const Parallel_2D& pv, const Parallel_2D& px);

// Integer lattice translations selected by the same cross-geometry cutoff.
// Retain these indices for HContainer storage and exp(+i k.R) Fourier phases.
std::vector<ModuleBase::Vector3<int>> cross_image_indices(
    const ModuleBase::Vector3<double>& displacement,
    const ModuleBase::Matrix3& lattice, double cutoff);

}
#endif
