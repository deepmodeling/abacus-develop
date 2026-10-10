#ifndef ABACUS_LR_STATE_OVERLAP_H
#define ABACUS_LR_STATE_OVERLAP_H
#include "source_base/vector3.h"
#include <vector>

// Cross-geometry AO/MO overlaps and reference projection, without tracking state.
// Distributed buffers are local column-major blocks; no matrices are gathered here.
namespace ModuleBase { class Matrix3; }
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

// Local column-major matrices on one BLACS grid: O = C_old^dagger S C_current.
// pc describes AO-by-band coefficients; ps describes S; po describes band-by-band O.
// Row/column MO indices refer to the old/current step respectively; S need not be Hermitian.
// All layouts share a BLACS grid and block size, including ranks with empty local blocks.
template <typename T>
std::vector<T> mo_overlap_dist(const std::vector<T>& old_coefficients,
    const std::vector<T>& ao, const std::vector<T>& current_coefficients,
    const Parallel_2D& pc, const Parallel_2D& ps, const Parallel_2D& po);
// X is stored as A = X^T (virtual rows, occupied columns). Return local
// A_projected = O_virt^dagger A_previous O_occ, without normalization.
// po describes the full MO overlap; px describes the electron-hole block.
// previous is the old-step local reference (possibly a JT mixture), not the current roots.
// Unit blocks align the occupied/virtual submatrix origins with px.
template <typename T>
std::vector<T> project_reference_dist(const std::vector<T>& previous,
    const std::vector<T>& overlap, const Parallel_2D& po,
    const Parallel_2D& px, int nocc, int nvirt);


// Integer lattice translations selected by the same cross-geometry cutoff.
// Retain these indices for HContainer storage and exp(+i k.R) Fourier phases.
std::vector<ModuleBase::Vector3<int>> cross_image_indices(
    const ModuleBase::Vector3<double>& displacement,
    const ModuleBase::Matrix3& lattice, double cutoff);

}
#endif
