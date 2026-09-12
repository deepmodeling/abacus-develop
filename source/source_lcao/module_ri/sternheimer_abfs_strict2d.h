#ifndef STERNHEIMER_ABFS_STRICT2D_H
#define STERNHEIMER_ABFS_STRICT2D_H

#include "source_lcao/module_ri/sternheimer_abfs_perturbation.h"

#include <cstddef>
#include <utility>

namespace ModuleRI
{

// Discrete strict-2D Coulomb contract (atomic units, distances in bohr):
// The open-z kernel follows from integrating the full Coulomb Fourier symbol:
// integral dkz/(2*pi) 4*pi*exp[i*kz*(z-z')]/(Q^2+kz^2)
//   = (2*pi/Q)*exp[-Q*abs(z-z')], Q>0.
//   rho_G(z_j) = (1/Nxy) sum_xy exp[-i(G+q).r_xy] rho_q(r_xy,z_j)
//   phi_q(r_xy,z_i) = sum_G exp[i(G+q).r_xy] (2*pi/Q)
//                     * dz * sum_j exp[-Q*dz*abs(i-j)] rho_G(z_j),
//   Q = |G+q|, dz = |a3|/nz, j = 0,...,nz-1.
// potential_r contains FULL complex Bloch densities on entry and potentials
// in Ha on return, not cell-periodic envelopes. Metadata is preserved and
// max_abs is updated. The planar 1/area Fourier normalization becomes 1/Nxy;
// there is no further cell-area factor and no Ry factor of two.
//
// All representable planar FFT modes are retained, with signed index
// i < (n+1)/2 ? i : i-n (negative even-grid Nyquist), exactly as in the ABFS
// periodic solver. q is used as supplied, without folding reciprocal shifts.
// The original fine grid has no duplicated endpoint: every z sample has the
// SAME dz weight, including both boundary samples (no half-endpoint weights).
// Each Q uses the open-z linear convolution, not a z FFT, minimum-image
// distance, cut Coulomb, or a sum of periodic-z images. Complexity per channel
// is O(nz*Nxy*log(Nxy) + nz*Nxy), with one reusable full-grid FFT workspace.
//
// Supported: grid.periodic=true for the in-plane Bloch sampling convention,
// nonsingular possibly skew a1,a2 in Cartesian xy, a3 parallel to Cartesian z,
// exactly qz=0, and non-Gamma planar q. Geometry roundoff tolerance is 1e-12
// relative to each lattice-vector norm. Gamma modulo planar reciprocal
// integers (tolerance 1e-14), any zero Q, tilted cells, nonfinite data and
// incompatible grid/channel sizes are explicitly rejected, without fallback.
// An explicit lattice determines dz and the metric; hx,hy,hz define the cell
// only when lattice_vectors is entirely zero, as in the existing grid helpers.
// Arithmetic overflow/nonfinite intermediates, FFT output or final magnitudes
// throw std::overflow_error. The failing channel and its max_abs are unchanged;
// earlier completed channels remain potentials (no whole-batch rollback).
//
// The kernel does NOT change ABF sampling or remove any periodic images
// already present in the input density. Physical admission therefore needs
// an unwrapped slab window with each compact ABF's z support contained inside
// it, or independent evidence that boundary density/image contamination is
// negligible. An atom near a z boundary can violate this prerequisite even
// though the discrete kernel is correct. Choose a suitable slab window in
// the caller; this API neither re-centers atoms nor modifies their densities.
//
// Exactness means the stated FINITE-GRID sum up to floating-point roundoff,
// not automatically exact continuum full-2D-Ewald matrix elements. Compare
// selected projected integrals below to same-q, same-ABFS full Ewald data and
// converge grid/window/support errors before admitting fineweak physics.
void solve_sternheimer_abf_strict2d_coulomb_in_place(
    std::vector<SternheimerABFBlochGridChannel>& density_channels,
    const SternheimerFDHamiltonian::Grid& grid,
    const SternheimerReducedKPoint& qpoint);

// Cheap selected matrix elements, in the supplied pair order (duplicates OK):
//   result[k] = dV * sum_r conj(densities[pairs[k].first].potential_r[r])
//                                * potentials[pairs[k].second].potential_r[r],
//   dV = |(a1 x a2).a3|/(nx*ny*nz).
// O(pairs.size()*grid.size()) work, no full auxiliary matrix or channel copy.
// Retain the original density samples before the in-place solve. Indices refer
// to the supplied vectors, not channel_index metadata; both fields must use
// the SAME q, grid, ABFS normalization and channel ordering. No symmetrization,
// eigenvalue clipping, Ewald rescaling, or physical pass/fail is performed.
// Only requested channels are inspected; unselected channel payloads may be
// absent. Density and potential vectors need not have the same channel count.
// Nonfinite products, running sums or volume-scaled results from finite inputs
// throw std::overflow_error; no partial result is returned.
// The caller compares these complex elements with the matching full Ewald
// entries; agreement/PSD alone does not establish response convergence.
std::vector<std::complex<double>> sternheimer_abf_strict2d_selected_coulomb_integrals(
    const std::vector<SternheimerABFBlochGridChannel>& densities,
    const std::vector<SternheimerABFBlochGridChannel>& potentials,
    const SternheimerFDHamiltonian::Grid& grid,
    const std::vector<std::pair<std::size_t, std::size_t>>& pairs);

} // namespace ModuleRI

#endif
