#ifndef SYMM_ROT_SPIN_H
#define SYMM_ROT_SPIN_H

#include "source_base/matrix3.h"

#include <array>
#include <complex>
#include <vector>

namespace ModuleSymmetry
{
/// @brief  SU(2) spin-1/2 representation of a real-space symmetry operation.
///
/// These utilities provide the spin (SU(2)) part of the symmetry operation needed for
/// nspin=4 (non-collinear / SOC) symmetrization. The orbital (real-spherical-harmonics)
/// part is handled separately by the existing Wigner-D / T(V) machinery; the full
/// representation of a symmetry operation on a spinor orbital is the tensor product
/// T(V) (x) U(V).
///
/// Convention notes:
///  - ABACUS uses the ROW-vector convention, i.e. a point transforms as r' = r * V,
///    so the cartesian rotation matrix `gmatc` stored by ABACUS equals the transpose of
///    the textbook (column-vector) rotation:  V_row = V_col^T.
///  - Spin is a pseudovector: it only sees the PROPER part of the operation. For an
///    improper operation (det(gmatc) = -1) we factor out the inversion and build U from
///    the proper rotation  R_proper = -gmatc  (which then has det = +1).
///  - With U built from the axis-angle of R_proper via
///        U = cos(theta/2) I - i sin(theta/2) (n . sigma),
///    the induced action on the Pauli (spin-density) vector is
///        rho'^i = sum_j W_ij rho^j,   with   W = R_proper^T  ( = R_col ),
///    where W_ij = (1/2) Tr(sigma_i U sigma_j U^dagger). This relation is the basis of
///    the unit tests and is convention-self-consistent regardless of the textbook
///    index ordering quoted in the formula document.
namespace SpinRotation
{
/// A 2x2 complex matrix stored row-major: {m00, m01, m10, m11}.
using Su2 = std::array<std::complex<double>, 4>;

/// @brief Build the SU(2) spin-1/2 matrix U corresponding to a cartesian symmetry
///        operation `gmatc` (row-vector convention, det = +-1).
///
/// Handles the special cases theta = 0 (U = I) and theta = pi (axis from the symmetric
/// part of the rotation, where the standard axis-angle formula is singular). For an
/// improper operation the proper part R_proper = -gmatc is used.
///
/// The returned U is defined up to the double-group sign (+-U); both signs give the same
/// similarity transform U D U^dagger, so this ambiguity is harmless for symmetrization.
Su2 so3_to_su2(const ModuleBase::Matrix3& gmatc, const double eps = 1e-6);

/// @brief The proper rotation acting on the spin-density 3-vector (Pauli x,y,z
///        components), i.e. the W matrix with rho'^i = sum_j W_ij rho^j.
///        Computed directly from gmatc as W = R_proper^T (R_proper = proper part).
ModuleBase::Matrix3 spin_so3(const ModuleBase::Matrix3& gmatc);

/// @brief The same W matrix computed independently from a given SU(2) matrix U via
///        W_ij = (1/2) Tr(sigma_i U sigma_j U^dagger). Used for verification.
ModuleBase::Matrix3 pauli_rotation_matrix(const Su2& U);

/// @brief Hermitian conjugate of a 2x2 SU(2) matrix.
Su2 dagger(const Su2& U);

/// @brief 2x2 complex matrix product A * B (both row-major).
Su2 mat2_mul(const Su2& A, const Su2& B);

/// @brief Rotate the 2x2 spin block of a spinor matrix element in place:
///        m_block' = U * m_block * U^dagger.
///        `block` is the 2x2 spin sub-matrix {uu, ud, du, dd} of a fixed orbital pair.
Su2 rotate_spin_block(const Su2& block, const Su2& U);

/// @brief Rotate the four Pauli components (rho^0, rho^x, rho^y, rho^z) of a single
///        real-space density point: rho^0 is unchanged, (rho^x, rho^y, rho^z) are mixed
///        by W = spin_so3(gmatc). The input/output are component values at one grid point.
void rotate_pauli_components(const ModuleBase::Matrix3& gmatc,
                             const double in[4],
                             double out[4]);
} // namespace SpinRotation
} // namespace ModuleSymmetry

#endif // SYMM_ROT_SPIN_H
