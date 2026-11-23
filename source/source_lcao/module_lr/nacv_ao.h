#ifndef ABACUS_LR_NACV_AO_H
#define ABACUS_LR_NACV_AO_H

#include <vector>

namespace ModuleBase
{
class matrix;
}
namespace container
{
class Tensor;
}
namespace hamilt
{
template <typename T> class HContainer;
}
namespace module_dm
{
template <typename TK, typename TR> class DensityMatrix;
}

namespace LR
{
template <typename T> struct GradientInputs;

/**
 * @brief Sum the spin-resolved interstate AO densities for an auxiliary CIS pair.
 * @param[in] inputs Validated orbitals and complete spin windows.
 * @param[in] ground True for a reference-to-excited transition.
 * @param[in] open_shell Whether the local amplitudes have alpha and beta blocks.
 * @param[in] bra Packed bra amplitudes; unused for the reference.
 * @param[in] ket Packed ket amplitudes, up block followed by down for UKS.
 * @param[in] spin_channel Closed-shell LR adaptation, or zero for UKS.
 * @return Owned spin-summed AO coefficients for the orbital connection.
 */
container::Tensor nac_density_ao(
    const LR::GradientInputs<double> &inputs,
    bool ground,
    bool open_shell,
    const double *bra,
    const double *ket,
    int spin_channel);

/**
 * @brief Construct sparse real-space density blocks for existing force contractions.
 * @param[in] inputs Borrowed geometry, neighbor lists and parallel AO layout.
 * @param[in] matrix Local AO coefficient tensor in inputs.pmat's layout.
 * @param[in] symmetrize True to apply the existing symmetric-density convention.
 * @return Density matrix with owned values and borrowed parallel-layout context.
 * @pre The referenced context remains valid while the returned matrix is used.
 */
module_dm::DensityMatrix<double, double> sparse_density(
    const LR::GradientInputs<double> &inputs, const container::Tensor &matrix, bool symmetrize);

/**
 * @brief Contract the antisymmetric moving-AO connection for an auxiliary state pair.
 * @param[in] inputs Borrowed current geometry and validated parallel AO layout.
 * @param[in] density Local interstate AO coefficient tensor, including spin factors.
 * @param[in] geometry Borrowed overlap derivatives for this geometry.
 * @return Owned atom-major nat-by-3 connection in bohr^-1.
 * @details Contracts the antisymmetric density with the bra-center overlap
 * derivative; the symmetric metric-response term is handled separately.
 */
ModuleBase::matrix antisymmetric_ao_connection(
    const LR::GradientInputs<double> &inputs, const container::Tensor &density, const LR::ForceGeometry &geometry);
} // namespace LR

#endif
