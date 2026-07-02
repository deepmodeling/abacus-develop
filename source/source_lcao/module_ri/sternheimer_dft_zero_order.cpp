#include "source_lcao/module_ri/sternheimer_dft_zero_order.h"

#include <cmath>
#include <stdexcept>

namespace ModuleRI
{

SternheimerDFTZeroOrderComparison compare_sternheimer_fd_zero_order_to_dft(
    const SternheimerFDZeroOrderStates& fd_states,
    const std::vector<double>& dft_eigenvalues,
    const double eigenvalue_tolerance)
{
    if (eigenvalue_tolerance < 0.0)
    {
        throw std::invalid_argument("Sternheimer DFT zero-order comparison requires a non-negative tolerance.");
    }

    if (fd_states.eigenvalues.size() != fd_states.residual_norms.size())
    {
        throw std::invalid_argument("Sternheimer DFT zero-order comparison has inconsistent FD state sizes.");
    }

    if (fd_states.eigenvalues.size() != dft_eigenvalues.size())
    {
        throw std::invalid_argument("Sternheimer DFT zero-order comparison requires matching band counts.");
    }

    SternheimerDFTZeroOrderComparison comparison;
    comparison.bands.reserve(fd_states.eigenvalues.size());

    for (std::size_t ib = 0; ib != fd_states.eigenvalues.size(); ++ib)
    {
        SternheimerDFTZeroOrderBandComparison band;
        band.band_index = static_cast<int>(ib);
        band.fd_eigenvalue = fd_states.eigenvalues[ib];
        band.dft_eigenvalue = dft_eigenvalues[ib];
        band.fd_minus_dft = band.fd_eigenvalue - band.dft_eigenvalue;
        band.fd_residual_norm = fd_states.residual_norms[ib];
        band.eigenvalue_within_tolerance = std::abs(band.fd_minus_dft) <= eigenvalue_tolerance;

        comparison.max_abs_fd_minus_dft = std::max(comparison.max_abs_fd_minus_dft, std::abs(band.fd_minus_dft));
        comparison.all_eigenvalues_within_tolerance
            = comparison.all_eigenvalues_within_tolerance && band.eigenvalue_within_tolerance;
        comparison.bands.push_back(band);
    }

    return comparison;
}

} // namespace ModuleRI
