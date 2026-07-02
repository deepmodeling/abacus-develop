#ifndef STERNHEIMER_DFT_ZERO_ORDER_H
#define STERNHEIMER_DFT_ZERO_ORDER_H

#include "source_lcao/module_ri/sternheimer_fd_solver.h"

#include <vector>

namespace ModuleRI
{

struct SternheimerDFTZeroOrderBandComparison
{
    int band_index = 0;
    double fd_eigenvalue = 0.0;
    double dft_eigenvalue = 0.0;
    double fd_minus_dft = 0.0;
    double fd_residual_norm = 0.0;
    bool eigenvalue_within_tolerance = false;
};

struct SternheimerDFTZeroOrderComparison
{
    std::vector<SternheimerDFTZeroOrderBandComparison> bands;
    double max_abs_fd_minus_dft = 0.0;
    bool all_eigenvalues_within_tolerance = true;
};

SternheimerDFTZeroOrderComparison compare_sternheimer_fd_zero_order_to_dft(
    const SternheimerFDZeroOrderStates& fd_states,
    const std::vector<double>& dft_eigenvalues,
    double eigenvalue_tolerance);

} // namespace ModuleRI

#endif
