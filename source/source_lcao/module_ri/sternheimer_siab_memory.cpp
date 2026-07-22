#include "source_lcao/module_ri/sternheimer_siab_memory.h"

#include "source_lcao/module_ri/sternheimer_abfs_perturbation.h"

#include <complex>
#include <cstddef>
#include <limits>
#include <stdexcept>

namespace ModuleRI
{

SternheimerSIABMemoryEstimate estimate_sternheimer_siab_dense_memory(const int grid_size,
                                                                     const int raw_auxiliary_channels,
                                                                     const int response_channels,
                                                                     const int primitive_count,
                                                                     const int reciprocal_count,
                                                                     const int occupied_states,
                                                                     const int frequency_count)
{
    if (grid_size <= 0 || raw_auxiliary_channels <= 0 || response_channels <= 0 || primitive_count <= 0
        || reciprocal_count <= 0 || occupied_states <= 0 || frequency_count <= 0)
    {
        throw std::invalid_argument("Sternheimer SIAB memory estimate requires positive dimensions.");
    }
    const auto multiply = [](const std::uint64_t left, const std::uint64_t right) {
        if (right != 0 && left > std::numeric_limits<std::uint64_t>::max() / right)
        {
            throw std::overflow_error("Sternheimer SIAB memory estimate overflowed uint64_t.");
        }
        return left * right;
    };
    const auto add = [](const std::uint64_t left, const std::uint64_t right) {
        if (left > std::numeric_limits<std::uint64_t>::max() - right)
        {
            throw std::overflow_error("Sternheimer SIAB memory estimate overflowed uint64_t.");
        }
        return left + right;
    };

    SternheimerSIABMemoryEstimate result;
    result.coulomb_metric_bytes = multiply(multiply(raw_auxiliary_channels, raw_auxiliary_channels), sizeof(double));
    result.transformed_potential_bytes = multiply(multiply(grid_size, response_channels), sizeof(double));
    const std::uint64_t transform_chunk = static_cast<std::uint64_t>(sternheimer_abfs_transform_grid_chunk);
    const std::uint64_t transform_columns
        = add(static_cast<std::uint64_t>(raw_auxiliary_channels),
              static_cast<std::uint64_t>(response_channels));
    result.channel_transform_workspace_bytes
        = add(multiply(multiply(transform_chunk, transform_columns), sizeof(double)),
              multiply(transform_chunk, sizeof(std::size_t)));
    result.reciprocal_primitive_bytes
        = multiply(multiply(primitive_count, reciprocal_count), sizeof(std::complex<double>));
    result.primitive_overlap_bytes
        = multiply(multiply(primitive_count, primitive_count), sizeof(std::complex<double>));
    const std::uint64_t row_bytes
        = add(7U * sizeof(double), multiply(primitive_count, sizeof(std::complex<double>)));
    const std::uint64_t row_count
        = multiply(multiply(occupied_states, response_channels), frequency_count);
    result.gathered_reference_row_bytes = multiply(3U, multiply(row_count, row_bytes));
    result.total_bytes = add(add(add(result.coulomb_metric_bytes, result.transformed_potential_bytes),
                                 result.channel_transform_workspace_bytes),
                             add(result.reciprocal_primitive_bytes,
                                 add(result.primitive_overlap_bytes, result.gathered_reference_row_bytes)));
    return result;
}

} // namespace ModuleRI
