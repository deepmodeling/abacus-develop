#ifndef STERNHEIMER_SIAB_MEMORY_H
#define STERNHEIMER_SIAB_MEMORY_H

#include <cstdint>

namespace ModuleRI
{

struct SternheimerSIABMemoryEstimate
{
    std::uint64_t coulomb_metric_bytes = 0;
    std::uint64_t transformed_potential_bytes = 0;
    std::uint64_t channel_transform_workspace_bytes = 0;
    std::uint64_t reciprocal_primitive_bytes = 0;
    std::uint64_t primitive_overlap_bytes = 0;
    std::uint64_t gathered_reference_row_bytes = 0;
    std::uint64_t total_bytes = 0;
};

SternheimerSIABMemoryEstimate estimate_sternheimer_siab_dense_memory(int grid_size,
                                                                     int raw_auxiliary_channels,
                                                                     int response_channels,
                                                                     int primitive_count,
                                                                     int reciprocal_count,
                                                                     int occupied_states,
                                                                     int frequency_count);

} // namespace ModuleRI

#endif
