#ifndef LIBRPA_BZ_SAMPLING_H
#define LIBRPA_BZ_SAMPLING_H

#include <stdexcept>

namespace RpaLriDetail
{
struct LibRpaStoredQIndex
{
    int coulomb_irreducible_index;
    int representative_scf_index;
};

inline int librpa_stored_coulomb_q_count(const int stored_kpoint_count)
{
    if (stored_kpoint_count <= 0)
    {
        throw std::invalid_argument("LibRPA BZ output requires stored k points.");
    }
    return stored_kpoint_count;
}

inline LibRpaStoredQIndex librpa_stored_q_index(const int stored_kpoint_index,
                                                const int stored_kpoint_count)
{
    if (stored_kpoint_index < 0 || stored_kpoint_index >= stored_kpoint_count)
    {
        throw std::invalid_argument("LibRPA BZ output found an invalid stored k-point index.");
    }
    const int one_based_index = stored_kpoint_index + 1;
    return {one_based_index, one_based_index};
}
} // namespace RpaLriDetail

#endif
