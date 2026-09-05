#include "vdwd3_data.h"

#include <cassert>
#include <cstddef>

namespace vdw
{
namespace d3
{
namespace data
{
namespace
{

#include "data/d3_reference.inc"

inline std::size_t pair_index(int high, int low)
{
    return static_cast<std::size_t>(high * (high - 1) / 2 + low - 1);
}

inline void check_atomic_number(int atomic_number)
{
    assert(atomic_number > 0 && atomic_number <= max_element);
}

} // namespace

int reference_count(int atomic_number)
{
    check_atomic_number(atomic_number);
    return static_cast<int>(kReferenceCount[atomic_number]);
}

double reference_cn(int atomic_number, int reference)
{
    check_atomic_number(atomic_number);
    assert(reference >= 0 && reference < reference_count(atomic_number));
    return kReferenceCn[atomic_number * max_reference + reference];
}

double reference_c6(int atomic_number_i, int reference_i, int atomic_number_j, int reference_j)
{
    check_atomic_number(atomic_number_i);
    check_atomic_number(atomic_number_j);
    assert(reference_i >= 0 && reference_i < reference_count(atomic_number_i));
    assert(reference_j >= 0 && reference_j < reference_count(atomic_number_j));

    int high = atomic_number_i;
    int low = atomic_number_j;
    int high_reference = reference_i;
    int low_reference = reference_j;
    if (high < low)
    {
        const int atomic_number = high;
        high = low;
        low = atomic_number;
        const int reference = high_reference;
        high_reference = low_reference;
        low_reference = reference;
    }

    const std::size_t pair = pair_index(high, low);
    const std::size_t offset = kReferenceC6Offset[pair]
                               + high_reference * reference_count(low) + low_reference;
    return kReferenceC6[offset];
}

double covalent_radius(int atomic_number)
{
    check_atomic_number(atomic_number);
    return kCovalentRadius[atomic_number];
}

double r4r2(int atomic_number)
{
    check_atomic_number(atomic_number);
    return kR4R2[atomic_number];
}

double vdw_radius(int atomic_number_i, int atomic_number_j)
{
    check_atomic_number(atomic_number_i);
    check_atomic_number(atomic_number_j);
    const int high = atomic_number_i > atomic_number_j ? atomic_number_i : atomic_number_j;
    const int low = atomic_number_i > atomic_number_j ? atomic_number_j : atomic_number_i;
    return kVdwRadius[pair_index(high, low)];
}

} // namespace data
} // namespace d3
} // namespace vdw
