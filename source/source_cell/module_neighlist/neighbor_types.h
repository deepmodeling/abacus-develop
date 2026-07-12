#ifndef NEIGHBOR_TYPES_H
#define NEIGHBOR_TYPES_H

#include <cstdint>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>

namespace ModuleNeighList
{

using GlobalAtomId = std::int64_t;
using LocalAtomIndex = std::int32_t;
using NeighborCount = std::int32_t;

inline int checked_int_size(const std::size_t value, const char* context)
{
    if (value > static_cast<std::size_t>(std::numeric_limits<int>::max()))
    {
        throw std::overflow_error(std::string(context) + " exceeds int range.");
    }
    return static_cast<int>(value);
}

inline LocalAtomIndex checked_local_atom_index(const std::size_t value, const char* context)
{
    if (value > static_cast<std::size_t>(std::numeric_limits<LocalAtomIndex>::max()))
    {
        throw std::overflow_error(std::string(context) + " exceeds local atom index range.");
    }
    return static_cast<LocalAtomIndex>(value);
}

inline std::size_t checked_size_product(const std::size_t lhs,
                                        const std::size_t rhs,
                                        const char* context)
{
    if (lhs != 0 && rhs > std::numeric_limits<std::size_t>::max() / lhs)
    {
        throw std::overflow_error(std::string(context) + " size product overflows.");
    }
    return lhs * rhs;
}

inline std::size_t checked_size_sum(const std::size_t lhs,
                                    const std::size_t rhs,
                                    const char* context)
{
    if (rhs > std::numeric_limits<std::size_t>::max() - lhs)
    {
        throw std::overflow_error(std::string(context) + " size sum overflows.");
    }
    return lhs + rhs;
}

} // namespace ModuleNeighList

#endif // NEIGHBOR_TYPES_H
