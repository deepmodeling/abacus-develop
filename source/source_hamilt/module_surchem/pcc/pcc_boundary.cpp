#include "pcc_boundary.h"

#include <algorithm>
#include <cctype>
#include <stdexcept>

namespace ModuleSccs
{

Boundary parse_boundary(const std::string& value)
{
    std::string normalized = value;
    std::transform(normalized.begin(), normalized.end(), normalized.begin(), [](const char character) {
        return static_cast<char>(std::tolower(static_cast<unsigned char>(character)));
    });
    if (normalized == "periodic")
    {
        return Boundary::Periodic;
    }
    if (normalized == "pcc_0d" || normalized == "pcc-0d")
    {
        return Boundary::Pcc0d;
    }
    if (normalized == "pcc_2d" || normalized == "pcc-2d")
    {
        return Boundary::Pcc2d;
    }
    throw std::invalid_argument("unknown PCC boundary: " + value);
}

} // namespace ModuleSccs
