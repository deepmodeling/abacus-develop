#include "vdw_xcname.h"

#include <algorithm>
#include <cctype>
#include <string>

namespace vdw
{
namespace
{

std::string trim_ascii(const std::string& value)
{
    const auto first = std::find_if_not(value.begin(), value.end(), [](unsigned char character) {
        return std::isspace(character) != 0;
    });
    if (first == value.end())
    {
        return std::string();
    }

    const auto last = std::find_if_not(value.rbegin(), value.rend(), [](unsigned char character) {
        return std::isspace(character) != 0;
    }).base();
    return std::string(first, last);
}

std::string normalize_component(std::string value)
{
    value = trim_ascii(value);
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char character) {
        return static_cast<char>(std::tolower(character));
    });

    if (value.compare(0, 3, "xc_") == 0)
    {
        value.erase(0, 3);
    }
    return value;
}

} // namespace

std::string normalize_xc_name(const std::string& input)
{
    const std::size_t plus = input.find('+');
    const std::size_t colon = input.find(':');
    const std::size_t separator = plus != std::string::npos ? plus : colon;

    if (separator == std::string::npos)
    {
        return normalize_component(input);
    }

    return normalize_component(input.substr(0, separator)) + ":"
           + normalize_component(input.substr(separator + 1));
}

} // namespace vdw
