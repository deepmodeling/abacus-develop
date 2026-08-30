#include "source_lcao/module_ri/sternheimer_runtime_options.h"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <stdexcept>
#include <string>

namespace ModuleRI
{

bool sternheimer_environment_flag(const char* name, const bool default_value)
{
    if (name == nullptr || name[0] == '\0')
    {
        throw std::invalid_argument("Sternheimer environment flag name must not be empty.");
    }

    const char* raw = std::getenv(name);
    if (raw == nullptr)
    {
        return default_value;
    }

    std::string value(raw);
    std::transform(value.begin(), value.end(), value.begin(), [](const unsigned char character) {
        return static_cast<char>(std::tolower(character));
    });
    if (value == "1" || value == "true" || value == "on" || value == "yes")
    {
        return true;
    }
    if (value == "0" || value == "false" || value == "off" || value == "no")
    {
        return false;
    }

    throw std::invalid_argument("Invalid boolean value for " + std::string(name) + ": " + value);
}

} // namespace ModuleRI
