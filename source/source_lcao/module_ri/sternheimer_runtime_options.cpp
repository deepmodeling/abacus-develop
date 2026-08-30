#include "source_lcao/module_ri/sternheimer_runtime_options.h"

#include <algorithm>
#include <cctype>
#include <cerrno>
#include <cstdlib>
#include <limits>
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

int sternheimer_channel_batch_width()
{
    constexpr const char* name = "ABACUS_STERNHEIMER_CHANNEL_BATCH_WIDTH";
    constexpr int default_width = 2;
    constexpr int maximum_width = 64;
    const char* raw = std::getenv(name);
    if (raw == nullptr)
    {
        return default_width;
    }

    errno = 0;
    char* end = nullptr;
    const long parsed = std::strtol(raw, &end, 10);
    if (errno != 0 || end == raw || end == nullptr || *end != '\0' || parsed <= 0 || parsed > maximum_width
        || parsed > std::numeric_limits<int>::max())
    {
        throw std::invalid_argument("Invalid channel batch width in " + std::string(name) + ": " + raw);
    }
    return static_cast<int>(parsed);
}

} // namespace ModuleRI
