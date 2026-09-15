#include "vdwd3_parameters.h"
#include "vdw_xcname.h"

#include <algorithm>
#include <cctype>
#include <cstddef>
#include <map>
#include <string>

namespace vdw
{
namespace d3
{
namespace
{

struct DampingParameterRecord
{
    const char* method_id;
    Damping damping;
    double s6;
    double s8;
    double s9;
    double rs6;
    double rs8;
    double a1;
    double a2;
    double alp;
};

struct MethodAlias
{
    const char* alias;
    const char* method_id;
};

#include "data/d3_damping_parameters.inc"
#include "data/d3_method_aliases.inc"

std::string lower_ascii(std::string value)
{
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char character) {
        return static_cast<char>(std::tolower(character));
    });
    return value;
}

std::string normalized(std::string value)
{
    value = lower_ascii(value);
    std::string result;
    result.reserve(value.size());
    for (unsigned char character : value)
    {
        if (std::isalnum(character) || character == '/')
        {
            result.push_back(static_cast<char>(character));
        }
    }
    return result;
}

std::string libxc_component(const std::string& component, const std::string& marker)
{
    const std::string lower = lower_ascii(component);
    const std::size_t position = lower.find(marker);
    if (position == std::string::npos)
    {
        return normalized(component);
    }
    return normalized(component.substr(position + marker.size()));
}

std::string canonicalize_libxc_pair(const std::string& input)
{
    const std::size_t separator = input.find(':');
    std::string first = input.substr(0, separator);
    std::string second = input.substr(separator + 1);
    const std::string first_lower = lower_ascii(first);
    if (first_lower.find("_c_") != std::string::npos)
    {
        std::swap(first, second);
    }

    const std::string exchange = libxc_component(first, "_x_");
    const std::string correlation = libxc_component(second, "_c_");
    if (exchange == correlation)
    {
        return exchange;
    }

    // Keys use the same separator-free spelling returned by
    // libxc_component(). Only mixed exchange/correlation names need an
    // explicit mapping; identical component names are handled above.
    static const std::map<std::string, std::string> special = {
        {"b88:p86", "bp86"},       {"b88:lyp", "blyp"},
        {"b88:pbe", "bpbe"},       {"pber:pbe", "revpbe"},
        {"rpbe:pbe", "rpbe"},      {"optx:lyp", "olyp"},
        {"mpw91:pw91", "mpwpw"},  {"ms2:regtpss", "ms2"},
        {"ms2h:regtpss", "ms2h"}, {"pbe:oppbe", "pbeop"},
        {"b88:opb88", "bop"},
    };
    const auto found = special.find(exchange + ":" + correlation);
    return found == special.end() ? std::string() : found->second;
}

} // namespace

std::string canonicalize_method_name(const std::string& input)
{
    const std::string xc_name = normalize_xc_name(input);
    std::string canonical;
    if (xc_name.find(':') != std::string::npos)
    {
        const std::string pair = canonicalize_libxc_pair(xc_name);
        if (!pair.empty())
        {
            canonical = pair;
        }
    }

    if (canonical.empty())
    {
        const std::size_t xc = xc_name.find("_xc_");
        canonical = xc == std::string::npos ? normalized(xc_name)
                                             : normalized(xc_name.substr(xc + 4));
    }

    // ABACUS calls its HSE06 implementation "HSE" in user input and
    // pseudopotential metadata. Keep that established spelling at the adapter
    // boundary while using s-dftd3's canonical method identifier internally.
    if (canonical == "hse")
    {
        canonical = "hse06";
    }
    return canonical;
}

bool lookup_parameters(const std::string& method,
                       Damping damping,
                       Parameters& parameters,
                       std::string& canonical_method)
{
    canonical_method = canonicalize_method_name(method);

    // Preserve the two historical LibXC spellings whose D3 parameter set is
    // represented by the `wb97x` method identifier in s-dftd3.  Keep the
    // canonical functional name unchanged for logging, and only remap the
    // lookup key for the damping variant where the legacy ABACUS mapping was
    // defined.
    std::string lookup_method = canonical_method;
    if (damping == Damping::Rational && lookup_method == "wb97xv")
    {
        lookup_method = "wb97x";
    }
    else if (damping == Damping::Zero && lookup_method == "wb97xd3")
    {
        lookup_method = "wb97x";
    }

    const char* method_id = nullptr;
    for (const MethodAlias& alias : kMethodAliases)
    {
        if (normalized(alias.alias) == lookup_method)
        {
            method_id = alias.method_id;
            break;
        }
    }
    if (method_id == nullptr)
    {
        return false;
    }

    for (const DampingParameterRecord& record : kDampingParameters)
    {
        if (record.damping == damping && std::string(record.method_id) == method_id)
        {
            parameters.damping = damping;
            parameters.s6 = record.s6;
            parameters.s8 = record.s8;
            parameters.s9 = record.s9;
            parameters.rs6 = record.rs6;
            parameters.rs8 = record.rs8;
            parameters.a1 = record.a1;
            parameters.a2 = record.a2;
            parameters.alp = record.alp;
            return true;
        }
    }
    return false;
}

} // namespace d3
} // namespace vdw
