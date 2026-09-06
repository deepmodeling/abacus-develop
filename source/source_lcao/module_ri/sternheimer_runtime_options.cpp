#include "source_lcao/module_ri/sternheimer_runtime_options.h"

#include "source_io/module_parameter/input_parameter.h"

#include <algorithm>
#include <cctype>
#include <cerrno>
#include <cstdlib>
#include <limits>
#include <stdexcept>
#include <string>

namespace ModuleRI
{

bool validate_sternheimer_molecular_coulomb(const Input_para& input)
{
    if (input.sternheimer_molecular_coulomb == "none")
    {
        if (!input.sternheimer_ao_potential_file.empty())
        {
            throw std::invalid_argument(
                "sternheimer_ao_potential_file requires sternheimer_molecular_coulomb isolated_ri.");
        }
        return false;
    }
    if (input.sternheimer_molecular_coulomb != "isolated_ri")
    {
        throw std::invalid_argument("sternheimer_molecular_coulomb must be none or isolated_ri.");
    }
    if (input.sternheimer_ao_potential_file.empty() || input.basis_type != "lcao" || !input.sternheimer_delta
        || !input.out_sternheimer_librpa || input.out_sternheimer_siab || input.sternheimer_q_index != 0
        || input.symmetry != "-1")
    {
        throw std::invalid_argument("isolated_ri requires an AO tensor, LCAO Delta-ST LibRPA output, q_index 0, "
                                    "symmetry -1 and no SIAB output. Gamma alone does not declare an isolated system.");
    }
    if (input.exx_singularity_correction != "limits" || input.out_librpa_2d_coulomb_method != "ewald"
        || input.out_librpa_3d_coulomb_method != "ewald")
    {
        throw std::invalid_argument("isolated_ri requires the complete free-space Center2/limits producer metric, "
                                    "not a periodic Ewald or direct-reciprocal metric.");
    }
    if (input.exx_rotate_abfs || input.exx_coul_moment || input.shrink_abfs_pca_thr >= 0.0 || input.cal_force)
    {
        throw std::invalid_argument(
            "isolated_ri does not yet support rotated, moment-split or shrunken ABFS, or forces; "
            "producer and response must use the same finalized radial basis.");
    }
    return true;
}

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
