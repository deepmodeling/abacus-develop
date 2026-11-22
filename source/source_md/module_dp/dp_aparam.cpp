#include "source_md/module_dp/dp_aparam.h"

#include "source_cell/module_neighlist/local_atom.h"

#include <stdexcept>

namespace ModuleMD
{
namespace DPDetail
{
std::vector<double> arrange_aparam(
    const std::vector<double> &input,
    const std::vector<LocalAtom> &owned,
    const std::vector<LocalAtom> &ghost,
    const std::vector<std::int64_t> &type_counts,
    int dimension,
    bool include_ghost)
{
    if (dimension < 0)
    {
        throw std::invalid_argument("DeePMD dim_aparam cannot be negative.");
    }
    const std::size_t width = static_cast<std::size_t>(dimension);
    const std::size_t offset_count = type_counts.size() + 1;
    std::vector<std::size_t> offsets(offset_count, 0);
    for (std::size_t type = 0; type < type_counts.size(); ++type)
    {
        if (type_counts[type] < 0)
        {
            throw std::invalid_argument("DeePMD type atom count cannot be negative.");
        }
        offsets[type + 1] = offsets[type] + static_cast<std::size_t>(type_counts[type]);
    }
    if (width == 0)
    {
        if (!input.empty())
        {
            throw std::invalid_argument("dp_aparam is set but the model has dim_aparam=0.");
        }
        return std::vector<double>();
    }
    if (input.size() != width && input.size() != offsets.back() * width)
    {
        throw std::invalid_argument("dp_aparam size must equal dim_aparam or nat * dim_aparam.");
    }

    const std::size_t count = owned.size() + (include_ghost ? ghost.size() : 0);
    const std::size_t output_size = count * width;
    std::vector<double> output(output_size, 0.0);
    for (std::size_t iat = 0; iat < count; ++iat)
    {
        const LocalAtom &atom = iat < owned.size() ? owned[iat] : ghost[iat - owned.size()];
        if (atom.type < 0 || static_cast<std::size_t>(atom.type) >= type_counts.size() || atom.type_index < 0 ||
            atom.type_index >= type_counts[atom.type])
        {
            throw std::invalid_argument("MDCell atomic parameter identity is invalid.");
        }
        const std::size_t global_index = offsets[atom.type] + static_cast<std::size_t>(atom.type_index);
        const std::size_t source = input.size() == width ? 0 : global_index * width;
        for (std::size_t component = 0; component < width; ++component)
        {
            output[iat * width + component] = input[source + component];
        }
    }
    return output;
}
} // namespace DPDetail
} // namespace ModuleMD
