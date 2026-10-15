#include "rpa_lri_coeff_v1.h"

#include "rpa_lri_detail.h"

#include <RI/global/Tensor.h>
#include <fstream>
#include <ios>
#include <stdexcept>

namespace RpaLriDetail
{
void write_lri_coef_v1(const std::string& filename,
                       const std::int32_t natom,
                       std::vector<LriCoefV1Block>& blocks)
{
    if (natom <= 0)
    {
        throw std::runtime_error("LibRPA v1 coefficient output requires a positive atom count.");
    }

    const std::int64_t nblocks
        = checked_i64_from_u64(static_cast<unsigned long long>(blocks.size()), "LibRPA v1 coefficient block count");
    const unsigned long long header_bytes = 3 * sizeof(std::int32_t) + 2 * sizeof(std::int64_t);
    const unsigned long long record_bytes
        = 5 * sizeof(std::int32_t) + sizeof(double) + sizeof(std::int64_t);
    unsigned long long offset = checked_add_u64(
        header_bytes,
        checked_mul_u64(static_cast<unsigned long long>(nblocks), record_bytes, "LibRPA v1 coefficient record table"),
        "LibRPA v1 coefficient data offset");

    for (auto& block: blocks)
    {
        if (block.tensor == nullptr || block.nw1 <= 0 || block.nw2 <= 0 || block.naux <= 0
            || block.tensor->shape.size() != 3 || block.tensor->shape[0] != static_cast<std::size_t>(block.naux)
            || block.tensor->shape[1] != static_cast<std::size_t>(block.nw1)
            || block.tensor->shape[2] != static_cast<std::size_t>(block.nw2))
        {
            throw std::runtime_error("LibRPA v1 coefficient output has an invalid block shape.");
        }
        block.offset = checked_i64_from_u64(offset, "LibRPA v1 coefficient data offset");
        const unsigned long long values = checked_mul_u64(
            checked_mul_u64(static_cast<unsigned long long>(block.nw1),
                            static_cast<unsigned long long>(block.nw2),
                            "LibRPA v1 coefficient block size"),
            static_cast<unsigned long long>(block.naux),
            "LibRPA v1 coefficient block size");
        offset = checked_add_u64(
            offset,
            checked_mul_u64(values, sizeof(double), "LibRPA v1 coefficient block size"),
            "LibRPA v1 coefficient data size");
    }

    std::ofstream ofs(filename.c_str(), std::ios::out | std::ios::binary | std::ios::trunc);
    if (!ofs.good())
    {
        throw std::runtime_error("Failed to open " + filename);
    }

    const std::int32_t marker = LIBRPA_LRICOEF_V1_MARKER;
    const std::int32_t ncell = 0;
    write_scalar(ofs, marker, filename);
    write_scalar(ofs, natom, filename);
    write_scalar(ofs, ncell, filename);
    write_scalar(ofs, nblocks, filename);
    write_scalar(ofs, nblocks, filename);
    for (const auto& block: blocks)
    {
        write_scalar(ofs, block.ia1, filename);
        write_scalar(ofs, block.ia2, filename);
        for (const std::int32_t coord: block.r)
        {
            write_scalar(ofs, coord, filename);
        }
        write_scalar(ofs, block.max_abs, filename);
        write_scalar(ofs, block.offset, filename);
    }
    for (const auto& block: blocks)
    {
        for (int i = 0; i < block.nw1; ++i)
        {
            for (int j = 0; j < block.nw2; ++j)
            {
                for (int mu = 0; mu < block.naux; ++mu)
                {
                    const double value = real_as_double((*block.tensor)(mu, i, j));
                    write_scalar(ofs, value, filename);
                }
            }
        }
    }
    ofs.close();
    if (!ofs.good())
    {
        throw std::runtime_error("Failed to finish writing " + filename);
    }
}
} // namespace RpaLriDetail
