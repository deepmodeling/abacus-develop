#ifndef HYBRID_IMAGE_DOMAIN_H
#define HYBRID_IMAGE_DOMAIN_H

#include <RI/global/Tensor.h>
#include <RI/global/Global_Func-2.h>
#include <array>
#include <algorithm>
#include <complex>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <vector>

namespace RI_2D_Comm
{
using ImageCell = std::array<int, 3>;

inline ImageCell centered_image(const ImageCell& cell, const ImageCell& period)
{
    ImageCell result;
    for (int d = 0; d < 3; ++d)
    {
        if (period[d] <= 0)
        {
            throw std::invalid_argument("Nonpositive RI image period");
        }
        const long long p = period[d];
        const long long c = cell[d];
        result[d] = static_cast<int>((c % p + 3 * p / 2) % p - p / 2);
    }
    return result;
}

// Use the existing LibRI max-norm threshold, with no additional approximation.
template <typename T>
ImageCell retained_image_extent(std::map<int, std::map<std::pair<int, ImageCell>, RI::Tensor<T>>>& blocks,
                                const double threshold)
{
    ImageCell extent{{0, 0, 0}};
    for (auto& row : blocks)
    {
        for (auto block = row.second.begin(); block != row.second.end();)
        {
            const double norm = block->second.norm(std::numeric_limits<double>::max());
            if (!std::isfinite(norm))
            {
                throw std::runtime_error("Nonfinite hybrid RI integral");
            }
            if (norm <= threshold)
            {
                block = row.second.erase(block);
                continue;
            }
            for (int d = 0; d < 3; ++d)
            {
                const long long value = block->first.second[d];
                const long long magnitude = value < 0 ? -value : value;
                if (magnitude > std::numeric_limits<int>::max() / 8)
                {
                    throw std::overflow_error("Hybrid RI integral displacement too large");
                }
                extent[d] = std::max(extent[d], static_cast<int>(magnitude));
            }
            ++block;
        }
    }
    return extent;
}

inline ImageCell hybrid_image_period(const ImageCell& c_extent, const ImageCell& v_extent)
{
    ImageCell period;
    for (int d = 0; d < 3; ++d)
    {
        const long long bound = static_cast<long long>(v_extent[d]) + 2LL * c_extent[d];
        // LibRI's centered modulo internally forms 3*period.
        const long long max_bound = (std::numeric_limits<int>::max() / 3LL - 1) / 2;
        if (c_extent[d] < 0 || v_extent[d] < 0 || bound > max_bound)
        {
            throw std::overflow_error("Hybrid RI image domain too large");
        }
        period[d] = static_cast<int>(2 * bound + 1);
    }
    return period;
}

template <typename T>
void unfold_hybrid_images(std::vector<std::map<int, std::map<std::pair<int, ImageCell>, RI::Tensor<T>>>>& density,
                          const ImageCell& physical_period, const ImageCell& image_period,
                          const std::array<double, 3>& mesh_origin,
                          const std::array<double, 3>& input_phase_per_cell)
{
    // Convert each raw source to the ordinary NAO frame BEFORE history mixing.
    // Fuse exp(-i k0.B) and exp(-i A_source.R) into one tensor multiplication.
    // input_phase_per_cell[d] is A_source.(lat0 * lattice_row[d]), in radians.
    // Pass zeros when only unfolding an ordinary density without a field phase.
    // On a regular full mesh, D(R+B) = exp(-i 2*pi*k0.B) D(R).
    // Any mesh point can serve as k0; B is a physical BvK translation.
    for (auto& spin : density)
    {
        for (auto& row : spin)
        {
            std::map<int, std::map<ImageCell, RI::Tensor<T>>> by_atom;
            for (const auto& block : row.second)
            {
                by_atom[block.first.first].emplace(block.first.second, block.second);
            }
            std::map<std::pair<int, ImageCell>, RI::Tensor<T>> unfolded;
            for (const auto& atom : by_atom)
            {
                for (int x = -image_period[0] / 2; x <= image_period[0] / 2; ++x)
                {
                    for (int y = -image_period[1] / 2; y <= image_period[1] / 2; ++y)
                    {
                        for (int z = -image_period[2] / 2; z <= image_period[2] / 2; ++z)
                        {
                            const ImageCell cell{{x, y, z}};
                            const ImageCell representative = centered_image(cell, physical_period);
                            const auto source = atom.second.find(representative);
                            if (source != atom.second.end())
                            {
                                double cycles = 0.0;
                                double input_phase = 0.0;
                                for (int d = 0; d < 3; ++d)
                                {
                                    cycles += mesh_origin[d] * (cell[d] - representative[d]);
                                    input_phase += input_phase_per_cell[d] * cell[d];
                                }
                                const double angle = -2.0 * std::acos(-1.0) * std::remainder(cycles, 1.0) - input_phase;
                                const std::complex<double> twist = std::exp(std::complex<double>(0.0, angle));
                                const T factor = RI::Global_Func::convert<T>(twist);
                                const auto image = source->second * factor;
                                unfolded.emplace(std::make_pair(atom.first, cell), image);
                            }
                        }
                    }
                }
            }
            row.second.swap(unfolded);
        }
    }
}

template <typename T>
void decode_hybrid_images(std::map<int, std::map<std::pair<int, ImageCell>, RI::Tensor<T>>>& blocks,
                          const ImageCell& period)
{
    for (auto& row : blocks)
    {
        std::map<std::pair<int, ImageCell>, RI::Tensor<T>> decoded;
        for (const auto& block : row.second)
        {
            const ImageCell cell = centered_image(block.first.second, period);
            const auto inserted = decoded.emplace(std::make_pair(block.first.first, cell), block.second);
            if (!inserted.second)
            {
                throw std::runtime_error("Aliased output in hybrid RI image domain");
            }
        }
        row.second.swap(decoded);
    }
}
} // namespace RI_2D_Comm
#endif
