#include "rpa_lri_bz_io.h"

#include "rpa_lri_detail.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>

namespace RpaLriDetail
{
namespace
{
bool finite_point(const std::array<double, 3>& point)
{
    return std::isfinite(point[0]) && std::isfinite(point[1]) && std::isfinite(point[2]);
}

std::size_t full_grid_size(const std::array<int, 3>& mesh)
{
    unsigned long long size = 1;
    for (const int value: mesh)
    {
        if (value <= 0)
        {
            throw std::runtime_error("bz_sample requires a positive Monkhorst-Pack mesh.");
        }
        size = checked_mul_u64(size, static_cast<unsigned long long>(value), "bz_sample grid size");
    }
    if (size > static_cast<unsigned long long>(std::numeric_limits<std::size_t>::max()))
    {
        throw std::runtime_error("bz_sample grid size exceeds size_t range.");
    }
    return static_cast<std::size_t>(size);
}

void write_point(std::ofstream& output, const BzSamplePoint& point, const int index, const double weight)
{
    output << std::setw(8) << index << std::setw(24) << std::scientific << std::setprecision(15) << weight
           << std::setw(24) << std::scientific << std::setprecision(15) << point.direct[0] << std::setw(24)
           << std::scientific << std::setprecision(15) << point.direct[1] << std::setw(24) << std::scientific
           << std::setprecision(15) << point.direct[2] << std::setw(24) << std::scientific << std::setprecision(15)
           << point.cartesian[0] << std::setw(24) << std::scientific << std::setprecision(15) << point.cartesian[1]
           << std::setw(24) << std::scientific << std::setprecision(15) << point.cartesian[2] << std::setw(8) << index
           << std::setw(8) << index << '\n';
}
} // namespace

void write_bz_sample(const std::string& filename, const BzSampleData& data)
{
    const std::size_t full_count = full_grid_size(data.mesh);
    const int nk_full = checked_i32_from_size(full_count, "bz_sample full-grid size");
    const int nk_reduced = checked_i32_from_size(data.reduced.size(), "bz_sample reduced-grid size");
    if (data.reduced.empty() || data.reduced.size() > full_count)
    {
        throw std::runtime_error("bz_sample has an invalid reduced-grid size.");
    }
    if (data.full_cartesian.size() != full_count || data.full_map.size() != full_count)
    {
        throw std::runtime_error("bz_sample full-grid metadata has an invalid size.");
    }

    double weight_sum = 0.0;
    for (const BzSamplePoint& point: data.reduced)
    {
        if (!std::isfinite(point.weight) || point.weight <= 0.0 || !finite_point(point.direct)
            || !finite_point(point.cartesian))
        {
            throw std::runtime_error("bz_sample contains a non-finite or non-positive point.");
        }
        weight_sum += point.weight;
    }
    if (!std::isfinite(weight_sum) || weight_sum <= 0.0)
    {
        throw std::runtime_error("bz_sample has a non-positive total k-point weight.");
    }

    std::vector<int> multiplicity(data.reduced.size(), 0);
    if (full_count > data.reduced.size())
    {
        for (std::size_t ik = 0; ik < full_count; ++ik)
        {
            const int representative = data.full_map[ik];
            if (representative < 0 || representative >= nk_reduced)
            {
                throw std::runtime_error("bz_sample representative is outside the reduced grid.");
            }
            if (!finite_point(data.full_cartesian[ik]))
            {
                throw std::runtime_error("bz_sample contains a non-finite full-grid point.");
            }
            ++multiplicity[representative];
        }
        for (std::size_t ik = 0; ik < data.reduced.size(); ++ik)
        {
            const double expected = static_cast<double>(multiplicity[ik]) / static_cast<double>(full_count);
            if (multiplicity[ik] == 0 || std::abs(data.reduced[ik].weight / weight_sum - expected) > 1e-6)
            {
                throw std::runtime_error("bz_sample weight disagrees with full-grid multiplicity.");
            }
        }
    }

    std::ofstream output(filename.c_str(), std::ios::out | std::ios::trunc);
    if (!output.good())
    {
        throw std::runtime_error("Failed to open " + filename);
    }
    output << data.mesh[0] << std::setw(6) << data.mesh[1] << std::setw(6) << data.mesh[2] << '\n';
    output << nk_reduced << std::setw(8) << nk_reduced << '\n';
    for (std::size_t ik = 0; ik < data.reduced.size(); ++ik)
    {
        write_point(output, data.reduced[ik], static_cast<int>(ik + 1), data.reduced[ik].weight / weight_sum);
    }
    if (full_count > data.reduced.size())
    {
        output << "full_kmap " << nk_full << '\n';
        for (std::size_t ik = 0; ik < full_count; ++ik)
        {
            const auto& point = data.full_cartesian[ik];
            output << std::setw(8) << ik + 1 << std::setw(24) << point[0] << std::setw(24) << point[1] << std::setw(24)
                   << point[2] << std::setw(8) << data.full_map[ik] + 1 << '\n';
        }
    }
    output.close();
    if (!output.good())
    {
        throw std::runtime_error("Failed to write " + filename);
    }
}
} // namespace RpaLriDetail
