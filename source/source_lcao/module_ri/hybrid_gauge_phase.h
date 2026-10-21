#ifndef HYBRID_GAUGE_PHASE_H
#define HYBRID_GAUGE_PHASE_H

#include "source_base/vector3.h"
#include "source_base/matrix3.h"
#include "source_lcao/module_ri/mix_dmk_2d.h"

#include <RI/global/Tensor.h>

#include <array>
#include <algorithm>
#include <complex>
#include <cmath>
#include <map>
#include <utility>
#include <vector>
#include <stdexcept>
#include <tuple>

namespace RI_2D_Comm
{
// The RI Fourier convention is exp(-ik.R). Its sHG input therefore needs
// exp(-i A_source.R) to enter the ordinary NAO frame. A_source must identify
// the density's time layer, not merely the current Hamiltonian's stage.
inline std::complex<double> hybrid_input_phase(
    const ModuleBase::Vector3<double>& vector_potential,
    const ModuleBase::Vector3<double>& lattice_translation)
{
    const double argument = -(vector_potential * lattice_translation);
    return std::exp(std::complex<double>(0.0, argument));
}

enum class HybridDensityTime { Left, Right };

// Restarting SCF mixing does not change the time layer of the current density.
inline HybridDensityTime hybrid_density_time(const int step, const int iteration)
{
    return step > 0 && iteration == 1 ? HybridDensityTime::Left : HybridDensityTime::Right;
}

using HybridGaugeComplexDensity = std::vector<
    std::map<int, std::map<std::pair<int, std::array<int, 3>>, RI::Tensor<std::complex<double>>>>>;

using HybridDensityBlock = std::tuple<int, int, std::array<int, 3>, std::vector<std::size_t>>;
using HybridDensityLayout = std::vector<std::vector<HybridDensityBlock>>;

inline HybridDensityLayout hybrid_density_layout(const HybridGaugeComplexDensity& density)
{
    HybridDensityLayout layout(density.size());
    for (std::size_t is = 0; is < density.size(); ++is)
    {
        for (const auto& row_atom : density[is])
        {
            for (const auto& block : row_atom.second)
            {
                if (block.second.empty())
                {
                    throw std::invalid_argument("Empty RI density block in common-frame mixing");
                }
                const std::vector<std::size_t> shape(block.second.shape.begin(), block.second.shape.end());
                layout[is].emplace_back(row_atom.first, block.first.first, block.first.second, shape);
            }
        }
    }
    return layout;
}

inline std::vector<std::vector<std::complex<double>>> flatten_hybrid_density(
    const HybridGaugeComplexDensity& density)
{
    std::vector<std::vector<std::complex<double>>> flat(density.size());
    for (std::size_t is = 0; is < density.size(); ++is)
    {
        for (const auto& row_atom : density[is])
        {
            for (const auto& block : row_atom.second)
            {
                const auto& tensor = block.second;
                if (tensor.empty())
                {
                    throw std::invalid_argument("Empty RI density block in common-frame mixing");
                }
                const std::size_t size = tensor.shape.get_shape_all();
                flat[is].insert(flat[is].end(), tensor.ptr(), tensor.ptr() + size);
            }
        }
    }
    return flat;
}

inline void assign_flat_hybrid_density(
    const std::vector<const std::vector<std::complex<double>>*>& flat,
    const HybridDensityLayout& expected_layout,
    HybridGaugeComplexDensity& density)
{
    const auto actual_layout = hybrid_density_layout(density);
    if (flat.size() != density.size() || actual_layout != expected_layout)
    {
        throw std::invalid_argument("RI density topology changed during common-frame mixing");
    }
    // Validate the whole transaction before writing any block.
    for (std::size_t is = 0; is < density.size(); ++is)
    {
        std::size_t expected_size = 0;
        for (const auto& row_atom : density[is])
        {
            for (const auto& block : row_atom.second)
            {
                expected_size += block.second.shape.get_shape_all();
            }
        }
        if (flat[is] == nullptr || flat[is]->size() != expected_size)
        {
            throw std::invalid_argument("RI density size changed during common-frame mixing");
        }
    }
    for (std::size_t is = 0; is < density.size(); ++is)
    {
        std::size_t offset = 0;
        for (auto& row_atom : density[is])
        {
            for (auto& block : row_atom.second)
            {
                const std::size_t size = block.second.shape.get_shape_all();
                // Avoid writing through a shallow copy into accepted-state storage.
                block.second = block.second.copy();
                std::copy(flat[is]->begin() + offset, flat[is]->begin() + offset + size, block.second.ptr());
                offset += size;
            }
        }
    }
}

// Tensor copies share storage in LibRI. Allocate a new phased tensor rather
// than changing entries in place, so cached source density remains untouched.
inline void apply_hybrid_input_phase(const ModuleBase::Matrix3& lattice,
                                     const double lat0,
                                     const ModuleBase::Vector3<double>& source_A,
                                     HybridGaugeComplexDensity& density)
{
    for (int idir = 0; idir < 3; ++idir)
    {
        if (!std::isfinite(source_A[idir]))
        {
            throw std::invalid_argument("Nonfinite RI density source phase");
        }
    }
    for (auto& spin_density : density)
    {
        for (auto& row_atom : spin_density)
        {
            for (auto& block : row_atom.second)
            {
                const auto& cell = block.first.second;
                const ModuleBase::Vector3<double> cell_index(cell[0], cell[1], cell[2]);
                const ModuleBase::Vector3<double> translation = (cell_index * lattice) * lat0;
                const std::complex<double> phase = hybrid_input_phase(source_A, translation);
                block.second = block.second * phase;
            }
        }
    }
}

// The production transaction: transform first, validate topology before mixing,
// then copy the mixed values back without changing shared source storage.
inline void mix_hybrid_density(const ModuleBase::Matrix3& lattice,
                               const double lat0,
                               const ModuleBase::Vector3<double>& source_A,
                               const bool restart,
                               Mix_DMk_2D<std::complex<double>>& mixer,
                               HybridDensityLayout& history_layout,
                               HybridGaugeComplexDensity& density)
{
    if (density.empty())
    {
        throw std::invalid_argument("Missing spin blocks in RI density mixing");
    }
    const auto layout = hybrid_density_layout(density);
    if (!restart && layout != history_layout)
    {
        throw std::invalid_argument("RI density topology changed without a mixing restart");
    }
    apply_hybrid_input_phase(lattice, lat0, source_A, density);
    const auto flat = flatten_hybrid_density(density);
    if (restart)
    {
        mixer.set_nks(static_cast<int>(density.size()));
        history_layout = layout;
    }
    mixer.mix(flat, restart);
    assign_flat_hybrid_density(mixer.get_DMk_out(), history_layout, density);
}
} // namespace RI_2D_Comm

#endif
