#include "hybrid_gauge_phase.h"
#include "source_base/matrix3.h"
#include "source_base/vector3.h"
#include "source_lcao/module_ri/mix_dmk_2d.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace RI_2D_Comm
{
std::complex<double> hybrid_input_phase(
    const ModuleBase::Vector3<double>& vector_potential,
    const ModuleBase::Vector3<double>& lattice_translation)
{
    const double argument = -(vector_potential * lattice_translation);
    return std::exp(std::complex<double>(0.0, argument));
}

HybridDensityTime hybrid_density_time(const int step, const int iteration)
{
    return step > 0 && iteration == 1 ? HybridDensityTime::Left : HybridDensityTime::Right;
}

HybridDensityLayout hybrid_density_layout(const HybridGaugeComplexDensity& density)
{
    HybridDensityLayout layout(density.size());
    for (std::size_t is = 0; is < density.size(); ++is)
    {
        for (const HybridSpinDensity::value_type& row_atom : density[is])
        {
            for (const HybridDensityRow::value_type& block : row_atom.second)
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

std::vector<std::vector<std::complex<double>>> flatten_hybrid_density(
    const HybridGaugeComplexDensity& density)
{
    std::vector<std::vector<std::complex<double>>> flat(density.size());
    for (std::size_t is = 0; is < density.size(); ++is)
    {
        for (const HybridSpinDensity::value_type& row_atom : density[is])
        {
            for (const HybridDensityRow::value_type& block : row_atom.second)
            {
                const RI::Tensor<std::complex<double>>& tensor = block.second;
                if (tensor.empty())
                {
                    throw std::invalid_argument("Empty RI density block in common-frame mixing");
                }
                const std::size_t size = tensor.shape.get_shape_all();
                const std::complex<double>* block_begin = tensor.ptr();
                const std::complex<double>* block_end = block_begin + size;
                flat[is].insert(flat[is].end(), block_begin, block_end);
            }
        }
    }
    return flat;
}

void assign_flat_hybrid_density(
    const std::vector<const std::vector<std::complex<double>>*>& flat,
    const HybridDensityLayout& expected_layout,
    HybridGaugeComplexDensity& density)
{
    const HybridDensityLayout actual_layout = hybrid_density_layout(density);
    if (flat.size() != density.size() || actual_layout != expected_layout)
    {
        throw std::invalid_argument("RI density topology changed during common-frame mixing");
    }
    // Validate the whole transaction before writing any block.
    for (std::size_t is = 0; is < density.size(); ++is)
    {
        std::size_t expected_size = 0;
        for (const HybridSpinDensity::value_type& row_atom : density[is])
        {
            for (const HybridDensityRow::value_type& block : row_atom.second)
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
        for (HybridSpinDensity::value_type& row_atom : density[is])
        {
            for (HybridDensityRow::value_type& block : row_atom.second)
            {
                const std::size_t size = block.second.shape.get_shape_all();
                // Avoid writing through a shallow copy into accepted-state storage.
                block.second = block.second.copy();
                const std::vector<std::complex<double>>::const_iterator block_begin = flat[is]->begin() + offset;
                const std::vector<std::complex<double>>::const_iterator block_end = block_begin + size;
                std::copy(block_begin, block_end, block.second.ptr());
                offset += size;
            }
        }
    }
}

void apply_hybrid_input_phase(const ModuleBase::Matrix3& lattice,
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
    for (HybridSpinDensity& spin_density : density)
    {
        for (HybridSpinDensity::value_type& row_atom : spin_density)
        {
            for (HybridDensityRow::value_type& block : row_atom.second)
            {
                const std::array<int, 3>& cell = block.first.second;
                const ModuleBase::Vector3<double> cell_index(cell[0], cell[1], cell[2]);
                const ModuleBase::Vector3<double> translation = (cell_index * lattice) * lat0;
                const std::complex<double> phase = hybrid_input_phase(source_A, translation);
                block.second = block.second * phase;
            }
        }
    }
}

std::array<double, 3> hybrid_input_phase_per_cell(
    const ModuleBase::Matrix3& lattice, const double lat0,
    const ModuleBase::Vector3<double>& source_A)
{
    std::array<double, 3> phase;
    for (int d = 0; d < 3; ++d)
    {
        if (!std::isfinite(source_A[d]))
        {
            throw std::invalid_argument("Nonfinite RI density source phase");
        }
        ModuleBase::Vector3<double> unit_cell(0.0, 0.0, 0.0);
        unit_cell[d] = 1.0;
        const ModuleBase::Vector3<double> translation = (unit_cell * lattice) * lat0;
        phase[d] = source_A * translation;
    }
    return phase;
}

void mix_common_frame_density(const bool restart,
                              Mix_DMk_2D<std::complex<double>>& mixer,
                              HybridDensityLayout& history_layout,
                              HybridGaugeComplexDensity& density)
{
    if (density.empty())
    {
        throw std::invalid_argument("Missing spin blocks in RI density mixing");
    }
    const HybridDensityLayout layout = hybrid_density_layout(density);
    if (!restart && layout != history_layout)
    {
        throw std::invalid_argument("RI density topology changed without a mixing restart");
    }
    const std::vector<std::vector<std::complex<double>>> flat = flatten_hybrid_density(density);
    if (restart)
    {
        mixer.set_nks(static_cast<int>(density.size()));
        history_layout = layout;
    }
    mixer.mix(flat, restart);
    assign_flat_hybrid_density(mixer.get_DMk_out(), history_layout, density);
}

void mix_hybrid_density(const ModuleBase::Matrix3& lattice,
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
    const HybridDensityLayout layout = hybrid_density_layout(density);
    if (!restart && layout != history_layout)
    {
        throw std::invalid_argument("RI density topology changed without a mixing restart");
    }
    apply_hybrid_input_phase(lattice, lat0, source_A, density);
    mix_common_frame_density(restart, mixer, history_layout, density);
}
} // namespace RI_2D_Comm
