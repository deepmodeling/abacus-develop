#ifndef HYBRID_GAUGE_PHASE_H
#define HYBRID_GAUGE_PHASE_H

#include "source_base/vector3.h"
#include "source_base/matrix3.h"
#include "source_lcao/module_ri/mix_dmk_2d.h"

#include <RI/global/Tensor.h>

#include <array>
#include <complex>
#include <map>
#include <utility>
#include <vector>
#include <tuple>

namespace RI_2D_Comm
{
// The RI Fourier convention is exp(-ik.R). Its sHG input therefore needs
// exp(-i A_source.R) to enter the ordinary NAO frame. A_source must identify
// the density's time layer, not merely the current Hamiltonian's stage.
std::complex<double> hybrid_input_phase(
    const ModuleBase::Vector3<double>& vector_potential,
    const ModuleBase::Vector3<double>& lattice_translation);

enum class HybridDensityTime { Left, Right };

// Restarting SCF mixing does not change the time layer of the current density.
HybridDensityTime hybrid_density_time(const int step, const int iteration);

using HybridDensityRow = std::map<std::pair<int, std::array<int, 3>>, RI::Tensor<std::complex<double>>>;
using HybridSpinDensity = std::map<int, HybridDensityRow>;
using HybridGaugeComplexDensity = std::vector<HybridSpinDensity>;

using HybridDensityBlock = std::tuple<int, int, std::array<int, 3>, std::vector<std::size_t>>;
using HybridDensityLayout = std::vector<std::vector<HybridDensityBlock>>;

HybridDensityLayout hybrid_density_layout(const HybridGaugeComplexDensity& density);

std::vector<std::vector<std::complex<double>>> flatten_hybrid_density(
    const HybridGaugeComplexDensity& density);

void assign_flat_hybrid_density(
    const std::vector<const std::vector<std::complex<double>>*>& flat,
    const HybridDensityLayout& expected_layout,
    HybridGaugeComplexDensity& density);

// Tensor copies share storage in LibRI. Allocate a new phased tensor rather
// than changing entries in place, so cached source density remains untouched.
void apply_hybrid_input_phase(const ModuleBase::Matrix3& lattice,
                              const double lat0,
                              const ModuleBase::Vector3<double>& source_A,
                              HybridGaugeComplexDensity& density);

// Project Cartesian A onto the lattice once, before expanding density images.
std::array<double, 3> hybrid_input_phase_per_cell(
    const ModuleBase::Matrix3& lattice, const double lat0,
    const ModuleBase::Vector3<double>& source_A);

// Mix density blocks that have already been converted to the ordinary NAO frame.
void mix_common_frame_density(const bool restart,
                              Mix_DMk_2D<std::complex<double>>& mixer,
                              HybridDensityLayout& history_layout,
                              HybridGaugeComplexDensity& density);

// The production transaction: transform first, validate topology before mixing,
// then copy the mixed values back without changing shared source storage.
void mix_hybrid_density(const ModuleBase::Matrix3& lattice,
                        const double lat0,
                        const ModuleBase::Vector3<double>& source_A,
                        const bool restart,
                        Mix_DMk_2D<std::complex<double>>& mixer,
                        HybridDensityLayout& history_layout,
                        HybridGaugeComplexDensity& density);
} // namespace RI_2D_Comm

#endif
