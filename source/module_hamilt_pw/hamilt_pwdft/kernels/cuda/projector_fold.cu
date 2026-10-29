#include "module_hamilt_pw/hamilt_pwdft/kernels/projector_fold.h"
#include "module_basis/module_pw/pw_basis_k.h"
#include "module_hsolver/kernels/math_kernel_op.h"
#include <base/macros/macros.h>
#include <algorithm>
#include <cstddef>
#include <cassert>
#include <cmath>
#include <limits>
#include <utility>
#include <vector>

namespace hamilt
{
namespace
{
using Complex = std::complex<double>;
using GPU = base_device::DEVICE_GPU;
using CPU = base_device::DEVICE_CPU;
const Complex one(1, 0), zero(0, 0);
constexpr int max_fields = 128;

cuDoubleComplex* device_view(Complex* pointer) { return reinterpret_cast<cuDoubleComplex*>(pointer); }
const cuDoubleComplex* device_view(const Complex* pointer) { return reinterpret_cast<const cuDoubleComplex*>(pointer); }

template <typename T>
class DeviceArray
{
  public:
    T* data = nullptr;
    DeviceArray() = default;
    DeviceArray(const DeviceArray&) = delete;
    DeviceArray& operator=(const DeviceArray&) = delete;
    ~DeviceArray() { base_device::memory::delete_memory_op<T, GPU>()(nullptr, data); }

    void reserve(size_t count)
    {
        if (count > capacity_)
        {
            base_device::memory::resize_memory_op<T, GPU>()(nullptr, data, count);
            capacity_ = count;
        }
    }
    void upload(const std::vector<T>& values)
    {
        reserve(values.size());
        if (!values.empty())
        {
            base_device::memory::synchronize_memory_op<T, GPU, CPU>()(
                nullptr, nullptr, data, values.data(), values.size());
        }
    }
  private:
    size_t capacity_ = 0;
};

__global__ void pack_projectors(const cuDoubleComplex* input, int leading, int reference, int channels,
                                int groups, int rows, const int* mapping, cuDoubleComplex* output)
{
    const size_t index = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (index >= size_t(groups) * rows * channels) { return; }
    const int row = index % rows, channel = (index / rows) % channels;
    const int group = index / (size_t(rows) * channels);
    const int g = mapping[group * rows + row];
    output[index] = g < 0 ? make_cuDoubleComplex(0, 0) : input[size_t(reference * channels + channel) * leading + g];
}

__global__ void pack_wavefunctions(const cuDoubleComplex* input, int leading, int fields, int groups,
                                  int rows, const int* mapping, cuDoubleComplex* output)
{
    const size_t index = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (index >= size_t(groups) * rows * fields) { return; }
    const int row = index % rows, field = (index / rows) % fields;
    const int group = index / (size_t(rows) * fields);
    const int g = mapping[group * rows + row];
    output[index] = g < 0 ? make_cuDoubleComplex(0, 0) : input[size_t(field) * leading + g];
}

__global__ void reorder_groups(const cuDoubleComplex* input, cuDoubleComplex* output, int groups, int channels,
                               int fields, bool to_dense)
{
    const size_t index = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    const size_t columns = size_t(channels) * fields;
    if (index >= columns * groups) { return; }
    const size_t group = index / columns, column = index % columns;
    const size_t dense_index = column * groups + group;
    if (to_dense) { output[dense_index] = input[index]; }
    else { output[index] = input[dense_index]; }
}

__global__ void scatter_projections(const cuDoubleComplex* input, cuDoubleComplex* output, int regular_count,
                                    int channels, int fields, int projectors, const int* atoms)
{
    const size_t index = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (index >= size_t(regular_count) * channels * fields) { return; }
    const int atom = atoms[index % regular_count];
    const size_t column = index / regular_count;
    const int channel = column % channels, field = column / channels;
    output[size_t(field) * projectors + atom * channels + channel] = input[index];
}

__global__ void gather_coefficients(const cuDoubleComplex* input, cuDoubleComplex* output, int input_leading,
                                    int regular_count, int channels, int fields, const int* atoms)
{
    const size_t index = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (index >= size_t(regular_count) * channels * fields) { return; }
    const int atom = atoms[index % regular_count];
    const size_t column = index / regular_count;
    const int channel = column % channels, field = column / channels;
    output[index] = input[size_t(atom * channels + channel) * input_leading + field];
}

__global__ void scatter_wavefunctions(const cuDoubleComplex* input, cuDoubleComplex* output, int plane_waves,
                                      int leading, int rows, int fields, const int* mapping)
{
    const size_t index = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (index >= size_t(plane_waves) * fields) { return; }
    const int g = index % plane_waves, field = index / plane_waves;
    const int packed = mapping[g], group = packed / rows, row = packed % rows;
    const auto value = input[size_t(group) * rows * fields + size_t(field) * rows + row];
    auto& target = output[size_t(field) * leading + g];
    target.x += value.x;
    target.y += value.y;
}

__global__ void pack_irregular_projectors(const cuDoubleComplex* input, cuDoubleComplex* output,
                                          int leading, int plane_waves, int count, const int* indices)
{
    const size_t index = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (index >= size_t(count) * plane_waves) { return; }
    const int projector = index / plane_waves, g = index % plane_waves;
    output[index] = input[size_t(indices[projector]) * leading + g];
}

__global__ void scatter_irregular_projections(const cuDoubleComplex* input, cuDoubleComplex* output,
                                             int count, int fields, int projectors, const int* indices)
{
    const size_t index = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (index >= size_t(count) * fields) { return; }
    const int projector = index % count, field = index / count;
    output[size_t(field) * projectors + indices[projector]] = input[index];
}

__global__ void gather_irregular_coefficients(const cuDoubleComplex* input, cuDoubleComplex* output,
                                             int leading, int count, int fields, const int* indices)
{
    const size_t index = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (index >= size_t(count) * fields) { return; }
    const int projector = index / fields, field = index % fields;
    output[index] = input[size_t(indices[projector]) * leading + field];
}

int positive_remainder(int value, int divisor) { return (value % divisor + divisor) % divisor; }

}

struct ProjectorFold::Impl
{
    struct KPoint
    {
        int plane_waves = 0, rows = 0;
        std::vector<int> packed_to_original, original_to_packed;
        std::vector<Complex> phases;
    };

    const ModulePW::PW_Basis_K& basis;
    const int atoms, channels;
    int nx = 1, ny = 1, nz = 1, reference_atom = 0;
    std::vector<ModuleBase::Vector3<double>> direct, cartesian;
    std::vector<int> regular_atoms, irregular_projectors;
    std::vector<KPoint> points;
    const KPoint* current = nullptr;
    DeviceArray<int> packed_to_original, original_to_packed, regular_indices, irregular_indices;
    DeviceArray<Complex> phases, common, irregular, wave, partial, dense, regular, offsite;

    Impl(int atom_count, int channel_count, const ModuleBase::Vector3<double>* direct_positions,
         const ModuleBase::Vector3<double>* cartesian_positions, const ModulePW::PW_Basis_K& reciprocal_basis)
        : basis(reciprocal_basis), atoms(atom_count), channels(channel_count),
          direct(direct_positions, direct_positions + atoms),
          cartesian(cartesian_positions, cartesian_positions + atoms), points(basis.nks)
    {
    }

    bool aligned(int atom, int x, int y, int z) const
    {
        const auto& position = direct[atom];
        return position.x * x == std::nearbyint(position.x * x)
            && position.y * y == std::nearbyint(position.y * y)
            && position.z * z == std::nearbyint(position.z * z);
    }

    bool select_grid()
    {
        // Binary grid periods preserve exact representability of aligned
        // fractional positions. Other atoms retain their original projectors.
        int best_count = 0, best_volume = std::numeric_limits<int>::max();
        for (int x = 1; x <= 16; x *= 2)
        for (int y = 1; y <= 16; y *= 2)
        for (int z = 1; z <= 16; z *= 2)
        {
            const int volume = x * y * z;
            if (volume > 512) { continue; }
            int count = 0;
            for (int atom = 0; atom < atoms; ++atom) { count += aligned(atom, x, y, z); }
            if (count > best_count || (count == best_count && volume < best_volume))
            {
                best_count = count; best_volume = volume;
                nx = x; ny = y; nz = z;
            }
        }
        if (best_count < 2 || 2 * best_count <= atoms) { return false; }
        for (int atom = 0; atom < atoms; ++atom)
        {
            if (aligned(atom, nx, ny, nz)) { regular_atoms.push_back(atom); }
            else
            {
                for (int channel = 0; channel < channels; ++channel)
                {
                    irregular_projectors.push_back(atom * channels + channel);
                }
            }
        }
        reference_atom = regular_atoms.front();
        return true;
    }

    void initialize_point(int ik)
    {
        auto& point = points[ik];
        point.plane_waves = basis.npwk[ik];
        const int groups = nx * ny * nz;
        std::vector<int> group(point.plane_waves), counts(groups, 0);
        for (int g = 0; g < point.plane_waves; ++g)
        {
            const int isz = basis.getigl2isz(ik, g), is = isz / basis.nz;
            const int ixy = basis.is2fftixy[is];
            int x = ixy / basis.fftny, y = ixy % basis.fftny, z = isz % basis.nz;
            if (x >= basis.nx / 2 + 1) { x -= basis.nx; }
            if (y >= basis.ny / 2 + 1) { y -= basis.ny; }
            if (z >= basis.nz / 2 + 1) { z -= basis.nz; }
            group[g] = (positive_remainder(x, nx) * ny + positive_remainder(y, ny)) * nz + positive_remainder(z, nz);
            ++counts[group[g]];
        }
        point.rows = *std::max_element(counts.begin(), counts.end());
        point.packed_to_original.assign(groups * point.rows, -1);
        point.original_to_packed.resize(point.plane_waves);
        std::fill(counts.begin(), counts.end(), 0);
        for (int g = 0; g < point.plane_waves; ++g)
        {
            const int index = group[g] * point.rows + counts[group[g]]++;
            point.original_to_packed[g] = index;
            point.packed_to_original[index] = g;
        }
        point.phases.resize(groups * regular_atoms.size());
        const double two_pi = 2 * std::acos(-1.0);
        for (int q = 0; q < groups; ++q)
        {
            const int x = q / (ny * nz), y = (q / nz) % ny, z = q % nz;
            for (size_t index = 0; index < regular_atoms.size(); ++index)
            {
                // P_a(G) = P_ref(G) exp[-2 pi i (G+k).(R_a-R_ref)].
                // The conjugate phase is constant within each reciprocal
                // residue group and supplies the final projection factor.
                const int atom = regular_atoms[index];
                const auto position = direct[atom] - direct[reference_atom];
                const auto cartesian_delta = cartesian[atom] - cartesian[reference_atom];
                const double angle = two_pi * (basis.kvec_c[ik] * cartesian_delta
                                               + x * position.x + y * position.y + z * position.z);
                point.phases[q * regular_atoms.size() + index] = std::polar(1.0, angle);
            }
        }
    }

    void prepare(int ik, const Complex* projectors, int leading)
    {
        assert(ik >= 0 && ik < basis.nks);
        if (points[ik].packed_to_original.empty()) { initialize_point(ik); }
        const auto& point = points[ik];
        const int groups = nx * ny * nz;
        const size_t work = size_t(groups) * point.rows + size_t(groups) * regular_atoms.size()
                            + size_t(irregular_projectors.size() / channels) * point.plane_waves;
        // Packing and phase transforms are amortized only when the contraction
        // work is less than half that of the complete projector matrix.
        current = 2 * work < size_t(atoms) * point.plane_waves ? &point : nullptr;
        if (current == nullptr) { return; }
        assert(leading >= point.plane_waves);
        packed_to_original.upload(point.packed_to_original);
        original_to_packed.upload(point.original_to_packed);
        phases.upload(point.phases);
        regular_indices.upload(regular_atoms);
        irregular_indices.upload(irregular_projectors);
        const size_t common_size = size_t(groups) * point.rows * channels;
        common.reserve(common_size);
        pack_projectors<<<(common_size + 255) / 256, 256>>>(device_view(projectors), leading, reference_atom,
            channels, groups, point.rows, packed_to_original.data, device_view(common.data));
        if (!irregular_projectors.empty())
        {
            const size_t size = size_t(point.plane_waves) * irregular_projectors.size();
            irregular.reserve(size);
            pack_irregular_projectors<<<(size + 255) / 256, 256>>>(device_view(projectors), device_view(irregular.data),
                leading, point.plane_waves, irregular_projectors.size(), irregular_indices.data);
        }
        cudaErrcheck(cudaGetLastError());
    }

    void reserve_workspace(int fields)
    {
        const int groups = nx * ny * nz;
        wave.reserve(size_t(groups) * current->rows * fields);
        partial.reserve(size_t(groups) * channels * fields);
        dense.reserve(size_t(groups) * channels * fields);
        regular.reserve(regular_atoms.size() * channels * fields);
        offsite.reserve(irregular_projectors.size() * fields);
    }

    void project_tile(int fields, int leading, const Complex* input, Complex* output)
    {
        const int groups = nx * ny * nz, rows = current->rows, regular_count = regular_atoms.size();
        const int off_count = irregular_projectors.size(), projectors = atoms * channels;
        const size_t wave_size = size_t(groups) * rows * fields;
        const size_t partial_size = size_t(groups) * channels * fields;
        const size_t regular_size = size_t(regular_count) * channels * fields;
        pack_wavefunctions<<<(wave_size + 255) / 256, 256>>>(device_view(input), leading, fields, groups,
            rows, packed_to_original.data, device_view(wave.data));
        hsolver::gemm_strided_batched_gpu('C', 'N', channels, fields, rows, one,
            common.data, rows, static_cast<long long>(rows) * channels,
            wave.data, rows, static_cast<long long>(rows) * fields, zero,
            partial.data, channels, static_cast<long long>(channels) * fields, groups);
        reorder_groups<<<(partial_size + 255) / 256, 256>>>(device_view(partial.data), device_view(dense.data),
            groups, channels, fields, true);
        hsolver::gemm_op<Complex, GPU>()(nullptr, 'N', 'N', regular_count, channels * fields, groups,
            &one, phases.data, regular_count, dense.data, groups, &zero, regular.data, regular_count);
        scatter_projections<<<(regular_size + 255) / 256, 256>>>(device_view(regular.data), device_view(output),
            regular_count, channels, fields, projectors, regular_indices.data);
        if (off_count != 0)
        {
            hsolver::gemm_op<Complex, GPU>()(nullptr, 'C', 'N', off_count, fields, current->plane_waves,
                &one, irregular.data, current->plane_waves, input, leading, &zero, offsite.data, off_count);
            const size_t size = size_t(off_count) * fields;
            scatter_irregular_projections<<<(size + 255) / 256, 256>>>(device_view(offsite.data), device_view(output),
                off_count, fields, projectors, irregular_indices.data);
        }
    }

    void add_tile(int fields, int coefficient_leading, int output_leading, const Complex* coefficients, Complex* output)
    {
        const int groups = nx * ny * nz, rows = current->rows, regular_count = regular_atoms.size();
        const int off_count = irregular_projectors.size();
        const size_t regular_size = size_t(regular_count) * channels * fields;
        const size_t partial_size = size_t(groups) * channels * fields;
        gather_coefficients<<<(regular_size + 255) / 256, 256>>>(device_view(coefficients), device_view(regular.data),
            coefficient_leading, regular_count, channels, fields, regular_indices.data);
        hsolver::gemm_op<Complex, GPU>()(nullptr, 'C', 'N', groups, channels * fields, regular_count,
            &one, phases.data, regular_count, regular.data, regular_count, &zero, dense.data, groups);
        reorder_groups<<<(partial_size + 255) / 256, 256>>>(device_view(dense.data), device_view(partial.data),
            groups, channels, fields, false);
        hsolver::gemm_strided_batched_gpu('N', 'N', rows, fields, channels, one,
            common.data, rows, static_cast<long long>(rows) * channels,
            partial.data, channels, static_cast<long long>(channels) * fields, zero,
            wave.data, rows, static_cast<long long>(rows) * fields, groups);
        const size_t wave_size = size_t(current->plane_waves) * fields;
        scatter_wavefunctions<<<(wave_size + 255) / 256, 256>>>(device_view(wave.data), device_view(output),
            current->plane_waves, output_leading, rows, fields, original_to_packed.data);
        if (off_count != 0)
        {
            const size_t size = size_t(off_count) * fields;
            gather_irregular_coefficients<<<(size + 255) / 256, 256>>>(device_view(coefficients), device_view(offsite.data),
                coefficient_leading, off_count, fields, irregular_indices.data);
            hsolver::gemm_op<Complex, GPU>()(nullptr, 'N', 'T', current->plane_waves, fields, off_count,
                &one, irregular.data, current->plane_waves, offsite.data, fields, &one, output, output_leading);
        }
    }
};

ProjectorFold::ProjectorFold(std::unique_ptr<Impl> implementation) : implementation_(std::move(implementation)) {}
ProjectorFold::~ProjectorFold() = default;

std::unique_ptr<ProjectorFold> ProjectorFold::create(int atoms, int channels,
    const ModuleBase::Vector3<double>* direct_positions, const ModuleBase::Vector3<double>* cartesian_positions,
    const ModulePW::PW_Basis_K& basis)
{
    if (atoms < 2 || channels <= 0 || basis.nks <= 0) { return {}; }
    auto implementation = std::make_unique<Impl>(atoms, channels, direct_positions, cartesian_positions, basis);
    if (!implementation->select_grid()) { return {}; }
    return std::unique_ptr<ProjectorFold>(new ProjectorFold(std::move(implementation)));
}

void ProjectorFold::prepare(int ik, const Complex* projectors, int projector_leading)
{
    implementation_->prepare(ik, projectors, projector_leading);
}

bool ProjectorFold::project(int columns, int leading, const Complex* input, Complex* output)
{
    auto& state = *implementation_;
    if (state.current == nullptr) { return false; }
    assert(leading >= state.current->plane_waves);
    state.reserve_workspace(std::min(columns, max_fields));
    for (int first = 0; first < columns; first += max_fields)
    {
        const int count = std::min(max_fields, columns - first);
        state.project_tile(count, leading, input + size_t(first) * leading,
                           output + size_t(first) * state.atoms * state.channels);
    }
    cudaErrcheck(cudaGetLastError());
    return true;
}

bool ProjectorFold::add(int columns, int leading, const Complex* coefficients, Complex* output)
{
    auto& state = *implementation_;
    if (state.current == nullptr) { return false; }
    assert(leading >= state.current->plane_waves);
    state.reserve_workspace(std::min(columns, max_fields));
    for (int first = 0; first < columns; first += max_fields)
    {
        const int count = std::min(max_fields, columns - first);
        state.add_tile(count, columns, leading, coefficients + first, output + size_t(first) * leading);
    }
    cudaErrcheck(cudaGetLastError());
    return true;
}

}
