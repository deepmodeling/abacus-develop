#include "source_lcao/module_ri/sternheimer_grid_transfer.h"

#include <fftw3.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <limits>
#include <memory>
#include <new>
#include <stdexcept>
#include <utility>
#include <vector>

namespace
{

using Grid = ModuleRI::SternheimerGridTransfer::Grid;
using Vector = ModuleRI::SternheimerGridTransfer::Vector;
using Complex = ModuleRI::SternheimerFDHamiltonian::Complex;
using Lattice = ModuleRI::SternheimerFDLatticeVectors;

struct GridInfo
{
    std::size_t count;
    Lattice lattice;
};

GridInfo validate_grid(const Grid& grid)
{
    if (!grid.periodic)
    {
        throw std::invalid_argument("Sternheimer grid transfer requires periodic grids.");
    }
    std::size_t count = 1;
    for (const int dimension : {grid.nx, grid.ny, grid.nz})
    {
        if (dimension <= 0
            || count > static_cast<std::size_t>(std::numeric_limits<int>::max())
                           / static_cast<std::size_t>(dimension))
        {
            throw std::invalid_argument("Sternheimer grid transfer has invalid or overflowing grid dimensions.");
        }
        count *= static_cast<std::size_t>(dimension);
    }
    for (const double spacing : {grid.hx, grid.hy, grid.hz})
    {
        if (!std::isfinite(spacing) || spacing <= 0.0)
        {
            throw std::invalid_argument("Sternheimer grid transfer requires finite positive grid spacings.");
        }
    }
    for (const double k : grid.kpoint)
    {
        if (!std::isfinite(k))
        {
            throw std::invalid_argument("Sternheimer grid transfer requires a finite Bloch k point.");
        }
    }
    const Lattice cell = ModuleRI::sternheimer_fd_grid_lattice_vectors(grid);
    for (const auto& row : cell)
    {
        for (const double entry : row)
        {
            if (!std::isfinite(entry))
            {
                throw std::invalid_argument("Sternheimer grid transfer requires a finite lattice.");
            }
        }
    }
    const double determinant
        = cell[0][0] * (cell[1][1] * cell[2][2] - cell[1][2] * cell[2][1])
          - cell[0][1] * (cell[1][0] * cell[2][2] - cell[1][2] * cell[2][0])
          + cell[0][2] * (cell[1][0] * cell[2][1] - cell[1][1] * cell[2][0]);
    if (!std::isfinite(determinant) || determinant == 0.0)
    {
        throw std::invalid_argument("Sternheimer grid transfer requires a nonsingular finite-volume lattice.");
    }
    return {count, cell};
}

struct PairInfo
{
    std::size_t coarse_count;
    std::size_t fine_count;
    bool identity;
};

PairInfo validate_pair(const Grid& coarse, const Grid& fine)
{
    const GridInfo coarse_info = validate_grid(coarse);
    const GridInfo fine_info = validate_grid(fine);
    if (fine.nx < coarse.nx || fine.ny < coarse.ny || fine.nz < coarse.nz)
    {
        throw std::invalid_argument("Sternheimer grid transfer fine dimensions must not be smaller than coarse.");
    }
    constexpr double tolerance = 1.0e-12;
    for (int i = 0; i != 3; ++i)
    {
        if (std::abs(coarse.kpoint[i] - fine.kpoint[i]) > tolerance)
        {
            throw std::invalid_argument("Sternheimer grid transfer requires the same explicit Bloch k point.");
        }
        for (int j = 0; j != 3; ++j)
        {
            const double left = coarse_info.lattice[i][j];
            const double right = fine_info.lattice[i][j];
            if (std::abs(left - right) > tolerance * std::max({1.0, std::abs(left), std::abs(right)}))
            {
                throw std::invalid_argument("Sternheimer grid transfer requires the same physical lattice.");
            }
        }
    }
    return {coarse_info.count, fine_info.count,
            coarse.nx == fine.nx && coarse.ny == fine.ny && coarse.nz == fine.nz};
}

std::size_t numeric_workspace_bytes(const PairInfo& grids, const bool spectral = false)
{
    if (grids.identity && !spectral)
    {
        return 0;
    }
    std::size_t total = 0;
    const auto append = [&total](const std::size_t count, const std::size_t bytes_per_item) {
        if (count > (std::numeric_limits<std::size_t>::max() - total) / bytes_per_item)
        {
            throw std::invalid_argument("Sternheimer grid transfer workspace size overflows.");
        }
        total += count * bytes_per_item;
    };
    append(grids.coarse_count, sizeof(fftw_complex) + sizeof(Complex));
    append(grids.fine_count, sizeof(fftw_complex) + sizeof(Complex));
    append(grids.coarse_count, sizeof(std::size_t));
    return total;
}

std::size_t index(const Grid& grid, const int ix, const int iy, const int iz)
{
    return (static_cast<std::size_t>(ix) * grid.ny + iy) * grid.nz + iz;
}

int embedded_index(const int coarse_index, const int coarse_dimension, const int fine_dimension)
{
    const int frequency = coarse_index <= coarse_dimension / 2 ? coarse_index : coarse_index - coarse_dimension;
    return frequency < 0 ? frequency + fine_dimension : frequency;
}

Lattice spectral_reciprocal_vectors(const Grid& grid)
{
    Lattice reciprocal = ModuleRI::sternheimer_fd_grid_dual_vectors(grid);
    const double two_pi = 2.0 * std::acos(-1.0);
    for (auto& row : reciprocal)
    {
        for (double& value : row)
        {
            value *= two_pi;
        }
    }
    const std::array<int, 3> dimensions{grid.nx, grid.ny, grid.nz};
    double maximum_squared = 0.0;
    for (int axis = 0; axis != 3; ++axis)
    {
        double maximum = 0.0;
        for (int direction = 0; direction != 3; ++direction)
        {
            maximum += std::abs(reciprocal[direction][axis])
                       * (std::abs(grid.kpoint[direction]) + dimensions[direction] / 2);
        }
        maximum_squared += maximum * maximum;
    }
    if (!std::isfinite(maximum_squared))
    {
        throw std::invalid_argument("Sternheimer grid transfer requires finite spectral derivative symbols.");
    }
    return reciprocal;
}

std::array<double, 3> spectral_wavevector(const Grid& grid, const Lattice& reciprocal, const std::size_t i)
{
    const std::array<int, 3> dimensions{grid.nx, grid.ny, grid.nz};
    const std::array<int, 3> indices{static_cast<int>(i / (static_cast<std::size_t>(grid.ny) * grid.nz)),
                                     static_cast<int>((i / grid.nz) % grid.ny),
                                     static_cast<int>(i % grid.nz)};
    std::array<double, 3> result{};
    for (int direction = 0; direction != 3; ++direction)
    {
        const int frequency = indices[direction] <= dimensions[direction] / 2
                                  ? indices[direction] : indices[direction] - dimensions[direction];
        // Keep the explicit k label: reducing k modulo a grid dimension would
        // preserve nodal phases but give the wrong analytic derivative.
        const double shifted = frequency + grid.kpoint[direction];
        for (int axis = 0; axis != 3; ++axis)
        {
            result[axis] += shifted * reciprocal[direction][axis];
        }
    }
    return result;
}

void validate_spectral_input(const Vector& input, const std::size_t count)
{
    if (input.size() != count)
    {
        throw std::invalid_argument("Sternheimer spectral transfer input size does not match its grid.");
    }
    for (const Complex value : input)
    {
        if (!std::isfinite(value.real()) || !std::isfinite(value.imag()))
        {
            throw std::invalid_argument("Sternheimer spectral transfer requires finite input fields.");
        }
    }
}

Vector bloch_phases(const Grid& grid, const std::size_t count)
{
    Vector phases(count);
    const double two_pi = 2.0 * std::acos(-1.0);
    // Integer multiples of a dimension do not change phases at that grid's nodes.
    const std::array<double, 3> k{std::remainder(grid.kpoint[0], static_cast<double>(grid.nx)),
                                 std::remainder(grid.kpoint[1], static_cast<double>(grid.ny)),
                                 std::remainder(grid.kpoint[2], static_cast<double>(grid.nz))};
    for (int ix = 0; ix != grid.nx; ++ix)
    {
        for (int iy = 0; iy != grid.ny; ++iy)
        {
            for (int iz = 0; iz != grid.nz; ++iz)
            {
                const double phase = two_pi * (k[0] * (static_cast<double>(ix) / grid.nx)
                                                + k[1] * (static_cast<double>(iy) / grid.ny)
                                                + k[2] * (static_cast<double>(iz) / grid.nz));
                phases[index(grid, ix, iy, iz)] = std::exp(Complex(0.0, phase));
            }
        }
    }
    return phases;
}

struct FFTGrid
{
    fftw_complex* data = nullptr;
    fftw_plan forward = nullptr;
    fftw_plan backward = nullptr;

    FFTGrid(const Grid& grid, const std::size_t count)
    {
        data = fftw_alloc_complex(count);
        if (data == nullptr)
        {
            throw std::bad_alloc();
        }
        for (std::size_t i = 0; i != count; ++i)
        {
            data[i][0] = 0.0;
            data[i][1] = 0.0;
        }
#ifdef _OPENMP
#pragma omp critical(sternheimer_fftw_plan_management)
#endif
        {
            forward = fftw_plan_dft_3d(grid.nx, grid.ny, grid.nz, data, data, FFTW_FORWARD, FFTW_ESTIMATE);
            backward = fftw_plan_dft_3d(grid.nx, grid.ny, grid.nz, data, data, FFTW_BACKWARD, FFTW_ESTIMATE);
        }
        if (forward == nullptr || backward == nullptr)
        {
            release();
            throw std::runtime_error("Failed to create Sternheimer grid transfer FFTW plans.");
        }
    }

    ~FFTGrid()
    {
        release();
    }

    FFTGrid(const FFTGrid&) = delete;
    FFTGrid& operator=(const FFTGrid&) = delete;

    void release()
    {
#ifdef _OPENMP
#pragma omp critical(sternheimer_fftw_plan_management)
#endif
        {
            if (forward != nullptr)
            {
                fftw_destroy_plan(forward);
            }
            if (backward != nullptr)
            {
                fftw_destroy_plan(backward);
            }
        }
        fftw_free(data);
        data = nullptr;
        forward = nullptr;
        backward = nullptr;
    }
};

} // namespace

namespace ModuleRI
{

struct SternheimerGridTransfer::Impl
{
    Grid coarse_grid;
    Grid fine_grid;
    PairInfo grids;
    Vector coarse_phases;
    Vector fine_phases;
    std::vector<std::size_t> coarse_to_fine;
    std::unique_ptr<FFTGrid> coarse_fft;
    std::unique_ptr<FFTGrid> fine_fft;

    Impl(const Grid& coarse, const Grid& fine)
        : coarse_grid(coarse), fine_grid(fine), grids(validate_pair(coarse, fine))
    {
        numeric_workspace_bytes(grids);
        if (grids.identity)
        {
            return;
        }
        initialize_fft();
    }

    void initialize_fft()
    {
        if (coarse_fft != nullptr)
        {
            return;
        }
        // Build once, including lazy same-grid spectral use. Locals ensure a
        // failed FFT allocation leaves no partial workspace to duplicate on retry.
        const Grid& coarse = coarse_grid;
        const Grid& fine = fine_grid;
        Vector new_coarse_phases = bloch_phases(coarse, grids.coarse_count);
        Vector new_fine_phases = bloch_phases(fine, grids.fine_count);
        std::vector<std::size_t> new_map(grids.coarse_count);
        for (int ix = 0; ix != coarse.nx; ++ix)
        {
            for (int iy = 0; iy != coarse.ny; ++iy)
            {
                for (int iz = 0; iz != coarse.nz; ++iz)
                {
                    new_map[index(coarse, ix, iy, iz)]
                        = index(fine, embedded_index(ix, coarse.nx, fine.nx),
                                embedded_index(iy, coarse.ny, fine.ny),
                                embedded_index(iz, coarse.nz, fine.nz));
                }
            }
        }
        auto new_coarse_fft = std::make_unique<FFTGrid>(coarse, grids.coarse_count);
        auto new_fine_fft = std::make_unique<FFTGrid>(fine, grids.fine_count);
        coarse_phases = std::move(new_coarse_phases);
        fine_phases = std::move(new_fine_phases);
        coarse_to_fine = std::move(new_map);
        coarse_fft = std::move(new_coarse_fft);
        fine_fft = std::move(new_fine_fft);
    }

    Lattice prepare_spectral(const std::size_t budget)
    {
        const Lattice reciprocal = spectral_reciprocal_vectors(coarse_grid);
        if (numeric_workspace_bytes(grids, true) > budget)
        {
            throw std::length_error("Sternheimer spectral transfer numeric workspace exceeds its budget.");
        }
        initialize_fft();
        return reciprocal;
    }

    static void forward_field(const Vector& input, const Vector& phases, FFTGrid& fft)
    {
        for (std::size_t i = 0; i != input.size(); ++i)
        {
            const Complex periodic = std::conj(phases[i]) * input[i];
            fft.data[i][0] = periodic.real();
            fft.data[i][1] = periodic.imag();
        }
        fftw_execute(fft.forward);
    }

    static void backward_field(FFTGrid& fft, const Vector& phases, Vector& output)
    {
        fftw_execute(fft.backward);
        output.resize(phases.size());
        for (std::size_t i = 0; i != output.size(); ++i)
        {
            output[i] = phases[i] * Complex(fft.data[i][0], fft.data[i][1]);
            if (!std::isfinite(output[i].real()) || !std::isfinite(output[i].imag()))
            {
                throw std::overflow_error("Sternheimer spectral transfer produced a nonfinite result.");
            }
        }
    }

    void interpolate_with_gradients(const Vector& input, Vector& values,
                                    std::array<Vector, 3>& gradients, const std::size_t budget)
    {
        for (const auto& field : gradients)
        {
            if (&values == &field)
            {
                throw std::invalid_argument("Sternheimer fine values and gradient outputs must be distinct vectors.");
            }
        }
        validate_spectral_input(input, grids.coarse_count);
        const Lattice reciprocal = prepare_spectral(budget);
        // Preserve the input spectrum before writing any possibly aliased output.
        forward_field(input, coarse_phases, *coarse_fft);
        const double normalization = 1.0 / static_cast<double>(grids.coarse_count);
        for (int component = -1; component != 3; ++component)
        {
            if (component == -1 && grids.identity)
            {
                values = input;
                continue;
            }
            for (std::size_t i = 0; i != grids.fine_count; ++i)
            {
                fine_fft->data[i][0] = 0.0;
                fine_fft->data[i][1] = 0.0;
            }
            for (std::size_t i = 0; i != grids.coarse_count; ++i)
            {
                const Complex coefficient(coarse_fft->data[i][0], coarse_fft->data[i][1]);
                const Complex factor = component == -1 ? Complex(normalization, 0.0)
                    : Complex(0.0, normalization * spectral_wavevector(coarse_grid, reciprocal, i)[component]);
                const Complex transformed = factor * coefficient;
                const std::size_t target = coarse_to_fine[i];
                fine_fft->data[target][0] = transformed.real();
                fine_fft->data[target][1] = transformed.imag();
            }
            backward_field(*fine_fft, fine_phases, component == -1 ? values : gradients[component]);
        }
    }

    void restrict_gradient_adjoint(const std::array<Vector, 3>& fields, Vector& output, const std::size_t budget)
    {
        for (const auto& field : fields)
        {
            validate_spectral_input(field, grids.fine_count);
        }
        const Lattice reciprocal = prepare_spectral(budget);
        for (std::size_t i = 0; i != grids.coarse_count; ++i)
        {
            coarse_fft->data[i][0] = 0.0;
            coarse_fft->data[i][1] = 0.0;
        }
        const double normalization = 1.0 / static_cast<double>(grids.fine_count);
        for (int component = 0; component != 3; ++component)
        {
            forward_field(fields[component], fine_phases, *fine_fft);
            for (std::size_t i = 0; i != grids.coarse_count; ++i)
            {
                const std::size_t source = coarse_to_fine[i];
                const Complex coefficient(fine_fft->data[source][0], fine_fft->data[source][1]);
                const Complex factor(0.0, -normalization * spectral_wavevector(coarse_grid, reciprocal, i)[component]);
                const Complex transformed = factor * coefficient;
                coarse_fft->data[i][0] += transformed.real();
                coarse_fft->data[i][1] += transformed.imag();
            }
        }
        // Consume all three input fields before an aliased output can overwrite one.
        backward_field(*coarse_fft, coarse_phases, output);
    }

    void apply_negative_laplacian(const Vector& input, Vector& output, const std::size_t budget)
    {
        validate_spectral_input(input, grids.coarse_count);
        const Lattice reciprocal = prepare_spectral(budget);
        forward_field(input, coarse_phases, *coarse_fft);
        const double normalization = 1.0 / static_cast<double>(grids.coarse_count);
        for (std::size_t i = 0; i != grids.coarse_count; ++i)
        {
            const auto q = spectral_wavevector(coarse_grid, reciprocal, i);
            const double factor = normalization * (q[0] * q[0] + q[1] * q[1] + q[2] * q[2]);
            coarse_fft->data[i][0] *= factor;
            coarse_fft->data[i][1] *= factor;
        }
        backward_field(*coarse_fft, coarse_phases, output);
    }

    void apply(const Vector& input, Vector& output, const bool interpolate)
    {
        const std::size_t input_size = interpolate ? grids.coarse_count : grids.fine_count;
        const std::size_t output_size = interpolate ? grids.fine_count : grids.coarse_count;
        if (input.size() != input_size)
        {
            throw std::invalid_argument("Sternheimer grid transfer input vector size does not match its grid.");
        }
        if (grids.identity)
        {
            output = input;
            return;
        }
        FFTGrid& source = interpolate ? *coarse_fft : *fine_fft;
        FFTGrid& target = interpolate ? *fine_fft : *coarse_fft;
        const Vector& source_phases = interpolate ? coarse_phases : fine_phases;
        const Vector& target_phases = interpolate ? fine_phases : coarse_phases;
        for (std::size_t i = 0; i != input_size; ++i)
        {
            const Complex periodic = std::conj(source_phases[i]) * input[i];
            source.data[i][0] = periodic.real();
            source.data[i][1] = periodic.imag();
        }
        fftw_execute(source.forward);
        for (std::size_t i = 0; i != output_size; ++i)
        {
            target.data[i][0] = 0.0;
            target.data[i][1] = 0.0;
        }
        // FFTW is unnormalized: 1/Nc gives J, while 1/Nf gives its dV-weighted adjoint.
        const double normalization = 1.0 / static_cast<double>(input_size);
        for (std::size_t i = 0; i != grids.coarse_count; ++i)
        {
            const std::size_t from = interpolate ? i : coarse_to_fine[i];
            const std::size_t to = interpolate ? coarse_to_fine[i] : i;
            target.data[to][0] = normalization * source.data[from][0];
            target.data[to][1] = normalization * source.data[from][1];
        }
        fftw_execute(target.backward);
        output.resize(output_size);
        for (std::size_t i = 0; i != output_size; ++i)
        {
            output[i] = target_phases[i] * Complex(target.data[i][0], target.data[i][1]);
        }
    }
};

SternheimerGridTransfer::SternheimerGridTransfer(const Grid& coarse_grid, const Grid& fine_grid)
    : impl_(std::make_unique<Impl>(coarse_grid, fine_grid))
{
}

SternheimerGridTransfer::~SternheimerGridTransfer() = default;
SternheimerGridTransfer::SternheimerGridTransfer(SternheimerGridTransfer&&) noexcept = default;
SternheimerGridTransfer& SternheimerGridTransfer::operator=(SternheimerGridTransfer&&) noexcept = default;

const SternheimerGridTransfer::Grid& SternheimerGridTransfer::coarse_grid() const
{
    return impl_->coarse_grid;
}

const SternheimerGridTransfer::Grid& SternheimerGridTransfer::fine_grid() const
{
    return impl_->fine_grid;
}

void SternheimerGridTransfer::interpolate(const Vector& coarse, Vector& fine)
{
    impl_->apply(coarse, fine, true);
}

void SternheimerGridTransfer::restrict_adjoint(const Vector& fine, Vector& coarse)
{
    impl_->apply(fine, coarse, false);
}

void SternheimerGridTransfer::interpolate_with_gradients(const Vector& coarse, Vector& fine_values,
                                                        std::array<Vector, 3>& fine_gradients,
                                                        const std::size_t max_workspace_bytes)
{
    impl_->interpolate_with_gradients(coarse, fine_values, fine_gradients, max_workspace_bytes);
}

void SternheimerGridTransfer::restrict_gradient_adjoint(const std::array<Vector, 3>& fine_fields, Vector& coarse,
                                                       const std::size_t max_workspace_bytes)
{
    impl_->restrict_gradient_adjoint(fine_fields, coarse, max_workspace_bytes);
}

void SternheimerGridTransfer::apply_negative_laplacian(const Vector& coarse, Vector& result,
                                                      const std::size_t max_workspace_bytes)
{
    impl_->apply_negative_laplacian(coarse, result, max_workspace_bytes);
}

std::size_t SternheimerGridTransfer::workspace_bytes_required(const Grid& coarse_grid, const Grid& fine_grid)
{
    return numeric_workspace_bytes(validate_pair(coarse_grid, fine_grid));
}

std::size_t SternheimerGridTransfer::gradient_workspace_bytes_required(const Grid& coarse_grid, const Grid& fine_grid)
{
    const PairInfo grids = validate_pair(coarse_grid, fine_grid);
    spectral_reciprocal_vectors(coarse_grid);
    return numeric_workspace_bytes(grids, true);
}

} // namespace ModuleRI
