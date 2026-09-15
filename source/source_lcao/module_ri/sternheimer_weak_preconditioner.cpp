#include "source_lcao/module_ri/sternheimer_weak_preconditioner.h"

#include <fftw3.h>

#include <array>
#include <cmath>
#include <complex>
#include <limits>
#include <new>
#include <stdexcept>
#include <utility>
#include <vector>

namespace
{

using Preconditioner = ModuleRI::SternheimerWeakSpectralPreconditioner;
using Grid = Preconditioner::Grid;
using Vector = Preconditioner::Vector;
using Complex = Preconditioner::Complex;
using Lattice = ModuleRI::SternheimerFDLatticeVectors;

std::size_t grid_count(const Grid& grid)
{
    if (!grid.periodic)
    {
        throw std::invalid_argument(
            "Weak spectral preconditioner requires a periodic grid.");
    }
    std::size_t count = 1;
    for (const int dimension : {grid.nx, grid.ny, grid.nz})
    {
        if (dimension <= 0
            || count > static_cast<std::size_t>(std::numeric_limits<int>::max())
                           / static_cast<std::size_t>(dimension))
        {
            throw std::invalid_argument(
                "Weak spectral preconditioner has invalid grid dimensions.");
        }
        count *= static_cast<std::size_t>(dimension);
    }
    for (const double spacing : {grid.hx, grid.hy, grid.hz})
    {
        if (!std::isfinite(spacing) || spacing <= 0.0)
        {
            throw std::invalid_argument(
                "Weak spectral preconditioner requires positive grid spacings.");
        }
    }
    for (const double component : grid.kpoint)
    {
        if (!std::isfinite(component))
        {
            throw std::invalid_argument(
                "Weak spectral preconditioner requires a finite Bloch point.");
        }
    }
    const Lattice lattice = ModuleRI::sternheimer_fd_grid_lattice_vectors(grid);
    for (const auto& vector : lattice)
    {
        for (const double component : vector)
        {
            if (!std::isfinite(component))
            {
                throw std::invalid_argument(
                    "Weak spectral preconditioner requires a finite lattice.");
            }
        }
    }
    const double determinant
        = lattice[0][0] * (lattice[1][1] * lattice[2][2]
                           - lattice[1][2] * lattice[2][1])
          - lattice[0][1] * (lattice[1][0] * lattice[2][2]
                             - lattice[1][2] * lattice[2][0])
          + lattice[0][2] * (lattice[1][0] * lattice[2][1]
                             - lattice[1][1] * lattice[2][0]);
    if (!std::isfinite(determinant) || determinant == 0.0)
    {
        throw std::invalid_argument(
            "Weak spectral preconditioner requires a nonsingular lattice.");
    }
    return count;
}

std::size_t grid_index(const Grid& grid, const int ix, const int iy, const int iz)
{
    return (static_cast<std::size_t>(ix) * grid.ny + iy) * grid.nz + iz;
}

Lattice reciprocal_vectors(const Grid& grid)
{
    Lattice reciprocal = ModuleRI::sternheimer_fd_grid_dual_vectors(grid);
    const double two_pi = 2.0 * std::acos(-1.0);
    for (auto& vector : reciprocal)
    {
        for (double& component : vector)
        {
            component *= two_pi;
        }
    }
    return reciprocal;
}

std::array<double, 3> wavevector(const Grid& grid,
                                 const Lattice& reciprocal,
                                 const int ix,
                                 const int iy,
                                 const int iz)
{
    const std::array<int, 3> dimensions{grid.nx, grid.ny, grid.nz};
    const std::array<int, 3> indices{ix, iy, iz};
    std::array<double, 3> result{};
    for (int direction = 0; direction != 3; ++direction)
    {
        const int frequency = indices[direction] <= dimensions[direction] / 2
                                  ? indices[direction]
                                  : indices[direction] - dimensions[direction];
        const double shifted = frequency + grid.kpoint[direction];
        for (int axis = 0; axis != 3; ++axis)
        {
            result[axis] += shifted * reciprocal[direction][axis];
        }
    }
    return result;
}

Vector bloch_phases(const Grid& grid, const std::size_t count)
{
    Vector phases(count);
    const double two_pi = 2.0 * std::acos(-1.0);
    const std::array<double, 3> k{
        std::remainder(grid.kpoint[0], static_cast<double>(grid.nx)),
        std::remainder(grid.kpoint[1], static_cast<double>(grid.ny)),
        std::remainder(grid.kpoint[2], static_cast<double>(grid.nz))};
    for (int ix = 0; ix != grid.nx; ++ix)
    {
        for (int iy = 0; iy != grid.ny; ++iy)
        {
            for (int iz = 0; iz != grid.nz; ++iz)
            {
                const double phase
                    = two_pi
                      * (k[0] * static_cast<double>(ix) / grid.nx
                         + k[1] * static_cast<double>(iy) / grid.ny
                         + k[2] * static_cast<double>(iz) / grid.nz);
                phases[grid_index(grid, ix, iy, iz)]
                    = std::exp(Complex(0.0, phase));
            }
        }
    }
    return phases;
}

bool finite(const Complex value)
{
    return std::isfinite(value.real()) && std::isfinite(value.imag());
}

void validate_input(const Vector& input, const std::size_t count)
{
    if (input.size() != count)
    {
        throw std::invalid_argument(
            "Weak spectral preconditioner input size does not match its grid.");
    }
    for (const Complex value : input)
    {
        if (!finite(value))
        {
            throw std::invalid_argument(
                "Weak spectral preconditioner requires finite input values.");
        }
    }
}

} // namespace

namespace ModuleRI
{

struct SternheimerWeakSpectralPreconditioner::Impl
{
    Impl(const Grid& input_grid,
         const double kinetic_prefactor,
         const double reference_eigenvalue,
         const double omega,
         const double regularization)
        : grid(input_grid), count(grid_count(grid)), phases(bloch_phases(grid, count)),
          inverse_denominator(count)
    {
        if (!std::isfinite(kinetic_prefactor) || kinetic_prefactor <= 0.0
            || !std::isfinite(reference_eigenvalue) || !std::isfinite(omega)
            || !std::isfinite(regularization) || regularization < 0.0)
        {
            throw std::invalid_argument(
                "Weak spectral preconditioner has an invalid kinetic shift.");
        }
        buffer = fftw_alloc_complex(count);
        if (buffer == nullptr)
        {
            throw std::bad_alloc();
        }
#ifdef _OPENMP
#pragma omp critical(sternheimer_fftw_plan_management)
#endif
        {
            forward = fftw_plan_dft_3d(grid.nx, grid.ny, grid.nz, buffer, buffer,
                                       FFTW_FORWARD, FFTW_ESTIMATE);
            backward = fftw_plan_dft_3d(grid.nx, grid.ny, grid.nz, buffer,
                                        buffer, FFTW_BACKWARD, FFTW_ESTIMATE);
        }
        if (forward == nullptr || backward == nullptr)
        {
            release();
            throw std::runtime_error(
                "Failed to initialize weak spectral preconditioner FFT plans.");
        }

        const Lattice reciprocal = reciprocal_vectors(grid);
        const double normalization = 1.0 / static_cast<double>(count);
        for (int ix = 0; ix != grid.nx; ++ix)
        {
            for (int iy = 0; iy != grid.ny; ++iy)
            {
                for (int iz = 0; iz != grid.nz; ++iz)
                {
                    const auto q = wavevector(grid, reciprocal, ix, iy, iz);
                    const double q2 = q[0] * q[0] + q[1] * q[1]
                                      + q[2] * q[2];
                    const Complex denominator(
                        kinetic_prefactor * q2 - reference_eigenvalue
                            + regularization,
                        omega);
                    const double magnitude = std::abs(denominator);
                    if (!std::isfinite(magnitude)
                        || magnitude <= std::numeric_limits<double>::min())
                    {
                        release();
                        throw std::domain_error(
                            "Weak spectral preconditioner found a singular Fourier mode.");
                    }
                    const Complex value = normalization / denominator;
                    if (!finite(value))
                    {
                        release();
                        throw std::overflow_error(
                            "Weak spectral preconditioner inverse symbol overflowed.");
                    }
                    inverse_denominator[grid_index(grid, ix, iy, iz)] = value;
                }
            }
        }
    }

    ~Impl()
    {
        release();
    }

    void release()
    {
#ifdef _OPENMP
#pragma omp critical(sternheimer_fftw_plan_management)
#endif
        {
            if (forward != nullptr)
            {
                fftw_destroy_plan(forward);
                forward = nullptr;
            }
            if (backward != nullptr)
            {
                fftw_destroy_plan(backward);
                backward = nullptr;
            }
        }
        fftw_free(buffer);
        buffer = nullptr;
    }

    void apply(const Vector& input, Vector& output) const
    {
        validate_input(input, count);
        for (std::size_t i = 0; i != count; ++i)
        {
            const Complex periodic = std::conj(phases[i]) * input[i];
            buffer[i][0] = periodic.real();
            buffer[i][1] = periodic.imag();
        }
        fftw_execute(forward);
        for (std::size_t i = 0; i != count; ++i)
        {
            const Complex coefficient(buffer[i][0], buffer[i][1]);
            const Complex value = inverse_denominator[i] * coefficient;
            buffer[i][0] = value.real();
            buffer[i][1] = value.imag();
        }
        fftw_execute(backward);
        Vector result(count);
        for (std::size_t i = 0; i != count; ++i)
        {
            result[i] = phases[i] * Complex(buffer[i][0], buffer[i][1]);
            if (!finite(result[i]))
            {
                throw std::overflow_error(
                    "Weak spectral preconditioner produced a nonfinite result.");
            }
        }
        output = std::move(result);
    }

    Grid grid;
    std::size_t count = 0;
    Vector phases;
    Vector inverse_denominator;
    mutable fftw_complex* buffer = nullptr;
    fftw_plan forward = nullptr;
    fftw_plan backward = nullptr;
};

SternheimerWeakSpectralPreconditioner::SternheimerWeakSpectralPreconditioner(
    const Grid& grid,
    const double kinetic_prefactor,
    const double reference_eigenvalue,
    const double omega,
    const double regularization)
    : impl_(new Impl(grid, kinetic_prefactor, reference_eigenvalue, omega,
                     regularization))
{}

SternheimerWeakSpectralPreconditioner::~SternheimerWeakSpectralPreconditioner()
    = default;

void SternheimerWeakSpectralPreconditioner::apply(const Vector& input,
                                                  Vector& output) const
{
    impl_->apply(input, output);
}

} // namespace ModuleRI
