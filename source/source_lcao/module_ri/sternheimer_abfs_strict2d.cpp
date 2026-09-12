#include "source_lcao/module_ri/sternheimer_abfs_strict2d.h"

#include "source_base/constants.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <fftw3.h>
#include <limits>
#include <new>
#include <stdexcept>

namespace
{
using Complex = std::complex<double>;
using Grid = ModuleRI::SternheimerFDHamiltonian::Grid;

struct Geometry
{
    std::size_t size;
    std::size_t nxy;
    double dz;
    double dv;
    ModuleRI::SternheimerFDLatticeVectors dual;
};

Geometry checked_geometry(const Grid& grid)
{
    if (!grid.periodic || grid.nx <= 0 || grid.ny <= 0 || grid.nz <= 0)
    {
        throw std::invalid_argument("Strict2D ABFS Coulomb requires positive dimensions and periodic xy sampling.");
    }
    for (const double spacing: {grid.hx, grid.hy, grid.hz})
    {
        if (!std::isfinite(spacing) || spacing <= 0.0)
        {
            throw std::invalid_argument("Strict2D ABFS Coulomb requires finite positive grid spacings.");
        }
    }
    // FFTW's many-transform strides/dimensions and the existing Grid use int.
    std::size_t size = 1;
    for (const int dimension: {grid.nx, grid.ny, grid.nz})
    {
        if (size > static_cast<std::size_t>(std::numeric_limits<int>::max()) / dimension)
        {
            throw std::invalid_argument("Strict2D ABFS Coulomb grid exceeds the supported FFT index range.");
        }
        size *= dimension;
    }
    const auto lattice = ModuleRI::sternheimer_fd_grid_lattice_vectors(grid);
    std::array<double, 3> lengths{};
    for (int i = 0; i < 3; ++i)
    {
        for (const double component: lattice[i])
        {
            if (!std::isfinite(component))
            {
                throw std::invalid_argument("Strict2D ABFS Coulomb requires a finite lattice.");
            }
        }
        lengths[i] = std::hypot(std::hypot(lattice[i][0], lattice[i][1]), lattice[i][2]);
        if (!std::isfinite(lengths[i]) || lengths[i] <= 0.0)
        {
            throw std::invalid_argument("Strict2D ABFS Coulomb requires nonzero finite lattice vectors.");
        }
    }
    constexpr double tolerance = 1.0e-12;
    if (std::abs(lattice[0][2]) > tolerance * lengths[0]
        || std::abs(lattice[1][2]) > tolerance * lengths[1]
        || std::hypot(lattice[2][0], lattice[2][1]) > tolerance * lengths[2])
    {
        throw std::invalid_argument("Strict2D ABFS Coulomb does not support tilted cells: require xy plane and z normal.");
    }
    const double area = std::abs(lattice[0][0] * lattice[1][1] - lattice[0][1] * lattice[1][0]);
    const double volume = area * std::abs(lattice[2][2]);
    if (!std::isfinite(area) || area <= tolerance * lengths[0] * lengths[1]
        || !std::isfinite(volume) || volume <= 0.0)
    {
        throw std::invalid_argument("Strict2D ABFS Coulomb requires a nonsingular planar cell.");
    }
    // Drop only accepted geometry roundoff, so Q has no out-of-plane part.
    Grid planar_grid = grid;
    planar_grid.lattice_vectors = lattice;
    planar_grid.lattice_vectors[0][2] = 0.0;
    planar_grid.lattice_vectors[1][2] = 0.0;
    planar_grid.lattice_vectors[2][0] = 0.0;
    planar_grid.lattice_vectors[2][1] = 0.0;
    return {size,
            static_cast<std::size_t>(grid.nx) * grid.ny,
            std::abs(lattice[2][2]) / grid.nz,
            volume / static_cast<double>(size),
            ModuleRI::sternheimer_fd_grid_dual_vectors(planar_grid)};
}

bool finite(const Complex& value)
{
    return std::isfinite(value.real()) && std::isfinite(value.imag());
}

struct PlanarFFT
{
    fftw_complex* data = nullptr;
    fftw_plan forward = nullptr;
    fftw_plan backward = nullptr;

    PlanarFFT(const Grid& grid, const std::size_t size)
    {
        data = fftw_alloc_complex(size);
        if (data == nullptr)
        {
            throw std::bad_alloc();
        }
        // Original storage is [(x*ny+y)*nz+z]: nz independent strided xy FFTs.
        const int dimensions[2] = {grid.nx, grid.ny};
#ifdef _OPENMP
#pragma omp critical(sternheimer_fftw_plan_management)
#endif
        {
            forward = fftw_plan_many_dft(2, dimensions, grid.nz, data, nullptr, grid.nz, 1,
                                         data, nullptr, grid.nz, 1, FFTW_FORWARD, FFTW_ESTIMATE);
            backward = fftw_plan_many_dft(2, dimensions, grid.nz, data, nullptr, grid.nz, 1,
                                          data, nullptr, grid.nz, 1, FFTW_BACKWARD, FFTW_ESTIMATE);
        }
        if (forward == nullptr || backward == nullptr)
        {
            release();
            throw std::runtime_error("Failed to initialize strict2D ABFS planar FFTW plans.");
        }
    }

    PlanarFFT(const PlanarFFT&) = delete;
    PlanarFFT& operator=(const PlanarFFT&) = delete;

    ~PlanarFFT()
    {
        release();
    }

    void release()
    {
#ifdef _OPENMP
#pragma omp critical(sternheimer_fftw_plan_management)
#endif
        {
            if (forward != nullptr) fftw_destroy_plan(forward);
            if (backward != nullptr) fftw_destroy_plan(backward);
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
void solve_sternheimer_abf_strict2d_coulomb_in_place(
    std::vector<SternheimerABFBlochGridChannel>& density_channels,
    const SternheimerFDHamiltonian::Grid& grid,
    const SternheimerReducedKPoint& qpoint)
{
    const Geometry geometry = checked_geometry(grid);
    for (const double coordinate: qpoint)
    {
        if (!std::isfinite(coordinate))
        {
            throw std::invalid_argument("Strict2D ABFS Coulomb requires a finite reduced q point.");
        }
    }
    if (qpoint[2] != 0.0)
    {
        throw std::invalid_argument("Strict2D ABFS Coulomb supports only qz=0, without reciprocal-z folding.");
    }
    if (std::abs(std::remainder(qpoint[0], 1.0)) <= 1.0e-14
        && std::abs(std::remainder(qpoint[1], 1.0)) <= 1.0e-14)
    {
        throw std::invalid_argument("Strict2D ABFS Coulomb does not support Gamma or its planar reciprocal equivalents.");
    }
    // Validate every channel before changing the first one.
    for (const auto& channel: density_channels)
    {
        if (channel.potential_r.size() != geometry.size)
        {
            throw std::invalid_argument("Strict2D ABFS Coulomb density size does not match the fine grid.");
        }
        for (const Complex value: channel.potential_r)
        {
            if (!finite(value))
            {
                throw std::invalid_argument("Strict2D ABFS Coulomb density contains a nonfinite value.");
            }
        }
    }
    std::vector<Complex> phase(geometry.nxy);
    std::vector<double> decay(geometry.nxy);
    std::vector<double> factor(geometry.nxy);
    for (int x = 0; x < grid.nx; ++x)
    {
        for (int y = 0; y < grid.ny; ++y)
        {
            const std::size_t p = static_cast<std::size_t>(x) * grid.ny + y;
            const int gx = x < grid.nx / 2 + grid.nx % 2 ? x : x - grid.nx;
            const int gy = y < grid.ny / 2 + grid.ny % 2 ? y : y - grid.ny;
            const double Qx = ModuleBase::TWO_PI * ((gx + qpoint[0]) * geometry.dual[0][0]
                                                  + (gy + qpoint[1]) * geometry.dual[1][0]);
            const double Qy = ModuleBase::TWO_PI * ((gx + qpoint[0]) * geometry.dual[0][1]
                                                  + (gy + qpoint[1]) * geometry.dual[1][1]);
            const double Q = std::hypot(Qx, Qy);
            if (!std::isfinite(Q) || Q <= 0.0)
            {
                throw std::invalid_argument("Strict2D ABFS Coulomb requires finite nonzero Q for EVERY planar mode.");
            }
            const double angle = ModuleBase::TWO_PI * (qpoint[0] * (static_cast<double>(x) / grid.nx)
                                                       + qpoint[1] * (static_cast<double>(y) / grid.ny));
            phase[p] = std::exp(Complex(0.0, angle));
            decay[p] = std::exp(-Q * geometry.dz);
            factor[p] = (ModuleBase::TWO_PI / Q) * geometry.dz / static_cast<double>(geometry.nxy);
            if (!finite(phase[p]) || !std::isfinite(factor[p]) || factor[p] <= 0.0)
            {
                throw std::invalid_argument("Strict2D ABFS Coulomb has an unrepresentable phase or kernel scale.");
            }
        }
    }
    if (density_channels.empty()) return;

    PlanarFFT fft(grid, geometry.size);
    std::vector<Complex> left(grid.nz);
    for (auto& channel: density_channels)
    {
        for (std::size_t p = 0; p < geometry.nxy; ++p)
        {
            for (int z = 0; z < grid.nz; ++z)
            {
                const std::size_t i = p * grid.nz + z;
                const Complex value = std::conj(phase[p]) * channel.potential_r[i];
                if (!finite(value))
                {
                    throw std::overflow_error("Strict2D ABFS Coulomb overflow in density Bloch phase removal.");
                }
                fft.data[i][0] = value.real();
                fft.data[i][1] = value.imag();
            }
        }
        fftw_execute(fft.forward);
        for (std::size_t p = 0; p < geometry.nxy; ++p)
        {
            const std::size_t base = p * grid.nz;
            const double a = decay[p];
            // left[z] includes the diagonal. right excludes it at evaluation;
            // this counts each source once, without periodic wrap or subtraction.
            Complex sum{};
            for (int z = 0; z < grid.nz; ++z)
            {
                const Complex rho(fft.data[base + z][0], fft.data[base + z][1]);
                if (!finite(rho))
                {
                    throw std::overflow_error("Strict2D ABFS Coulomb overflow in forward FFT.");
                }
                sum = rho + a * sum;
                if (!finite(sum))
                {
                    throw std::overflow_error("Strict2D ABFS Coulomb overflow in forward recurrence.");
                }
                left[z] = sum;
            }
            Complex right{};
            for (int z = grid.nz - 1; z >= 0; --z)
            {
                const Complex rho(fft.data[base + z][0], fft.data[base + z][1]);
                const Complex convolution = left[z] + a * right;
                const Complex potential = factor[p] * convolution;
                right = rho + a * right;
                if (!finite(convolution) || !finite(right) || !finite(potential))
                {
                    throw std::overflow_error("Strict2D ABFS Coulomb overflow in backward recurrence or kernel scaling.");
                }
                fft.data[base + z][0] = potential.real();
                fft.data[base + z][1] = potential.imag();
            }
        }
        fftw_execute(fft.backward);
        double max_abs = 0.0;
        for (std::size_t p = 0; p < geometry.nxy; ++p)
        {
            for (int z = 0; z < grid.nz; ++z)
            {
                const std::size_t i = p * grid.nz + z;
                const Complex periodic_value(fft.data[i][0], fft.data[i][1]);
                if (!finite(periodic_value))
                {
                    throw std::overflow_error("Strict2D ABFS Coulomb overflow in inverse FFT.");
                }
                const Complex value = phase[p] * periodic_value;
                const double magnitude = std::abs(value);
                if (!finite(value) || !std::isfinite(magnitude))
                {
                    throw std::overflow_error("Strict2D ABFS Coulomb overflow in potential Bloch phase or magnitude.");
                }
                fft.data[i][0] = value.real();
                fft.data[i][1] = value.imag();
                max_abs = std::max(max_abs, magnitude);
            }
        }
        // Validate the complete channel in the workspace before replacing its
        // density. Earlier completed channels remain potentials on failure.
        for (std::size_t i = 0; i < geometry.size; ++i)
        {
            channel.potential_r[i] = Complex(fft.data[i][0], fft.data[i][1]);
        }
        channel.max_abs = max_abs;
    }
}

std::vector<std::complex<double>> sternheimer_abf_strict2d_selected_coulomb_integrals(
    const std::vector<SternheimerABFBlochGridChannel>& densities,
    const std::vector<SternheimerABFBlochGridChannel>& potentials,
    const SternheimerFDHamiltonian::Grid& grid,
    const std::vector<std::pair<std::size_t, std::size_t>>& pairs)
{
    const Geometry geometry = checked_geometry(grid);
    std::vector<Complex> result;
    result.reserve(pairs.size());
    for (const auto& pair: pairs)
    {
        if (pair.first >= densities.size() || pair.second >= potentials.size())
        {
            throw std::invalid_argument("Strict2D ABFS selected Coulomb integral has an invalid channel index.");
        }
        const auto& rho = densities[pair.first].potential_r;
        const auto& phi = potentials[pair.second].potential_r;
        if (rho.size() != geometry.size || phi.size() != geometry.size)
        {
            throw std::invalid_argument("Strict2D ABFS selected Coulomb integral size does not match the fine grid.");
        }
        Complex sum{};
        for (std::size_t i = 0; i < geometry.size; ++i)
        {
            if (!finite(rho[i]) || !finite(phi[i]))
            {
                throw std::invalid_argument("Strict2D ABFS selected Coulomb integral contains a nonfinite field.");
            }
            const Complex term = std::conj(rho[i]) * phi[i];
            if (!finite(term))
            {
                throw std::overflow_error("Strict2D ABFS selected Coulomb integral overflow in field product.");
            }
            sum += term;
            if (!finite(sum))
            {
                throw std::overflow_error("Strict2D ABFS selected Coulomb integral overflow in accumulation.");
            }
        }
        const Complex integral = geometry.dv * sum;
        if (!finite(integral))
        {
            throw std::overflow_error("Strict2D ABFS selected Coulomb integral overflow in volume scaling.");
        }
        result.push_back(integral);
    }
    return result;
}
} // namespace ModuleRI
