#include "source_lcao/module_ri/sternheimer_fd_preconditioner.h"

#include "source_base/constants.h"

#include <fftw3.h>

#include <array>
#include <cmath>
#include <complex>
#include <new>
#include <stdexcept>
#include <vector>

namespace
{

struct FiniteDifferenceWeights
{
    int radius = 0;
    std::array<double, 5> second{};
    std::array<double, 4> first{};
};

FiniteDifferenceWeights finite_difference_weights(const int order)
{
    FiniteDifferenceWeights weights;
    weights.radius = order / 2;
    switch (order)
    {
    case 2:
        weights.second = {{-2.0, 1.0, 0.0, 0.0, 0.0}};
        weights.first = {{1.0 / 2.0, 0.0, 0.0, 0.0}};
        break;
    case 4:
        weights.second = {{-5.0 / 2.0, 4.0 / 3.0, -1.0 / 12.0, 0.0, 0.0}};
        weights.first = {{2.0 / 3.0, -1.0 / 12.0, 0.0, 0.0}};
        break;
    case 6:
        weights.second = {{-49.0 / 18.0, 3.0 / 2.0, -3.0 / 20.0, 1.0 / 90.0, 0.0}};
        weights.first = {{3.0 / 4.0, -3.0 / 20.0, 1.0 / 60.0, 0.0}};
        break;
    case 8:
        weights.second = {{-205.0 / 72.0, 8.0 / 5.0, -1.0 / 5.0, 8.0 / 315.0, -1.0 / 560.0}};
        weights.first = {{4.0 / 5.0, -1.0 / 5.0, 4.0 / 105.0, -1.0 / 280.0}};
        break;
    default:
        throw std::logic_error("Unsupported Sternheimer finite-difference order.");
    }
    return weights;
}

int signed_fft_index(const int index, const int dimension)
{
    return index <= dimension / 2 ? index : index - dimension;
}

struct KineticSymbolData
{
    std::array<std::array<double, 3>, 3> laplacian_coefficients{};
    std::array<std::vector<double>, 3> first_symbols;
    std::array<std::vector<double>, 3> second_symbols;
};

KineticSymbolData make_kinetic_symbol_data(const ModuleRI::SternheimerFDHamiltonian& hamiltonian)
{
    const auto& grid = hamiltonian.grid();
    const FiniteDifferenceWeights weights = finite_difference_weights(hamiltonian.finite_difference_order());
    const ModuleRI::SternheimerFDLatticeVectors dual = ModuleRI::sternheimer_fd_grid_dual_vectors(grid);
    const std::array<double, 3> dimensions{static_cast<double>(grid.nx),
                                           static_cast<double>(grid.ny),
                                           static_cast<double>(grid.nz)};
    KineticSymbolData data;
    for (int left = 0; left != 3; ++left)
    {
        for (int right = 0; right != 3; ++right)
        {
            for (int component = 0; component != 3; ++component)
            {
                data.laplacian_coefficients[left][right]
                    += dual[left][component] * dual[right][component]
                       * dimensions[left] * dimensions[right];
            }
        }
    }

    const std::array<int, 3> grid_dimensions{grid.nx, grid.ny, grid.nz};
    for (int direction = 0; direction != 3; ++direction)
    {
        data.first_symbols[direction].resize(static_cast<std::size_t>(grid_dimensions[direction]));
        data.second_symbols[direction].resize(static_cast<std::size_t>(grid_dimensions[direction]));
        for (int index = 0; index != grid_dimensions[direction]; ++index)
        {
            const double reduced_mode = static_cast<double>(signed_fft_index(index, grid_dimensions[direction]))
                                        + (grid.periodic ? grid.kpoint[direction] : 0.0);
            const double angle = ModuleBase::TWO_PI * reduced_mode / grid_dimensions[direction];
            double first_symbol = 0.0;
            double second_symbol = weights.second[0];
            for (int offset = 1; offset <= weights.radius; ++offset)
            {
                second_symbol += 2.0 * weights.second[static_cast<std::size_t>(offset)]
                                 * std::cos(offset * angle);
                first_symbol += 2.0 * weights.first[static_cast<std::size_t>(offset - 1)]
                                * std::sin(offset * angle);
            }
            data.first_symbols[direction][static_cast<std::size_t>(index)] = first_symbol;
            data.second_symbols[direction][static_cast<std::size_t>(index)] = second_symbol;
        }
    }
    return data;
}

double kinetic_symbol(const KineticSymbolData& data,
                      const ModuleRI::SternheimerFDHamiltonian& hamiltonian,
                      const int ix,
                      const int iy,
                      const int iz)
{
    const std::array<int, 3> indices{ix, iy, iz};

    double laplacian_symbol = 0.0;
    for (int direction = 0; direction != 3; ++direction)
    {
        laplacian_symbol
            += data.laplacian_coefficients[direction][direction]
               * data.second_symbols[direction][static_cast<std::size_t>(indices[direction])];
    }
    for (int left = 0; left != 3; ++left)
    {
        for (int right = left + 1; right != 3; ++right)
        {
            laplacian_symbol -= 2.0 * data.laplacian_coefficients[left][right]
                                * data.first_symbols[left][static_cast<std::size_t>(indices[left])]
                                * data.first_symbols[right][static_cast<std::size_t>(indices[right])];
        }
    }
    return -hamiltonian.kinetic_prefactor() * laplacian_symbol;
}

} // namespace

namespace ModuleRI
{

struct SternheimerFDSpectralPreconditioner::Impl
{
    explicit Impl(const SternheimerFDHamiltonian& hamiltonian,
                  const double reference_eigenvalue,
                  const double omega,
                  const double regularization)
        : grid(hamiltonian.grid()),
          kinetic_prefactor(hamiltonian.kinetic_prefactor()),
          finite_difference_order(hamiltonian.finite_difference_order()),
          reference_eigenvalue(reference_eigenvalue),
          omega(omega),
          regularization(regularization),
          inverse_denominator(static_cast<std::size_t>(grid.size()))
    {
        if (!std::isfinite(reference_eigenvalue) || !std::isfinite(omega)
            || !std::isfinite(regularization) || regularization < 0.0)
        {
            throw std::invalid_argument("Invalid Sternheimer FD spectral preconditioner shift.");
        }
        buffer = fftw_alloc_complex(static_cast<std::size_t>(grid.size()));
        if (buffer == nullptr)
        {
            throw std::bad_alloc();
        }
#ifdef _OPENMP
#pragma omp critical(sternheimer_fftw_plan_management)
#endif
        {
            forward = fftw_plan_dft_3d(grid.nx, grid.ny, grid.nz,
                                       buffer, buffer, FFTW_FORWARD, FFTW_ESTIMATE);
            backward = fftw_plan_dft_3d(grid.nx, grid.ny, grid.nz,
                                        buffer, buffer, FFTW_BACKWARD, FFTW_ESTIMATE);
        }
        if (forward == nullptr || backward == nullptr)
        {
            destroy();
            throw std::runtime_error("Failed to initialize the Sternheimer spectral preconditioner FFT.");
        }

        has_bloch_phase = grid.periodic
                          && (grid.kpoint[0] != 0.0 || grid.kpoint[1] != 0.0 || grid.kpoint[2] != 0.0);
        if (has_bloch_phase)
        {
            bloch_phase.resize(static_cast<std::size_t>(grid.size()));
            for (int ix = 0; ix != grid.nx; ++ix)
            {
                for (int iy = 0; iy != grid.ny; ++iy)
                {
                    for (int iz = 0; iz != grid.nz; ++iz)
                    {
                        const int index = (ix * grid.ny + iy) * grid.nz + iz;
                        const double phase_angle
                            = ModuleBase::TWO_PI
                              * (grid.kpoint[0] * static_cast<double>(ix) / grid.nx
                                 + grid.kpoint[1] * static_cast<double>(iy) / grid.ny
                                 + grid.kpoint[2] * static_cast<double>(iz) / grid.nz);
                        bloch_phase[static_cast<std::size_t>(index)]
                            = std::exp(std::complex<double>(0.0, phase_angle));
                    }
                }
            }
        }

        const KineticSymbolData symbol_data = make_kinetic_symbol_data(hamiltonian);
        const double normalization = 1.0 / static_cast<double>(grid.size());
        for (int ix = 0; ix != grid.nx; ++ix)
        {
            for (int iy = 0; iy != grid.ny; ++iy)
            {
                for (int iz = 0; iz != grid.nz; ++iz)
                {
                    const int index = (ix * grid.ny + iy) * grid.nz + iz;
                    const std::complex<double> denominator(
                        kinetic_symbol(symbol_data, hamiltonian, ix, iy, iz)
                            - reference_eigenvalue + regularization,
                        omega);
                    if (std::abs(denominator) < 1.0e-14)
                    {
                        destroy();
                        throw std::runtime_error("Sternheimer spectral preconditioner found a singular mode.");
                    }
                    inverse_denominator[static_cast<std::size_t>(index)] = normalization / denominator;
                }
            }
        }
    }

    ~Impl()
    {
        destroy();
    }

    void destroy()
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
        if (buffer != nullptr)
        {
            fftw_free(buffer);
            buffer = nullptr;
        }
    }

    SternheimerFDHamiltonian::Grid grid;
    double kinetic_prefactor = 0.0;
    int finite_difference_order = 0;
    double reference_eigenvalue = 0.0;
    double omega = 0.0;
    double regularization = 0.0;
    bool has_bloch_phase = false;
    std::vector<std::complex<double>> bloch_phase;
    fftw_complex* buffer = nullptr;
    fftw_plan forward = nullptr;
    fftw_plan backward = nullptr;
    std::vector<std::complex<double>> inverse_denominator;
};

SternheimerFDSpectralPreconditioner::SternheimerFDSpectralPreconditioner(
    const SternheimerFDHamiltonian& hamiltonian,
    const double reference_eigenvalue,
    const double omega,
    const double regularization)
    : impl_(new Impl(hamiltonian, reference_eigenvalue, omega, regularization))
{
}

SternheimerFDSpectralPreconditioner::~SternheimerFDSpectralPreconditioner() = default;

void SternheimerFDSpectralPreconditioner::apply(const SternheimerFDHamiltonian::Vector& input,
                                                 SternheimerFDHamiltonian::Vector& output) const
{
    if (input.size() != static_cast<std::size_t>(impl_->grid.size()))
    {
        throw std::invalid_argument("Sternheimer spectral preconditioner input size does not match the grid.");
    }

    for (int ix = 0; ix != impl_->grid.nx; ++ix)
    {
        for (int iy = 0; iy != impl_->grid.ny; ++iy)
        {
            for (int iz = 0; iz != impl_->grid.nz; ++iz)
            {
                const int index = (ix * impl_->grid.ny + iy) * impl_->grid.nz + iz;
                std::complex<double> value = input[static_cast<std::size_t>(index)];
                if (impl_->has_bloch_phase)
                {
                    value *= std::conj(impl_->bloch_phase[static_cast<std::size_t>(index)]);
                }
                impl_->buffer[index][0] = value.real();
                impl_->buffer[index][1] = value.imag();
            }
        }
    }
    fftw_execute(impl_->forward);
    for (int index = 0; index != impl_->grid.size(); ++index)
    {
        const std::complex<double> value(impl_->buffer[index][0], impl_->buffer[index][1]);
        const std::complex<double> preconditioned
            = impl_->inverse_denominator[static_cast<std::size_t>(index)] * value;
        impl_->buffer[index][0] = preconditioned.real();
        impl_->buffer[index][1] = preconditioned.imag();
    }
    fftw_execute(impl_->backward);

    output.resize(input.size());
    for (int ix = 0; ix != impl_->grid.nx; ++ix)
    {
        for (int iy = 0; iy != impl_->grid.ny; ++iy)
        {
            for (int iz = 0; iz != impl_->grid.nz; ++iz)
            {
                const int index = (ix * impl_->grid.ny + iy) * impl_->grid.nz + iz;
                const std::complex<double> value(impl_->buffer[index][0], impl_->buffer[index][1]);
                output[static_cast<std::size_t>(index)]
                    = impl_->has_bloch_phase
                          ? impl_->bloch_phase[static_cast<std::size_t>(index)] * value
                          : value;
            }
        }
    }
}

bool SternheimerFDSpectralPreconditioner::is_compatible(
    const SternheimerFDHamiltonian& hamiltonian,
    const double reference_eigenvalue,
    const double omega,
    const double regularization) const
{
    const auto& grid = hamiltonian.grid();
    return impl_->grid.nx == grid.nx && impl_->grid.ny == grid.ny && impl_->grid.nz == grid.nz
           && impl_->grid.hx == grid.hx && impl_->grid.hy == grid.hy && impl_->grid.hz == grid.hz
           && impl_->grid.periodic == grid.periodic && impl_->grid.kpoint == grid.kpoint
           && impl_->grid.lattice_vectors == grid.lattice_vectors
           && impl_->kinetic_prefactor == hamiltonian.kinetic_prefactor()
           && impl_->finite_difference_order == hamiltonian.finite_difference_order()
           && impl_->reference_eigenvalue == reference_eigenvalue && impl_->omega == omega
           && impl_->regularization == regularization;
}

std::shared_ptr<SternheimerFDSpectralPreconditioner>
get_thread_local_sternheimer_fd_spectral_preconditioner(const SternheimerFDHamiltonian& hamiltonian,
                                                         const double reference_eigenvalue,
                                                         const double omega,
                                                         const double regularization)
{
    static thread_local std::shared_ptr<SternheimerFDSpectralPreconditioner> cached;
    if (!cached || !cached->is_compatible(hamiltonian, reference_eigenvalue, omega, regularization))
    {
        cached = std::make_shared<SternheimerFDSpectralPreconditioner>(
            hamiltonian, reference_eigenvalue, omega, regularization);
    }
    return cached;
}

} // namespace ModuleRI
