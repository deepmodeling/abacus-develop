#include "source_lcao/module_ri/sternheimer_fd_preconditioner.h"

#include <cmath>
#include <complex>
#include <gtest/gtest.h>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

TEST(SternheimerFDPreconditioner, InvertsPeriodicFDKineticShiftOnBlochMode)
{
    using Hamiltonian = ModuleRI::SternheimerFDHamiltonian;
    using Complex = Hamiltonian::Complex;
    using Vector = Hamiltonian::Vector;

    Hamiltonian::Grid grid{8, 6, 4, 0.41, 0.53, 0.67, true};
    grid.kpoint = {0.17, -0.11, 0.08};
    Vector mode(static_cast<std::size_t>(grid.size()));
    const double pi = std::acos(-1.0);
    constexpr int mx = 2;
    constexpr int my = -1;
    constexpr int mz = 1;
    for (int ix = 0; ix != grid.nx; ++ix)
    {
        for (int iy = 0; iy != grid.ny; ++iy)
        {
            for (int iz = 0; iz != grid.nz; ++iz)
            {
                const double phase = 2.0 * pi
                                     * ((mx + grid.kpoint[0]) * static_cast<double>(ix) / grid.nx
                                        + (my + grid.kpoint[1]) * static_cast<double>(iy) / grid.ny
                                        + (mz + grid.kpoint[2]) * static_cast<double>(iz) / grid.nz);
                const int index = (ix * grid.ny + iy) * grid.nz + iz;
                mode[static_cast<std::size_t>(index)]
                    = std::exp(Complex(0.0, phase));
            }
        }
    }

    constexpr double reference_eigenvalue = -0.83;
    constexpr double omega = 0.37;
    constexpr double eta = 0.0;
    for (const int order : {2, 4, 6, 8})
    {
        Hamiltonian hamiltonian(grid, std::vector<double>(grid.size(), 0.0), 1.0, nullptr, order);
        ModuleRI::SternheimerFDSpectralPreconditioner preconditioner(
            hamiltonian, reference_eigenvalue, omega, eta);
        Vector preconditioned;
        preconditioner.apply(mode, preconditioned);

        Vector recovered;
        hamiltonian.apply_kinetic(preconditioned, recovered);
        const Complex shift(-reference_eigenvalue + eta, omega);
        for (std::size_t ir = 0; ir != recovered.size(); ++ir)
        {
            recovered[ir] += shift * preconditioned[ir];
            EXPECT_NEAR(recovered[ir].real(), mode[ir].real(), 2.0e-11) << "FD order " << order;
            EXPECT_NEAR(recovered[ir].imag(), mode[ir].imag(), 2.0e-11) << "FD order " << order;
        }
    }
}

TEST(SternheimerFDPreconditioner, RejectsNegativeRegularization)
{
    using Hamiltonian = ModuleRI::SternheimerFDHamiltonian;
    Hamiltonian::Grid grid{4, 4, 4, 0.5, 0.5, 0.5, true};
    Hamiltonian hamiltonian(grid, std::vector<double>(grid.size(), 0.0), 1.0, nullptr, 8);

    EXPECT_THROW(ModuleRI::SternheimerFDSpectralPreconditioner(hamiltonian, -0.5, 0.2, -0.1),
                 std::invalid_argument);
}

TEST(SternheimerFDPreconditioner, BatchMatchesIndependentScalarApplications)
{
    using Hamiltonian = ModuleRI::SternheimerFDHamiltonian;
    using Complex = Hamiltonian::Complex;
    using Vector = Hamiltonian::Vector;
    using Matrix = Hamiltonian::Matrix;
    Hamiltonian::Grid grid{8, 6, 4, 0.41, 0.53, 0.67, true};
    grid.kpoint = {0.17, -0.11, 0.08};
    const Hamiltonian hamiltonian(grid, std::vector<double>(grid.size(), 0.0), 1.0, nullptr, 8);
    const ModuleRI::SternheimerFDSpectralPreconditioner preconditioner(hamiltonian, -0.83, 0.37, 0.02);

    Matrix input(4, Vector(static_cast<std::size_t>(grid.size())));
    for (std::size_t column = 0; column != input.size(); ++column)
    {
        for (int ir = 0; ir != grid.size(); ++ir)
        {
            input[column][static_cast<std::size_t>(ir)]
                = Complex(0.01 * static_cast<double>((ir + 2 * column) % 17),
                          -0.02 * static_cast<double>((3 * ir + column) % 13));
        }
    }
    Matrix expected(input.size());
    for (std::size_t column = 0; column != input.size(); ++column)
    {
        preconditioner.apply(input[column], expected[column]);
    }

    Matrix actual;
    preconditioner.apply_batch(input, actual);

    EXPECT_EQ(actual, expected);
    preconditioner.apply_batch({}, actual);
    EXPECT_TRUE(actual.empty());
    input[1].pop_back();
    EXPECT_THROW(preconditioner.apply_batch(input, actual), std::invalid_argument);
}

TEST(SternheimerFDPreconditioner, ParallelizesIndependentBatchTransforms)
{
#ifdef _OPENMP
    using Hamiltonian = ModuleRI::SternheimerFDHamiltonian;
    using Complex = Hamiltonian::Complex;
    using Vector = Hamiltonian::Vector;
    using Matrix = Hamiltonian::Matrix;
    Hamiltonian::Grid grid{8, 6, 4, 0.41, 0.53, 0.67, true};
    grid.kpoint = {0.17, -0.11, 0.08};
    const Hamiltonian hamiltonian(grid, std::vector<double>(grid.size(), 0.0), 1.0, nullptr, 8);
    const ModuleRI::SternheimerFDSpectralPreconditioner preconditioner(hamiltonian, -0.83, 0.37, 0.02);
    Matrix input(4, Vector(static_cast<std::size_t>(grid.size()), Complex(0.25, -0.5)));
    Matrix output;

    const int previous_threads = omp_get_max_threads();
    const int previous_dynamic = omp_get_dynamic();
    omp_set_dynamic(0);
    omp_set_num_threads(4);
    int threads_used = 0;
    preconditioner.apply_batch(input, output, &threads_used);
    omp_set_num_threads(previous_threads);
    omp_set_dynamic(previous_dynamic);

    EXPECT_EQ(output.size(), input.size());
    EXPECT_GT(threads_used, 1);
#endif
}

TEST(SternheimerFDPreconditioner, ReusesCompatibleThreadLocalWorkspace)
{
    using Hamiltonian = ModuleRI::SternheimerFDHamiltonian;
    Hamiltonian::Grid grid{8, 6, 4, 0.5, 0.5, 0.5, true};
    Hamiltonian hamiltonian(grid, std::vector<double>(grid.size(), 0.0), 1.0, nullptr, 8);

    const auto first = ModuleRI::get_thread_local_sternheimer_fd_spectral_preconditioner(
        hamiltonian, -0.5, 0.2, 0.1);
    const auto same = ModuleRI::get_thread_local_sternheimer_fd_spectral_preconditioner(
        hamiltonian, -0.5, 0.2, 0.1);
    const auto shifted = ModuleRI::get_thread_local_sternheimer_fd_spectral_preconditioner(
        hamiltonian, -0.5, 0.3, 0.1);

    EXPECT_EQ(first.get(), same.get());
    EXPECT_NE(same.get(), shifted.get());
}
