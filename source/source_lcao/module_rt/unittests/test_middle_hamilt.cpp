#include "source_base/global_variable.h"
#include "source_base/matrix2d_view.h"
#include "source_base/parallel_blkcyc.h"
#include "source_lcao/module_rt/middle_hamilt.h"

#include <gtest/gtest.h>

#include <complex>
#include <vector>

extern int ictxt;

namespace
{
constexpr double threshold = 1.0e-8;
} // namespace

TEST(MiddleHamiltTest, testMiddleHamilt)
{
    using Complex = std::complex<double>;
    const int nlocal = 4;
    Parallel_Blkcyclic layout;
    ASSERT_EQ(layout.set(nlocal, nlocal, 1, ictxt), 0);

    std::vector<Complex> current_hamiltonian(nlocal * nlocal);
    std::vector<Complex> previous_hamiltonian(nlocal * nlocal);
    std::vector<Complex> current_overlap(nlocal * nlocal);
    std::vector<Complex> previous_overlap(nlocal * nlocal);
    Matrix2dView<Complex> hamiltonian(
        current_hamiltonian.data(), nlocal, nlocal, layout.desc());
    Matrix2dView<Complex> overlap(
        current_overlap.data(), nlocal, nlocal, layout.desc());
    const Matrix2dView<const Complex> previous_hamiltonian_view(
        previous_hamiltonian.data(), nlocal, nlocal, layout.desc());
    const Matrix2dView<const Complex> previous_overlap_view(
        previous_overlap.data(), nlocal, nlocal, layout.desc());

    for (int i = 0; i < nlocal; ++i)
    {
        for (int j = 0; j < nlocal; ++j)
        {
            const double index = static_cast<double>(i * nlocal + j);
            if (i == j)
            {
                hamiltonian.at_local(i, j) = Complex(1.0, 0.0);
                overlap.at_local(i, j) = Complex(1.0, 0.0);
                previous_hamiltonian[i * nlocal + j]
                    = Complex(1.0 + 0.2 * index, 0.0);
                previous_overlap[i * nlocal + j]
                    = Complex(1.0 + 0.2 * index, 0.0);
            }
            else
            {
                previous_hamiltonian[i * nlocal + j]
                    = Complex(0.2 * index, 0.0);
                previous_overlap[i * nlocal + j]
                    = Complex(0.2 * index, 0.0);
            }
        }
    }

    const bool print_matrix = false;
    module_rt::half_Hmatrix(previous_hamiltonian_view,
                            previous_overlap_view,
                            hamiltonian,
                            overlap,
                            GlobalV::ofs_running,
                            print_matrix);

    for (int i = 0; i < nlocal; ++i)
    {
        for (int j = 0; j < nlocal; ++j)
        {
            const double index = static_cast<double>(i * nlocal + j);
            const double previous
                = i == j ? 1.0 + 0.2 * index : 0.2 * index;
            const double initial = i == j ? 1.0 : 0.0;
            const double expected = 0.5 * initial + 0.5 * previous;
            EXPECT_NEAR(hamiltonian.at_local(i, j).real(),
                        expected,
                        threshold);
            EXPECT_NEAR(hamiltonian.at_local(i, j).imag(), 0.0, threshold);
            EXPECT_NEAR(overlap.at_local(i, j).real(), expected, threshold);
            EXPECT_NEAR(overlap.at_local(i, j).imag(), 0.0, threshold);
        }
    }
}
