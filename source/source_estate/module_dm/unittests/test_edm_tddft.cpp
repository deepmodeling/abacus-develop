#include "gtest/gtest.h"
#include "source_base/module_device/types.h"
#include "source_estate/module_dm/edm_tddft.h"

#include <complex>
#include <vector>

TEST(EdmTddft, ComplexMatricesInSuppliedLayout)
{
    Parallel_Blkcyclic layout;
    const int info = layout.init(4, 4, 1, MPI_COMM_WORLD, false);
    ASSERT_EQ(info, 0);
    const int nrow = layout.nrloc();
    const int ncol = layout.ncloc();
    std::vector<std::complex<double>> h(layout.size());
    std::vector<std::complex<double>> s(layout.size());
    std::vector<std::complex<double>> dmk(layout.size());
    for (int col = 0; col < ncol; ++col)
    {
        const int j = layout.icglb(col);
        for (int row = 0; row < nrow; ++row)
        {
            const int i = layout.irglb(row);
            const int index = col * nrow + row;
            const double real = i + j + 1;
            const double imag = i - j;
            h[index] = std::complex<double>(real, imag);
            if (i == j)
            {
                s[index] = 2.0;
                dmk[index] = i + 1;
            }
        }
    }
    const auto original_h = h;
    const auto original_s = s;
    const auto original_dmk = dmk;
    const std::size_t rows = nrow;
    const std::size_t cols = ncol;
    const ModuleBase::MatrixBlock<std::complex<double>> h_mat{h.data(), rows, cols, layout.desc()};
    const ModuleBase::MatrixBlock<std::complex<double>> s_mat{s.data(), rows, cols, layout.desc()};
    ModuleBase::ComplexMatrix edmk;
    for (const bool tensor_lapack : {false, true})
    {
        if (tensor_lapack)
        {
            module_dm::edm_tddft_lapack<base_device::DEVICE_CPU>(layout, h_mat, s_mat, dmk.data(), edmk);
        }
        else
        {
            module_dm::edm_tddft(layout, h_mat, s_mat, dmk.data(), edmk);
        }
        EXPECT_EQ(edmk.nr, ncol);
        EXPECT_EQ(edmk.nc, nrow);
        // S = 2I, D = diag(1,2,3,4): E_ij = H_ji * (D_ii + D_jj) / 4.
        // Complex off-diagonal H distinguishes the required transpose from an adjoint.
        for (int col = 0; col < ncol; ++col)
        {
            const int j = layout.icglb(col);
            for (int row = 0; row < nrow; ++row)
            {
                const int i = layout.irglb(row);
                const double real = i + j + 1;
                const double imag = j - i;
                const std::complex<double> expected = std::complex<double>(real, imag) * (i + j + 2.0) / 4.0;
                EXPECT_NEAR(edmk(col, row).real(), expected.real(), 1e-12);
                EXPECT_NEAR(edmk(col, row).imag(), expected.imag(), 1e-12);
            }
        }
        EXPECT_EQ(h, original_h);
        EXPECT_EQ(s, original_s);
        EXPECT_EQ(dmk, original_dmk);
    }
}
