#include "gtest/gtest.h"
#include "source_base/module_device/types.h"
#include "source_base/parallel_comm.h"
#include "source_estate/module_dm/edm_tddft.h"

#include <array>
#include <complex>
#include <vector>

namespace
{
using Complex = std::complex<double>;
using Matrix = std::array<Complex, 16>;

// The fixture matrices are row-major; distributed buffers use BLACS column-major storage.
void check_edm(const MPI_Comm comm,
               const Matrix& h_full,
               const Matrix& s_full,
               const Matrix& dmk_full,
               const Matrix& expected)
{
    Parallel_Blkcyclic layout;
    const int info = layout.init(4, 4, 1, comm, false);
    ASSERT_EQ(info, 0);
    const int nrow = layout.nrloc();
    const int ncol = layout.ncloc();
    std::vector<Complex> h(layout.size());
    std::vector<Complex> s(layout.size());
    std::vector<Complex> dmk(layout.size());
    for (int col = 0; col < ncol; ++col)
    {
        const int j = layout.icglb(col);
        for (int row = 0; row < nrow; ++row)
        {
            const int i = layout.irglb(row);
            const int index = col * nrow + row;
            const int global_index = i * 4 + j;
            h[index] = h_full[global_index];
            s[index] = s_full[global_index];
            dmk[index] = dmk_full[global_index];
        }
    }
    const std::vector<Complex> original_h = h;
    const std::vector<Complex> original_s = s;
    const std::vector<Complex> original_dmk = dmk;
    const std::size_t rows = nrow;
    const std::size_t cols = ncol;
    const ModuleBase::MatrixBlock<Complex> h_mat{h.data(), rows, cols, layout.desc()};
    const ModuleBase::MatrixBlock<Complex> s_mat{s.data(), rows, cols, layout.desc()};
    for (const bool tensor_lapack : {false, true})
    {
        ModuleBase::ComplexMatrix edmk(ncol, nrow, true);
        const char* kernel = tensor_lapack ? "Tensor/LAPACK" : "ScaLAPACK";
        SCOPED_TRACE(kernel);
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
        for (int col = 0; col < ncol; ++col)
        {
            const int j = layout.icglb(col);
            for (int row = 0; row < nrow; ++row)
            {
                const int i = layout.irglb(row);
                const Complex reference = expected[i * 4 + j];
                EXPECT_NEAR(edmk(col, row).real(), reference.real(), 1e-12);
                EXPECT_NEAR(edmk(col, row).imag(), reference.imag(), 1e-12);
            }
        }
        EXPECT_EQ(h, original_h);
        EXPECT_EQ(s, original_s);
        EXPECT_EQ(dmk, original_dmk);
    }
}

void check_diagonal_overlap(const MPI_Comm comm)
{
    Matrix h{};
    Matrix s{};
    Matrix dmk{};
    Matrix expected{};
    for (int i = 0; i < 4; ++i)
    {
        for (int j = 0; j < 4; ++j)
        {
            const int index = i * 4 + j;
            const double real = i + j + 1;
            const double imag = i - j;
            h[index] = Complex(real, imag);
            // S = 2I, D = diag(1,2,3,4): E_ij = H_ji * (D_ii + D_jj) / 4.
            expected[index] = std::conj(h[index]) * (i + j + 2.0) / 4.0;
            if (i == j)
            {
                s[index] = 2.0;
                dmk[index] = i + 1;
            }
        }
    }
    check_edm(comm, h, s, dmk, expected);
}

void check_noncommuting_overlap(const MPI_Comm comm)
{
    const Matrix h = {{
        {1, 0}, {2, -1}, {3, -2}, {4, -3},
        {2, 1}, {3, 0}, {4, -1}, {5, -2},
        {3, 2}, {4, 1}, {5, 0}, {6, -1},
        {4, 3}, {5, 2}, {6, 1}, {7, 0}
    }};
    // S = I + u*u^dagger, u = (1,i,0,0): eigenvalues 3,1,1,1.
    // S^-1 = I - u*u^dagger/3, and H*S^-1 != S^-1*H.
    const Matrix s = {{
        {2, 0}, {0, -1}, {0, 0}, {0, 0},
        {0, 1}, {2, 0}, {0, 0}, {0, 0},
        {0, 0}, {0, 0}, {1, 0}, {0, 0},
        {0, 0}, {0, 0}, {0, 0}, {1, 0}
    }};
    const Matrix dmk = {{
        {4, 0}, {0, 1}, {1, 0}, {0, 0},
        {0, -1}, {5, 0}, {0, 0}, {1, 0},
        {1, 0}, {0, 0}, {6, 0}, {0, -1},
        {0, 0}, {1, 0}, {0, 1}, {7, 0}
    }};
    // Independently expanded with the analytic inverse above:
    // 6E = (H*M)^T*D + D*(M*H)^T, M = 3*S^-1; T is not an adjoint.
    Matrix expected = {{
        {20, 0}, {56, 1}, {90, 44}, {134, 53},
        {56, -1}, {66, 0}, {116, 39}, {164, 70},
        {90, -44}, {116, -39}, {192, 0}, {258, 9},
        {134, -53}, {164, -70}, {258, -9}, {318, 0}
    }};
    for (Complex& value : expected)
    {
        value /= 6.0;
    }
    check_edm(comm, h, s, dmk, expected);
}
} // namespace

TEST(EdmTddft, ComplexMatricesInSuppliedLayout)
{
    check_diagonal_overlap(MPI_COMM_WORLD);
}

TEST(EdmTddft, NoncommutingOverlap)
{
    check_noncommuting_overlap(MPI_COMM_WORLD);
}

TEST(EdmTddft, SubcommunicatorRoots)
{
    MPICommGroup groups(MPI_COMM_WORLD);
    if (groups.gsize != 4)
    {
        GTEST_SKIP() << "The subgroup case uses the existing four-rank test registration.";
    }
    groups.divide_group_comm(2, true);
    if (groups.my_group == 1 && groups.rank_in_group == 0)
    {
        EXPECT_NE(groups.grank, 0);
    }
    // Each helper destroys its BLACS layout before groups frees the communicator.
    check_diagonal_overlap(groups.group_comm);
    check_noncommuting_overlap(groups.group_comm);
}
