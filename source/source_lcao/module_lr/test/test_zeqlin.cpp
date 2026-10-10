#include <gtest/gtest.h>
#include "../zeq_solver.hpp"
#include "../cal_edm.h"
#include "../cal_edm.h"
#include "mpi.h"
#include "../zeqlin_solv.h"

#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>

#include "source_lcao/module_lr/utils/lr_util.h"

// The distributed Z-vector solvers are checked against the replicated LAPACK solve
// (`LR_Util::lapack_linear_solver`) on matrices generated identically on every rank.
namespace
{
    // deterministic pseudo-random entry in [-1, 1], identical on every rank
    double entry(const int i, const int j, const int salt)
    {
        return std::sin(0.7 * i + 1.3 * j + 0.37 * salt + 0.1 * i * j);
    }
    double conj_if_complex(const double v) { return v; }
    std::complex<double> conj_if_complex(const std::complex<double>& v) { return std::conj(v); }
    template <typename T> T make_value(const double re, const double im);
    template <> double make_value<double>(const double re, const double im) { return re; }
    template <> std::complex<double> make_value<std::complex<double>>(const double re, const double im)
    {
        return std::complex<double>(re, im);
    }

    /// column-major n x n Hermitian matrix; positive definite when `shift` > n
    template <typename T>
    std::vector<T> hermitian_matrix(const int n, const double shift)
    {
        std::vector<T> a(static_cast<std::size_t>(n) * n);
        for (int j = 0; j < n; ++j)
        {
            for (int i = 0; i <= j; ++i)
            {
                const T v = make_value<T>(entry(i, j, 0), entry(i, j, 1));
                a[j * n + i] = v;
                a[i * n + j] = conj_if_complex(v);
            }
            a[j * n + j] = T(entry(j, j, 2) + shift);
        }
        return a;
    }

    template <typename T>
    std::vector<T> rhs_matrix(const int n, const int nrhs)
    {
        std::vector<T> b(static_cast<std::size_t>(n) * nrhs);
        for (int j = 0; j < nrhs; ++j)
        {
            for (int i = 0; i < n; ++i) { b[j * n + i] = make_value<T>(entry(i, j, 3), entry(i, j, 4)); }
        }
        return b;
    }

    template <typename T>
    std::vector<T> to_local(const std::vector<T>& full, const Parallel_2D& pv)
    {
        const int nrow = pv.get_global_row_size();
        std::vector<T> loc(pv.get_local_size());
        for (int lc = 0; lc < pv.get_col_size(); ++lc)
        {
            for (int lr = 0; lr < pv.get_row_size(); ++lr)
            {
                loc[lc * pv.get_row_size() + lr] = full[pv.local2global_col(lc) * nrow + pv.local2global_row(lr)];
            }
        }
        return loc;
    }

    void expect_near(const double a, const double b) { EXPECT_NEAR(a, b, 1e-10); }
    void expect_near(const std::complex<double>& a, const std::complex<double>& b)
    {
        EXPECT_NEAR(a.real(), b.real(), 1e-10);
        EXPECT_NEAR(a.imag(), b.imag(), 1e-10);
    }

    /// solve the distributed system with `solver` and compare with LAPACK on `a_ref`
    template <typename T>
    void check_against_lapack(void (*solver)(T*, T*, const Parallel_2D&, const Parallel_2D&),
        const std::vector<T>& a_in, const std::vector<T>& a_ref, const int n, const int nrhs, const int nb)
    {
        const std::vector<T> b = rhs_matrix<T>(n, nrhs);
        std::vector<T> x_ref(b.size());
        LR_Util::lapack_linear_solver(a_ref.data(), x_ref.data(), b.data(), n, nrhs);

        Parallel_2D pa;
        LR_Util::setup_2d_division(pa, nb, n, n);
        Parallel_2D pb;
        LR_Util::setup_2d_division(pb, nb, n, nrhs, pa.blacs_ctxt);
        std::vector<T> a_loc = to_local(a_in, pa);
        std::vector<T> x_loc = to_local(b, pb);
        solver(a_loc.data(), x_loc.data(), pa, pb);

        std::vector<T> x(b.size(), T(0));
        LR_Util::gather_2d_to_full(pb, x_loc.data(), x.data(), false, n, nrhs);
        for (std::size_t i = 0; i < x.size(); ++i) { expect_near(x[i], x_ref[i]); }
    }
}

class ZeqLinearSolverTest : public testing::Test
{
  public:
    // (n, nrhs, nb): n not a multiple of nb, a single rhs, and nb = 1
    const std::vector<std::vector<int>> sizes{ {37, 3, 4}, {20, 1, 3}, {9, 2, 1} };
};

TEST_F(ZeqLinearSolverTest, ScalapackDouble)
{
    for (const auto& s : sizes)
    {
        const std::vector<double> a = hermitian_matrix<double>(s[0], s[0] + 1.0);
        check_against_lapack<double>(&LR::scalapack_linear_solver<double>, a, a, s[0], s[1], s[2]);
    }
}

TEST_F(ZeqLinearSolverTest, ScalapackComplex)
{
    typedef std::complex<double> C;
    for (const auto& s : sizes)
    {
        const std::vector<C> a = hermitian_matrix<C>(s[0], s[0] + 1.0);
        check_against_lapack<C>(&LR::scalapack_linear_solver<C>, a, a, s[0], s[1], s[2]);
    }
}

TEST_F(ZeqLinearSolverTest, ScalapackNonSymmetric)
{
    // LU must not assume symmetry: perturb with an antisymmetric part
    for (const auto& s : sizes)
    {
        const int n = s[0];
        std::vector<double> a = hermitian_matrix<double>(n, n + 1.0);
        for (int j = 0; j < n; ++j)
        {
            for (int i = 0; i < j; ++i)
            {
                a[j * n + i] += 0.3 * entry(i, j, 5);
                a[i * n + j] -= 0.3 * entry(i, j, 5);
            }
        }
        check_against_lapack<double>(&LR::scalapack_linear_solver<double>, a, a, n, s[1], s[2]);
    }
}

TEST_F(ZeqLinearSolverTest, ScalapackCholDouble)
{
    for (const auto& s : sizes)
    {
        const std::vector<double> a = hermitian_matrix<double>(s[0], s[0] + 1.0);
        check_against_lapack<double>(&LR::scalapack_cholesky_linear_solver<double>, a, a, s[0], s[1], s[2]);
    }
}

TEST_F(ZeqLinearSolverTest, ScalapackCholComplex)
{
    typedef std::complex<double> C;
    for (const auto& s : sizes)
    {
        const std::vector<C> a = hermitian_matrix<C>(s[0], s[0] + 1.0);
        check_against_lapack<C>(&LR::scalapack_cholesky_linear_solver<C>, a, a, s[0], s[1], s[2]);
    }
}

TEST_F(ZeqLinearSolverTest, ScalapackCholSolvesHermitianPart)
{
    // as for ELPA: a slightly non-symmetric input is solved as (A + A^T)/2
    for (const auto& s : sizes)
    {
        const int n = s[0];
        const std::vector<double> a_sym = hermitian_matrix<double>(n, n + 1.0);
        std::vector<double> a = a_sym;
        for (int j = 0; j < n; ++j)
        {
            for (int i = 0; i < j; ++i)
            {
                a[j * n + i] += 1e-2 * entry(i, j, 5);
                a[i * n + j] -= 1e-2 * entry(i, j, 5);
            }
        }
        check_against_lapack<double>(&LR::scalapack_cholesky_linear_solver<double>, a, a_sym, n, s[1], s[2]);
    }
}

TEST_F(ZeqLinearSolverTest, ScalapackCholRejectsIndefinite)
{
    // every rank: p?potrf's INFO is global, so unlike ELPA this must throw everywhere, not hang
    const int n = 20;
    const int nb = 3;
    const std::vector<double> a = hermitian_matrix<double>(n, -(n + 1.0));   // negative definite
    Parallel_2D pa;
    LR_Util::setup_2d_division(pa, nb, n, n);
    Parallel_2D pb;
    LR_Util::setup_2d_division(pb, nb, n, 1, pa.blacs_ctxt);
    std::vector<double> a_loc = to_local(a, pa);
    std::vector<double> b_loc = to_local(rhs_matrix<double>(n, 1), pb);
    EXPECT_THROW(LR::scalapack_cholesky_linear_solver<double>(a_loc.data(), b_loc.data(), pa, pb),
                 std::runtime_error);
}

#ifdef __ELPA
TEST_F(ZeqLinearSolverTest, ElpaDouble)
{
    for (const auto& s : sizes)
    {
        const std::vector<double> a = hermitian_matrix<double>(s[0], s[0] + 1.0);
        check_against_lapack<double>(&LR::elpa_linear_solver<double>, a, a, s[0], s[1], s[2]);
    }
}

TEST_F(ZeqLinearSolverTest, ElpaComplex)
{
    typedef std::complex<double> C;
    for (const auto& s : sizes)
    {
        const std::vector<C> a = hermitian_matrix<C>(s[0], s[0] + 1.0);
        check_against_lapack<C>(&LR::elpa_linear_solver<C>, a, a, s[0], s[1], s[2]);
    }
}

TEST_F(ZeqLinearSolverTest, ElpaSolvesHermitianPart)
{
    // a slightly non-symmetric input must be solved as its symmetric part (A + A^T)/2,
    // not as whatever its upper triangle happens to be
    for (const auto& s : sizes)
    {
        const int n = s[0];
        const std::vector<double> a_sym = hermitian_matrix<double>(n, n + 1.0);
        std::vector<double> a = a_sym;
        for (int j = 0; j < n; ++j)
        {
            for (int i = 0; i < j; ++i)
            {
                a[j * n + i] += 1e-2 * entry(i, j, 5);
                a[i * n + j] -= 1e-2 * entry(i, j, 5);
            }
        }
        check_against_lapack<double>(&LR::elpa_linear_solver<double>, a, a_sym, n, s[1], s[2]);
    }
}

TEST_F(ZeqLinearSolverTest, ElpaRejectsIndefinite)
{
    // single rank only: with several ranks ELPA's Cholesky deadlocks on an indefinite matrix
    // (only the rank owning the failing block returns), see `elpa_linear_solver`
    int nproc = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    if (nproc > 1) { return; }
    const int n = 20;
    const int nb = 3;
    const std::vector<double> a = hermitian_matrix<double>(n, -(n + 1.0));   // negative definite
    Parallel_2D pa;
    LR_Util::setup_2d_division(pa, nb, n, n);
    Parallel_2D pb;
    LR_Util::setup_2d_division(pb, nb, n, 1, pa.blacs_ctxt);
    std::vector<double> a_loc = to_local(a, pa);
    std::vector<double> b_loc = to_local(rhs_matrix<double>(n, 1), pb);
    EXPECT_THROW(LR::elpa_linear_solver<double>(a_loc.data(), b_loc.data(), pa, pb), std::runtime_error);
}
#endif


TEST(ZeqCGReview, SolvesIdentityWithoutDumpingProductionVector)
{
    double rhs[] = {2.0, -3.0};
    double z[] = {0.0, 0.0};
    const std::function<void(const double*, double*)> identity =
        [](const double* x, double* y) { y[0] = x[0]; y[1] = x[1]; };
    testing::internal::CaptureStdout();
    LR::solve_Z_CG(z, rhs, 2, 1, identity, false);
    const std::string output = testing::internal::GetCapturedStdout();
    EXPECT_NEAR(z[0], rhs[0], 1e-12);
    EXPECT_NEAR(z[1], rhs[1], 1e-12);
    EXPECT_EQ(output.find("Final Z-vector:"), std::string::npos);
}

TEST(ZeqCGReview, ZeroRhsIsConverged)
{
    double rhs[] = {0.0, 0.0};
    double z[] = {9.0, -1.0};
    const std::function<void(const double*, double*)> identity =
        [](const double* x, double* y) { y[0] = x[0]; y[1] = x[1]; };
    LR::solve_Z_CG(z, rhs, 2, 1, identity, false);
    EXPECT_DOUBLE_EQ(z[0], 0.0);
    EXPECT_DOUBLE_EQ(z[1], 0.0);
}

TEST(ZeqCGReview, RejectsAnUnsolvableSystem)
{
    double rhs[] = {1.0, 2.0};
    double z[] = {0.0, 0.0};
    const std::function<void(const double*, double*)> zero =
        [](const double*, double* y) { y[0] = 0.0; y[1] = 0.0; };
    EXPECT_THROW(LR::solve_Z_CG(z, rhs, 2, 1, zero, false), std::runtime_error);
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
    MPI_Finalize();
    return result;
}
