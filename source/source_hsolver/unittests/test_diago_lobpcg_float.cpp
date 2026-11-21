#include "source_hsolver/diago_lobpcg.h"
#include "source_hsolver/rr_solver.h"
#include "source_base/module_container/base/third_party/lapack.h"
#include "source_base/parallel_cell.h"
#include "source_base/parallel_comm.h"
#include "source_base/parallel_device.h"
#include "source_basis/module_pw/test/test_tool.h"

#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include <complex>
#include <vector>

using Locking = hsolver::LobpcgLocking;
class LobpcgFloat : public ::testing::TestWithParam<Locking> {};
namespace
{
using Complex = std::complex<float>;

class DenseOperator : public hsolver::HSOperator<Complex>
{
  public:
    DenseOperator(const std::vector<Complex>& matrix, const int n, const hsolver::diag_comm_info& comm)
        : n_(n), comm_(comm), local_(n / comm.nproc), first_(comm.rank * local_), h_(local_ * n)
    {
        for (int i = 0; i < local_; ++i)
        {
            for (int j = 0; j < n_; ++j)
            {
                h_[i + j * local_] = matrix[first_ + i + j * n_];
            }
        }
    }
    void update_k(const int) override {}
    void hpsi(const Complex* x, Complex* hx, const int ld, const int nvec) const override
    {
        for (int band = 0; band < nvec; ++band)
        {
            std::vector<Complex> full(n_, Complex(0));
            for (int j = 0; j < n_; ++j)
            {
                for (int i = 0; i < local_; ++i)
                {
                    full[j] += h_[i + j * local_] * x[i + band * ld];
                }
            }
#ifdef __MPI
            if (comm_.nproc > 1)
            {
                Parallel_Common::reduce_data(full.data(), n_, comm_.comm);
            }
#endif
            for (int i = 0; i < local_; ++i)
            {
                hx[i + band * ld] = full[first_ + i];
            }
        }
    }
    void spsi(const Complex* x, Complex* sx, const int ld, const int nvec) const override
    {
        std::copy(x, x + static_cast<std::size_t>(ld) * nvec, sx);
    }

  private:
    int n_;
    hsolver::diag_comm_info comm_;
    int local_;
    int first_;
    std::vector<Complex> h_;
};

hsolver::diag_comm_info comm_info()
{
#ifdef __MPI
    const auto world = ModuleBase::world_comm_domain();
    return hsolver::diag_comm_info(POOL_WORLD, world.rank(), world.size());
#else
    return hsolver::diag_comm_info(0, 1);
#endif
}

TEST_P(LobpcgFloat, ConvergesAndMatchesLapack)
{
    const auto comm = comm_info();
    const int rank = comm.rank;
    const int size = comm.nproc;
    const int n = 40;
    const int nband = 4;
    ASSERT_EQ(n % size, 0);
    const int local = n / size;
    std::vector<Complex> h(n * n, Complex(0));
    for (int j = 0; j < n; ++j)
    {
        for (int i = 0; i < n; ++i)
        {
            if (i == j)
            {
                h[i + j * n] = Complex(j + 1.0f, 0);
            }
            else
            {
                const Complex value(0.02f / (std::abs(i - j) + 1.0f), 0.003f * (i < j ? 1 : -1));
                h[i + j * n] = value;
                h[j + i * n] = std::conj(value);
            }
        }
    }
    DenseOperator op(h, n, comm);
    std::vector<float> reference(n);
    if (rank == 0)
    {
        char jobz = 'N';
        char uplo = 'U';
        int lwork = 2 * n;
        std::vector<Complex> work(lwork);
        std::vector<float> rwork(3 * n - 2);
        int info = 0;
        ct::lapackConnector::heev(jobz, uplo, n, h.data(), n, reference.data(), work.data(), lwork, rwork.data(), info);
        ASSERT_EQ(info, 0);
    }
#ifdef __MPI
    if (size > 1)
    {
        Parallel_Common::bcast_data(reference.data(), n, comm.comm, 0);
    }
#endif
    std::vector<Complex> psi(local * nband);
    for (int j = 0; j < nband; ++j)
    {
        for (int i = 0; i < local; ++i)
        {
            psi[i + j * local] = Complex(std::sin(0.17f * (i + 1) * (j + 2)) + 0.3f, 0.1f);
        }
    }
    std::vector<float> prec(local, 1.0f);
    std::vector<float> eigenvalues(nband);
    const hsolver::LobpcgBasisLayout layout = {local, 1};
    const Locking locking = LobpcgFloat::GetParam();
    hsolver::rr::ParallelRREigensolver solver_rr;
    hsolver::rr::GpuRREigensolver solver_rr_gpu;
    hsolver::DiagoLOBPCG<Complex> solver(prec.data(), nband, local, nband + 10, layout, comm, 0, 0, false, locking, solver_rr, solver_rr_gpu);
    ASSERT_TRUE(solver.diag(op, psi.data(), local, eigenvalues.data(), 1e-4, 150, {}));
    for (int band = 0; band < nband; ++band)
    {
        EXPECT_NEAR(eigenvalues[band], reference[band], 3e-3);
        std::vector<Complex> hx(local);
        op.hpsi(psi.data() + band * local, hx.data(), local, 1);
        float residual = 0.0f;
        for (int i = 0; i < local; ++i)
        {
            residual += std::norm(hx[i] - eigenvalues[band] * psi[i + band * local]);
        }
#ifdef __MPI
        if (size > 1)
        {
            Parallel_Common::reduce_data(&residual, 1, comm.comm);
        }
#endif
        EXPECT_LT(std::sqrt(residual), 3e-2);
    }
}

class SpinorOperator : public hsolver::HSOperator<Complex>
{
  public:
    SpinorOperator(const int npw, const int stride, const int first)
        : npw_(npw), stride_(stride), first_(first)
    {
    }
    void update_k(const int) override {}
    void hpsi(const Complex* x, Complex* hx, const int ld, const int nvec) const override
    {
        std::fill(hx, hx + ld * nvec, Complex(0));
        for (int band = 0; band < nvec; ++band)
        {
            for (int spin = 0; spin < 2; ++spin)
            {
                for (int ig = 0; ig < npw_; ++ig)
                {
                    const int row = band * ld + spin * stride_ + ig;
                    const float value = spin * 40 + first_ + ig - 1.5f;
                    hx[row] = value * x[row];
                }
                for (int ig = npw_; ig < stride_; ++ig)
                {
                    EXPECT_EQ(x[band * ld + spin * stride_ + ig], Complex(0));
                }
            }
        }
    }
    void spsi(const Complex* x, Complex* sx, const int ld, const int nvec) const override
    {
        std::copy(x, x + ld * nvec, sx);
    }

  private:
    int npw_;
    int stride_;
    int first_;
};

TEST_P(LobpcgFloat, SpinorPadding)
{
    const auto comm = comm_info();
    ASSERT_EQ(40 % comm.nproc, 0);
    const int npw = 40 / comm.nproc;
    const int stride = npw + 5 + comm.rank;
    const int first = comm.rank * npw;
    const int dim = 2 * stride;
    const int nband = 6;
    const int nmax = nband + 10;
    SpinorOperator op(npw, stride, first);
    std::vector<Complex> psi(dim * nband, Complex(0.3f, 0.1f));
    std::vector<float> prec(dim, 1.0f);
    for (int spin = 0; spin < 2; ++spin)
    {
        for (int ig = 0; ig < npw; ++ig)
        {
            const int physical = spin * 40 + first + ig;
            const int row = spin * stride + ig;
            prec[row] = std::max(1.0f, physical - 1.5f);
            for (int j = 0; j < nband; ++j)
            {
                psi[row + j * dim] = Complex(std::sin(0.17f * (physical + 1) * (j + 2)), 0.1f);
            }
        }
    }
    const hsolver::LobpcgBasisLayout layout = {npw, 2};
    const Locking locking = LobpcgFloat::GetParam();
    hsolver::rr::ParallelRREigensolver solver_rr;
    hsolver::rr::GpuRREigensolver solver_rr_gpu;
    hsolver::DiagoLOBPCG<Complex> solver(prec.data(), nband, dim, nmax, layout, comm, 0, 0, false, locking, solver_rr, solver_rr_gpu);
    std::vector<float> eigenvalues(nband);
    ASSERT_TRUE(solver.diag(op, psi.data(), dim, eigenvalues.data(), 1e-5, 200, {}));
    std::vector<Complex> hx(dim * nband);
    op.hpsi(psi.data(), hx.data(), dim, nband);
    for (int j = 0; j < nband; ++j)
    {
        EXPECT_NEAR(eigenvalues[j], j - 1.5f, 1e-3);
        double residual = 0;
        for (int i = 0; i < dim; ++i)
        {
            residual += std::norm(hx[i + j * dim] - eigenvalues[j] * psi[i + j * dim]);
        }
#ifdef __MPI
        Parallel_Common::reduce_data(&residual, 1, comm.comm);
#endif
        EXPECT_LT(std::sqrt(residual), 1e-2);
    }
    hsolver::LinearAlgebra<Complex, base_device::DEVICE_CPU> algebra(comm);
    const std::vector<std::complex<double>> gram = algebra.gram(dim, dim, nband, psi.data());
    for (int j = 0; j < nband; ++j)
    {
        for (int i = 0; i < nband; ++i)
        {
            const std::complex<double> expected = i == j ? 1.0 : 0.0;
            EXPECT_LT(std::abs(gram[i + j * nband] - expected), 2e-5);
        }
    }
}

TEST_P(LobpcgFloat, ExactInitialRitzPairs)
{
    const auto comm = comm_info();
    const int n = 40;
    const int local = n / comm.nproc;
    const int nband = 4;
    std::vector<Complex> matrix(n * n, Complex(0));
    for (int i = 0; i < n; ++i)
    {
        matrix[i + i * n] = i + 1.0f;
    }
    DenseOperator op(matrix, n, comm);
    std::vector<Complex> psi(local * nband, Complex(0));
    for (int j = 0; j < nband; ++j)
    {
        const int row = j - comm.rank * local;
        if (row >= 0 && row < local)
        {
            psi[row + j * local] = 1.0f;
        }
    }
    const std::vector<Complex> initial(psi);
    std::vector<float> prec(local, 1.0f);
    std::vector<float> eigenvalues(nband);
    const hsolver::LobpcgBasisLayout layout = {local, 1};
    const Locking locking = LobpcgFloat::GetParam();
    hsolver::rr::ParallelRREigensolver solver_rr;
    hsolver::rr::GpuRREigensolver solver_rr_gpu;
    hsolver::DiagoLOBPCG<Complex> solver(prec.data(), nband, local, nband, layout, comm, 0, 0, false, locking, solver_rr, solver_rr_gpu);
    ASSERT_TRUE(solver.diag(op, psi.data(), local, eigenvalues.data(), 1e-5, 100, {}));
    EXPECT_EQ(solver.iterations(), 0);
    for (int j = 0; j < nband; ++j)
    {
        EXPECT_NEAR(eigenvalues[j], j + 1.0f, 1e-6);
        for (int i = 0; i < local; ++i)
        {
            EXPECT_NEAR(std::abs(psi[i + j * local]), std::abs(initial[i + j * local]), 1e-6);
        }
    }
}
} // namespace

int main(int argc, char** argv)
{
    int myrank = 0;
#ifdef __MPI
    int nproc = 1;
    int nproc_in_pool = 0;
    int kpar = 1;
    int mypool = 0;
    int rank_in_pool = 0;
    setupmpi(argc, argv, nproc, myrank);
    divide_pools(nproc, myrank, nproc_in_pool, kpar, mypool, rank_in_pool);
#endif
    testing::InitGoogleTest(&argc, argv);
    testing::TestEventListeners& listeners = testing::UnitTest::GetInstance()->listeners();
    if (myrank != 0)
    {
        delete listeners.Release(listeners.default_result_printer());
    }
    const int result = RUN_ALL_TESTS();
#ifdef __MPI
    finishmpi();
#endif
    return result;
}

INSTANTIATE_TEST_SUITE_P(LockingModes, LobpcgFloat, ::testing::Values(Locking::Hard, Locking::Soft));
