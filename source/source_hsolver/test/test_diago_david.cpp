#include "source_hsolver/diago_david.h"
#include "source_hsolver/diago_dav_subspace.h"
#include "source_hsolver/diago_cg.h"
#include "source_base/parallel_cell.h"
#include "source_base/parallel_comm.h"
#include "source_base/parallel_global.h"
#include "source_base/parallel_reduce.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <gtest/gtest.h>
#include <vector>

namespace
{
// A distributed diagonal H with S=I. ld must remain the physical local size,
// including zero: artificial vector padding must never reach the operator.
template <typename T>
class DiagonalOperator final : public hsolver::HSOperator<T>
{
  public:
    DiagonalOperator(const int offset, const int size) : offset_(offset), size_(size) {}
    void update_k(const int ik) override {}
    void hpsi(const T* x, T* hx, const int ld, const int nvec) const override
    {
        EXPECT_EQ(ld, size_);

        for (int band = 0; band < nvec; ++band)
        {
            for (int row = 0; row < size_; ++row)
            {
                hx[band * ld + row] = (offset_ + row + 1.0) * x[band * ld + row];
            }
        }
    }
    void spsi(const T* x, T* sx, const int ld, const int nvec) const override
    {
        EXPECT_EQ(ld, size_);
        const int count = size_ * nvec;
        std::copy_n(x, count, sx);
    }

  private:
    const int offset_;
    const int size_;
};

template <typename T>
T initial_value(const int row, const int band)
{
    return std::sin((row + 1) * (band + 1) * 0.37) + 0.2 * std::cos((row + 2) * 0.51);
}

template <>
std::complex<double> initial_value(const int row, const int band)
{
    const double real = initial_value<double>(row, band);
    const double imag = 0.3 * std::cos((row + 1) * (band + 2) * 0.29);
    return {real, imag};
}

template <typename T>
class DavidsonOwnershipTest : public ::testing::Test {};
using ScalarTypes = ::testing::Types<double, std::complex<double>>;
TYPED_TEST_SUITE(DavidsonOwnershipTest, ScalarTypes);

TYPED_TEST(DavidsonOwnershipTest, PhysicalEmptyAndUnevenBlocks)
{
    using T = TypeParam;
    const auto domain = ModuleBase::world_comm_domain();
    const int rank = domain.rank();
    const int processes = domain.size();
    const int dimension = 12;
    const int bands = 2;
#ifdef __MPI
    const hsolver::diag_comm_info comm = {POOL_WORLD, rank, processes};
#else
    const hsolver::diag_comm_info comm = {rank, processes};
#endif
    for (int layout = 0; layout < 3; ++layout)
    {
        std::vector<int> sizes(processes, 0);
        for (int process = 0; process < processes; ++process)
        {
            sizes[process] = dimension / processes;
            if (process < dimension % processes)
            {
                ++sizes[process];
            }
        }
        if (layout == 1 && processes > 1)
        {
            std::fill(sizes.begin(), sizes.end(), 0);
            sizes[1] = 1;
            sizes[processes - 1] += dimension - 1;
        }
        if (layout == 2)
        {
            std::fill(sizes.begin(), sizes.end(), 0);
            sizes[0] = dimension;
        }
        int offset = 0;
        for (int process = 0; process < rank; ++process)
        {
            offset += sizes[process];
        }
        const int local_size = sizes[rank];
        const int storage_size = std::max(local_size * bands, 1);
        std::vector<double> precondition(local_size, 1.0);
        for (int method = 0; method < 5; ++method)
        {
            SCOPED_TRACE(::testing::Message() << "layout=" << layout << " method=" << method << " rank=" << rank);
            std::vector<T> vectors(storage_size, T(999));
            for (int band = 0; band < bands; ++band)
            {
                for (int row = 0; row < local_size; ++row)
                {
                    vectors[band * local_size + row] = initial_value<T>(offset + row, band);
                }
            }
            double energies[bands] = {};
            const std::vector<double> thresholds(bands, 1e-13);
            DiagonalOperator<T> op(offset, local_size);
            if (method == 0)
            {
                hsolver::DiagoDavid<T> solver(precondition.data(), bands, local_size, 3, false, comm);
                solver.diag(op, local_size, vectors.data(), energies, thresholds, 300, 5, 0);
            }
            else if (method == 1)
            {
                hsolver::Diago_DavSubspace<T> solver(precondition, bands, local_size, 3, 1e-13, 300, comm, 0, 1);
                solver.diag(op, vectors.data(), local_size, energies, thresholds, false);
            }
            else
            {
                const bool initial_subspace = method == 3;
                hsolver::DiagoCG<T> solver("pw", "nscf", initial_subspace, comm, 1e-13, 300);
                const double* precondition_data = precondition.data();
                if (method == 4)
                {
                    precondition_data = nullptr;
                }
                solver.diag(op, local_size, bands, local_size, vectors.data(), energies, thresholds, precondition_data);
            }
            if (local_size == 0)
            {
                EXPECT_EQ(vectors[0], T(999));
            }
            for (int band = 0; band < bands; ++band)
            {
                EXPECT_NEAR(energies[band], band + 1.0, 1e-9);
                double residual = 0.0;
                double norm = 0.0;
                for (int row = 0; row < local_size; ++row)
                {
                    const double amplitude = std::abs(vectors[band * local_size + row]);
                    const double error = (offset + row + 1.0 - energies[band]) * amplitude;
                    residual += error * error;
                    norm += amplitude * amplitude;
                }
                Parallel_Reduce::reduce_pool(residual);
                Parallel_Reduce::reduce_pool(norm);
                EXPECT_LT(std::sqrt(residual), 2e-6);
                EXPECT_NEAR(norm, 1.0, 1e-10);
            }
        }
    }
}

TYPED_TEST(DavidsonOwnershipTest, ExactAndApproximateWarmStart)
{
    using T = TypeParam;
    const auto domain = ModuleBase::world_comm_domain();
    const int rank = domain.rank();
    const int processes = domain.size();
    const int bands = 2;
    const int dimension = std::max(4, 2 * processes);
#ifdef __MPI
    const hsolver::diag_comm_info comm = {POOL_WORLD, rank, processes};
#else
    const hsolver::diag_comm_info comm = {rank, processes};
#endif

    std::vector<int> sizes(processes, 0);
    for (int process = 0; process < processes; ++process)
    {
        sizes[process] = dimension / processes;
        if (process < dimension % processes)
        {
            ++sizes[process];
        }
    }
    int offset = 0;
    for (int process = 0; process < rank; ++process)
    {
        offset += sizes[process];
    }
    const int local_size = sizes[rank];
    const int storage_size = std::max(local_size * bands, 1);

    std::vector<double> precondition(std::max(local_size, 1), 1.0);
    for (int method = 0; method < 2; ++method)
    {
        SCOPED_TRACE(::testing::Message() << "method=" << method << " rank=" << rank);
        std::vector<T> vectors(storage_size, T(0));
        for (int band = 0; band < bands; ++band)
        {
            for (int row = 0; row < local_size; ++row)
            {
                const int global_row = offset + row;
                if (band == 0 && global_row == 0)
                {
                    vectors[band * local_size + row] = T(1);
                }
                if (band == 1 && (global_row == 1 || global_row == 2))
                {
                    vectors[band * local_size + row] = T(1.0 / std::sqrt(2.0));
                }
            }
        }

        double energies[bands] = {};
        const std::vector<double> thresholds(bands, 1e-8);
        DiagonalOperator<T> op(offset, local_size);
        if (method == 0)
        {
            hsolver::DiagoDavid<T> solver(precondition.data(), bands, local_size, 3, false, comm);
            solver.diag(op, local_size, vectors.data(), energies, thresholds, 100, 5, 0);
        }
        else
        {
            hsolver::Diago_DavSubspace<T> solver(precondition, bands, local_size, 3, 1e-8, 100, comm, 0, 1);
            solver.diag(op, vectors.data(), local_size, energies, thresholds, false);
        }

        for (int band = 0; band < bands; ++band)
        {
            EXPECT_NEAR(energies[band], band + 1.0, 1e-5);
            double residual = 0.0;
            double norm = 0.0;
            for (int row = 0; row < local_size; ++row)
            {
                const double amplitude = std::abs(vectors[band * local_size + row]);
                const double error = (offset + row + 1.0 - energies[band]) * amplitude;
                residual += error * error;
                norm += amplitude * amplitude;
            }
            Parallel_Reduce::reduce_pool(residual);
            Parallel_Reduce::reduce_pool(norm);
            EXPECT_LT(std::sqrt(residual), 1e-4);
            EXPECT_NEAR(norm, 1.0, 1e-6);
        }
    }
}

} // namespace

int main(int argc, char** argv)
{
    int processes = 1;
    int threads = 1;
    int rank = 0;
    Parallel_Global::read_pal_param(argc, argv, processes, threads, rank);
#ifdef __MPI
    POOL_WORLD = MPI_COMM_WORLD;
    KP_WORLD = MPI_COMM_NULL;
    INT_BGROUP = MPI_COMM_NULL;
    BP_WORLD = MPI_COMM_NULL;
    GRID_WORLD = MPI_COMM_NULL;
    DIAG_WORLD = MPI_COMM_NULL;
#endif
    ::testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
#ifdef __MPI
    Parallel_Global::finalize_mpi();
#endif
    return result;
}
