#include "source_base/parallel_global.h"
#include "source_hsolver/hsolver_pw_tddft.h"
#include "source_hsolver/test/linear_test_utils.h"

#include <algorithm>
#include <cmath>
#include <gtest/gtest.h>
#include <sstream>

namespace
{
using T = linear_test::Complex;
using Device = base_device::DEVICE_CPU;

class DenseHamiltonian final : public hsolver::HSOperator<T>
{
  public:
    std::vector<std::vector<T>> matrices;
    std::vector<int> dimensions;

    explicit DenseHamiltonian(const std::vector<int>& sizes) : dimensions(sizes)
    {
        for (const int n: sizes)
        {
            matrices.push_back(linear_test::hermitian<T>(n));
        }
    }
    void update_k(const int ik) override
    {
        current_k_ = ik;
    }
    void hpsi(const T* x, T* y, const int ld, const int nvec) const override
    {
        const int n = dimensions.at(current_k_);
        const std::vector<T>& h = matrices.at(current_k_);
        for (int band = 0; band < nvec; ++band)
        {
            for (int i = 0; i < n; ++i)
            {
                y[band * ld + i] = T(0);
                for (int j = 0; j < n; ++j)
                {
                    y[band * ld + i] += h[j * n + i] * x[band * ld + j];
                }
            }
        }
    }
    void spsi(const T* x, T* y, const int ld, const int nvec) const override
    {
        FAIL() << "PW CN propagation must not request an overlap operator";
    }

  private:
    int current_k_ = -1;
};

void initialize_basis(const std::vector<int>& sizes, ModulePW::PW_Basis_K* basis)
{
    basis->nks = sizes.size();
    basis->npwk_max = *std::max_element(sizes.begin(), sizes.end());
    basis->tpiba = 0.8;
    basis->tpiba2 = 0.64;
    // These arrays are owned and released by PW_Basis/PW_Basis_K.
    basis->npwk = new int[basis->nks];
    basis->kvec_c = new ModuleBase::Vector3<double>[basis->nks];
    const int size = basis->nks * basis->npwk_max;
    basis->gcar = new ModuleBase::Vector3<double>[size];
    basis->gk2 = new double[size]();
    for (int ik = 0; ik < basis->nks; ++ik)
    {
        basis->npwk[ik] = sizes[ik];
        basis->kvec_c[ik] = ModuleBase::Vector3<double>(0.1 * ik, 0.07, -0.03);
        for (int i = 0; i < sizes[ik]; ++i)
        {
            const int index = ik * basis->npwk_max + i;
            basis->gcar[index] = ModuleBase::Vector3<double>(0.2 * i, 0.1 * (i % 3), -0.15);
            const ModuleBase::Vector3<double> q = basis->gcar[index] + basis->kvec_c[ik];
            basis->gk2[index] = q.norm2();
        }
    }
}

void check_dense_step(const DenseHamiltonian& op, const psi::Psi<T>& previous, const psi::Psi<T>& current, const double dt)
{
    const int bands = current.get_nbands();
    for (int ik = 0; ik < current.get_nk(); ++ik)
    {
        previous.fix_k(ik);
        current.fix_k(ik);
        const int n = op.dimensions[ik];
        std::vector<T> lhs(n * n);
        std::vector<T> rhs_matrix(n * n);
        for (int j = 0; j < n; ++j)
        {
            for (int i = 0; i < n; ++i)
            {
                const T identity = i == j ? T(1) : T(0);
                lhs[j * n + i] = identity + T(0, dt / 4.0) * op.matrices[ik][j * n + i];
                rhs_matrix[j * n + i] = identity - T(0, dt / 4.0) * op.matrices[ik][j * n + i];
            }
        }
        std::vector<T> packed(n * bands);
        for (int band = 0; band < bands; ++band)
        {
            for (int i = 0; i < n; ++i)
            {
                packed[band * n + i] = previous(band, i);
            }
        }
        const std::vector<T> rhs = linear_test::multiply(rhs_matrix, packed, n, bands);
        const std::vector<T> reference = linear_test::lapack_solve(lhs, rhs, n, bands);
        for (int band = 0; band < bands; ++band)
        {
            for (int i = 0; i < n; ++i)
            {
                EXPECT_LT(std::abs(current(band, i) - reference[band * n + i]), 1e-10);
            }
            for (int other = 0; other < bands; ++other)
            {
                T overlap(0);
                for (int i = 0; i < n; ++i)
                {
                    overlap += std::conj(current(other, i)) * current(band, i);
                }
                const T expected = band == other ? T(1) : T(0);
                EXPECT_LT(std::abs(overlap - expected), 1e-10);
            }
        }
    }
}

TEST(PWTDDFT, DenseCNCorrectorsAndConservation)
{
    const std::vector<int> sizes{9, 6};
    const int bands = 3;
    const int ld = 12;
    const double dt = 0.27;
    const ModuleBase::Vector3<double> shift(0.13, -0.07, 0.02);
#ifdef __MPI
    const hsolver::diag_comm_info comm(MPI_COMM_SELF, 0, 1);
#else
    const hsolver::diag_comm_info comm(0, 1);
#endif
    for (const std::string method: {"bicgstab", "cgs"})
    {
        for (const std::string precond: {"none", "kinetic"})
        {
            SCOPED_TRACE(method + "/" + precond);
            ModulePW::PW_Basis_K basis;
            initialize_basis(sizes, &basis);
            DenseHamiltonian op(sizes);
            std::ostringstream log;
            hsolver::HSolverPWTDDFT<T, Device> solver(basis, method, precond, 1e-13, 100, true, comm, log);
            psi::Psi<T> previous(2, bands, ld, sizes, true);
            for (int ik = 0; ik < 2; ++ik)
            {
                previous.fix_k(ik);
                for (int band = 0; band < bands; ++band)
                {
                    for (int i = 0; i < ld; ++i)
                    {
                        const double phase = 2.0 * std::acos(-1.0) * i * band / sizes[ik];
                        previous(band, i) = i < sizes[ik] ? std::polar(1.0 / std::sqrt(sizes[ik]), phase) : T(0);
                    }
                }
            }
            psi::Psi<T> current(previous);
            for (int step = 1; step <= 12; ++step)
            {
                solver.solve(op, previous, &current, dt, shift, step, 1, false, log);
                check_dense_step(op, previous, current, dt);
                if (step == 1)
                {
                    // A corrector rebuilds the RHS from the fixed previous state, even with a changed H and guess.
                    for (int ik = 0; ik < 2; ++ik)
                    {
                        const int n = sizes[ik];
                        current.fix_k(ik);
                        for (int i = 0; i < n; ++i)
                        {
                            op.matrices[ik][i * n + i] += 0.15 * (i + 1);
                        }
                        for (int band = 0; band < bands; ++band)
                        {
                            for (int i = 0; i < n; ++i)
                            {
                                current(band, i) *= T(0.8, 0.1);
                            }
                        }
                    }
                    solver.solve(op, previous, &current, dt, shift, step, 2, false, log);
                    check_dense_step(op, previous, current, dt);
                }
                previous = current;
            }
            // Evaluate expectations with a different endpoint H, not the last propagation H.
            for (int ik = 0; ik < 2; ++ik)
            {
                for (int i = 0; i < sizes[ik]; ++i)
                {
                    op.matrices[ik][i * sizes[ik] + i] += 0.3;
                }
            }
            ModuleBase::matrix energies(2, bands);
            solver.cal_band_energy(op, current, &energies);
            for (int ik = 0; ik < 2; ++ik)
            {
                current.fix_k(ik);
                const int n = sizes[ik];
                for (int band = 0; band < bands; ++band)
                {
                    T expected(0);
                    for (int j = 0; j < n; ++j)
                    {
                        for (int i = 0; i < n; ++i)
                        {
                            expected += std::conj(current(band, i)) * op.matrices[ik][j * n + i] * current(band, j);
                        }
                    }
                    EXPECT_NEAR(energies(ik, band), expected.real(), 1e-12);
                }
            }
        }
    }
}
} // namespace

int main(int argc, char** argv)
{
    int nproc = 1;
    int nthread = 1;
    int rank = 0;
    Parallel_Global::read_pal_param(argc, argv, nproc, nthread, rank);
#ifdef __MPI
    // This numerical test creates no application pools; the shared finalizer owns these handles.
    POOL_WORLD = MPI_COMM_NULL;
    KP_WORLD = MPI_COMM_NULL;
    INT_BGROUP = MPI_COMM_NULL;
    BP_WORLD = MPI_COMM_NULL;
    GRID_WORLD = MPI_COMM_NULL;
    DIAG_WORLD = MPI_COMM_NULL;
#endif
    testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
#ifdef __MPI
    Parallel_Global::finalize_mpi();
#endif
    return result;
}
