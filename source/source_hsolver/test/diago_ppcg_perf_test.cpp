#include "../diag_comm_info.h"
#include "../diago_bpcg.h"
#include "../diago_cg.h"
#include "../diago_david.h"
#include "../diago_iter_assist.h"
#include "../diago_ppcg.h"
#include "diago_mock.h"
#include "source_base/kernels/math_kernel_op.h"
#include "source_base/module_external/lapack_connector.h"
#include "source_psi/psi.h"

#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <mpi.h>
#include <random>
#include <vector>

#ifdef __MPI
#include "source_base/parallel_comm.h"
#endif

using T = std::complex<double>;

// LAPACK reference eigenvalues (values only)
void lapack_eigenvalues(int npw, const std::vector<T>& hm, double* e)
{
    std::vector<T> tmp = hm;
    int lwork = 2 * npw;
    std::vector<T> work(lwork);
    std::vector<double> rwork(3 * npw - 2);
    int info = 0;
    char jobz = 'N', uplo = 'U';
    zheev_(&jobz, &uplo, &npw, tmp.data(), &npw, e, work.data(), &lwork, rwork.data(), &info);
}

// Unified H|psi> via gemm
auto make_hpsi_func(const std::vector<T>& hmat, int dim)
{
    return [hmat, dim](T* psi_in, T* hpsi_out, const int ld_psi, const int nvec) {
        T one = T(1.0);
        T zero = T(0.0);
        base_device::DEVICE_CPU* ctx = {};
        ModuleBase::gemm_op<T, base_device::DEVICE_CPU>()('N',
                                                          'N',
                                                          dim,
                                                          nvec,
                                                          dim,
                                                          &one,
                                                          hmat.data(),
                                                          dim,
                                                          psi_in,
                                                          ld_psi,
                                                          &zero,
                                                          hpsi_out,
                                                          ld_psi);
    };
}

// S|psi> = |psi> (identity, norm-conserving)
auto spsi_identity = [](T* psi_in, T* spsi_out, const int ld_psi, const int nvec) {
    for (int i = 0; i < ld_psi * nvec; ++i)
    {
        spsi_out[i] = psi_in[i];
    }
};

struct PerfResult
{
    std::string name;
    double time = 0.0;
    double max_err = 0.0;
    bool converged = false;
};

// -------------------- PPCG --------------------
PerfResult test_ppcg(int nband,
                     int npw,
                     double ethr,
                     const psi::Psi<T>& psi0,
                     const std::function<void(T*, T*, const int, const int)>& hpsi_func,
                     double* precondition,
                     const std::vector<double>& e_ref)
{
    psi::Psi<T> psi(psi0);
    psi.fix_k(0);
    std::vector<double> en(nband, 0.0);

    hsolver::DiagoPPCG<T> ppcg(precondition);
    ppcg.init_iter(nband, nband, npw, npw);
    hsolver::DiagoIterAssist<T>::SCF_ITER = 1; // first SCF step
    std::vector<double> ethr_band(nband, ethr);

    double t1 = MPI_Wtime();
    ppcg.diag(hpsi_func, psi.get_pointer(), en.data(), ethr_band);
    double t2 = MPI_Wtime();

    double err = 0.0;
    for (int i = 0; i < nband; ++i)
    {
        err = std::max(err, std::abs(en[i] - e_ref[i]));
    }
    return {"PPCG", t2 - t1, err, err < 1e-2};
}

// -------------------- BPCG --------------------
PerfResult test_bpcg(int nband,
                     int npw,
                     double ethr,
                     const psi::Psi<T>& psi0,
                     const std::function<void(T*, T*, const int, const int)>& hpsi_func,
                     double* precondition,
                     const std::vector<double>& e_ref)
{
    psi::Psi<T> psi(psi0);
    psi.fix_k(0);
    std::vector<double> en(nband, 0.0);

    hsolver::DiagoBPCG<T> bpcg(precondition);
    bpcg.init_iter(nband, nband, npw, npw);
    hsolver::DiagoIterAssist<T>::SCF_ITER = 1;
    std::vector<double> ethr_band(nband, ethr);

    double t1 = MPI_Wtime();
    bpcg.diag(hpsi_func, psi.get_pointer(), en.data(), ethr_band);
    double t2 = MPI_Wtime();

    double err = 0.0;
    for (int i = 0; i < nband; ++i)
    {
        err = std::max(err, std::abs(en[i] - e_ref[i]));
    }
    return {"BPCG", t2 - t1, err, err < 1e-2};
}

// -------------------- CG --------------------
PerfResult test_cg(int nband,
                   int npw,
                   double ethr,
                   int maxiter,
                   const psi::Psi<T>& psi0,
                   const std::function<void(T*, T*, const int, const int)>& hpsi_func,
                   double* precondition,
                   const std::vector<double>& e_ref)
{
    psi::Psi<T> psi(psi0);
    psi.fix_k(0);
    std::vector<double> en(nband, 0.0);

    hsolver::DiagoCG<T> cg("pw", "scf");
    hsolver::DiagoIterAssist<T>::PW_DIAG_NMAX = maxiter;
    hsolver::DiagoIterAssist<T>::PW_DIAG_THR = ethr;
    std::vector<double> ethr_band(nband, ethr);

    double t1 = MPI_Wtime();
    cg.diag(hpsi_func, spsi_identity, npw, nband, npw, psi.get_pointer(), en.data(), ethr_band, precondition);
    double t2 = MPI_Wtime();

    double err = 0.0;
    for (int i = 0; i < nband; ++i)
    {
        err = std::max(err, std::abs(en[i] - e_ref[i]));
    }
    return {"CG", t2 - t1, err, err < 1e-2};
}

// -------------------- Davidson --------------------
PerfResult test_david(int nband,
                      int npw,
                      double ethr,
                      int maxiter,
                      const psi::Psi<T>& psi0,
                      const std::function<void(T*, T*, const int, const int)>& hpsi_func,
                      double* precondition,
                      const std::vector<double>& e_ref)
{
    psi::Psi<T> psi(psi0);
    psi.fix_k(0);
    std::vector<double> en(nband, 0.0);

#ifdef __MPI
    int rank = 0, nproc = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    const hsolver::diag_comm_info comm_info(MPI_COMM_WORLD, rank, nproc);
#else
    const hsolver::diag_comm_info comm_info(0, 1);
#endif

    hsolver::DiagoDavid<T> david(precondition, nband, npw, 4, false, comm_info);
    std::vector<double> ethr_band(nband, ethr);

    double t1 = MPI_Wtime();
    david.diag(hpsi_func, spsi_identity, npw, psi.get_pointer(), en.data(), ethr_band, maxiter);
    double t2 = MPI_Wtime();

    double err = 0.0;
    for (int i = 0; i < nband; ++i)
    {
        err = std::max(err, std::abs(en[i] - e_ref[i]));
    }
    return {"Davidson", t2 - t1, err, err < 1e-2};
}

// ============================================================
int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0, nproc = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

#ifdef __MPI
    BP_WORLD = MPI_COMM_WORLD;
#endif

    // ---------- test parameters ----------
    int nband = 20;
    int npw = 500;
    int sparsity = 0; // 0 = dense
    double ethr = 1e-5;
    int maxiter = 300;
    // -------------------------------------

    // generate Hamiltonian, precondition and initial guess
    HPsi<T> hpsi_gen(nband, npw, sparsity);
    DIAGOTEST::hmatrix = hpsi_gen.hamilt();
    DIAGOTEST::npw = npw;
    DIAGOTEST::npw_local = new int[1];
    DIAGOTEST::npw_local[0] = npw;
    DIAGOTEST::hmatrix_local = DIAGOTEST::hmatrix;

    double* precondition = hpsi_gen.precond();

    // LAPACK reference
    std::vector<double> e_ref(npw);
    lapack_eigenvalues(npw, DIAGOTEST::hmatrix, e_ref.data());

    // initial psi guess (perturbed eigenvectors)
    psi::Psi<T> psi0(1, nband, npw, npw, true);
    std::default_random_engine p(1);
    std::uniform_int_distribution<unsigned> u(1, 10);
    for (int i = 0; i < nband; ++i)
    {
        for (int j = 0; j < npw; ++j)
        {
            double r = static_cast<double>(u(p)) / 10.0;
            psi0(0, i, j) = DIAGOTEST::hmatrix[j * npw + i] * r;
        }
    }

    auto hpsi_func = make_hpsi_func(DIAGOTEST::hmatrix_local, npw);

    // run benchmarks
    PerfResult r_ppcg = test_ppcg(nband, npw, ethr, psi0, hpsi_func, precondition, e_ref);
    PerfResult r_bpcg = test_bpcg(nband, npw, ethr, psi0, hpsi_func, precondition, e_ref);
    PerfResult r_cg = test_cg(nband, npw, ethr, maxiter, psi0, hpsi_func, precondition, e_ref);
    PerfResult r_david = test_david(nband, npw, ethr, maxiter, psi0, hpsi_func, precondition, e_ref);

    if (rank == 0)
    {
        std::cout << "\n========================================\n";
        std::cout << "  Diagonalization Performance Comparison\n";
        std::cout << "  nband=" << nband << ", npw=" << npw << ", sparsity=" << sparsity << "\n";
        std::cout << "========================================\n";
        std::cout << std::setw(10) << "Method" << std::setw(14) << "Time(s)" << std::setw(14) << "MaxError"
                  << std::setw(8) << "OK" << "\n";
        std::cout << "----------------------------------------\n";
        auto print = [](const PerfResult& r) {
            std::cout << std::setw(10) << r.name << std::setw(14) << std::scientific << std::setprecision(3) << r.time
                      << std::setw(14) << r.max_err << std::setw(8) << (r.converged ? "Yes" : "No") << "\n";
        };
        print(r_ppcg);
        print(r_bpcg);
        print(r_cg);
        print(r_david);
        std::cout << "========================================\n\n";
    }

    delete[] DIAGOTEST::npw_local;
    MPI_Finalize();
    return 0;
}
