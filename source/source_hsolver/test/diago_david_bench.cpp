/**
 * Davidson benchmark: measures runtime and iterations for configurable test cases.
 * Outputs CSV lines: npw,nband,sparsity,mpi_procs,omp_threads,iterations,time_ms,max_error
 */
#include "../diag_comm_info.h"
#include "../diago_david.h"
#include "../diago_iter_assist.h"
#include "diago_mock.h"
#include "source_base/kernels/math_kernel_op.h"
#include "source_basis/module_pw/test/test_tool.h"
#include "source_base/module_external/lapack_connector.h"
#include "source_hamilt/hamilt.h"
#include "source_pw/module_pwdft/hamilt_pw.h"
#include "source_psi/psi.h"
#include "source_base/parallel_comm.h"

#include <chrono>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <random>
#include <string>
#include <vector>

namespace
{

void lapackEigen(const int npw, std::vector<std::complex<double>>& hm, double* e)
{
    int lwork = 2 * npw;
    std::vector<std::complex<double>> work(lwork);
    std::vector<double> rwork(3 * npw - 2);
    int info = 0;
    char jobz = 'V';
    char uplo = 'U';
    zheev_(&jobz, &uplo, &npw, hm.data(), &npw, e, work.data(), &lwork, rwork.data(), &info);
    if (info != 0)
    {
        std::cerr << "zheev failed with info=" << info << std::endl;
    }
}

} // namespace

int main(int argc, char** argv)
{
    int nproc = 1, myrank = 0;

#ifdef __MPI
    int nproc_in_pool, kpar = 1, mypool, rank_in_pool;
    setupmpi(argc, argv, nproc, myrank);
    divide_pools(nproc, myrank, nproc_in_pool, kpar, mypool, rank_in_pool);
    MPI_Comm_split(MPI_COMM_WORLD, myrank, 0, &BP_WORLD);
    POOL_WORLD = MPI_COMM_WORLD;  // Required by DiagoDavid internal assertions
    GlobalV::NPROC_IN_POOL = nproc;
#else
    MPI_Init(&argc, &argv);
    POOL_WORLD = MPI_COMM_WORLD;
#endif

    int npw = (argc > 1) ? std::atoi(argv[1]) : 100;
    int nband = (argc > 2) ? std::atoi(argv[2]) : 10;
    int sparsity = (argc > 3) ? std::atoi(argv[3]) : 6;
    double ethr = (argc > 4) ? std::atof(argv[4]) : 1e-7;

    int omp_threads = 1;
    const char* omp_env = std::getenv("OMP_NUM_THREADS");
    if (omp_env)
    {
        omp_threads = std::atoi(omp_env);
    }

    double max_error = 0.0;

    // Generate test problem
    HPsi<std::complex<double>> hpsi_mock(nband, npw, sparsity);
    DIAGOTEST::hmatrix = hpsi_mock.hamilt();
    DIAGOTEST::npw = npw;

    // Reference eigenvalues
    std::vector<double> e_lapack(npw, 0.0);
    auto h_lapack = DIAGOTEST::hmatrix;
    lapackEigen(npw, h_lapack, e_lapack.data());
#ifdef __MPI
    MPI_Bcast(e_lapack.data(), npw, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

    // Initial psi with perturbation
    psi::Psi<std::complex<double>> psi;
    psi.resize(1, nband, npw);
    std::default_random_engine engine(7);
    std::uniform_real_distribution<double> dist(0.2, 1.0);
    for (int ib = 0; ib < nband; ++ib)
    {
        for (int ig = 0; ig < npw; ++ig)
        {
            psi(ib, ig) = h_lapack[ig + ib * npw] * dist(engine);
        }
    }

    // MPI distribution: each process keeps full data for correct benchmark
    psi::Psi<std::complex<double>> psi_local;
    DIAGOTEST::npw_local = new int[nproc];
    double* precondition_local = nullptr;
#ifdef __MPI
    DIAGOTEST::cal_division(DIAGOTEST::npw);
    DIAGOTEST::hmatrix_local = DIAGOTEST::hmatrix;
    for (int i = 0; i < nproc; i++) {
        DIAGOTEST::npw_local[i] = DIAGOTEST::npw;
    }
    psi_local = psi;
    precondition_local = new double[DIAGOTEST::npw];
    for (int ig = 0; ig < DIAGOTEST::npw; ++ig)
    {
        precondition_local[ig] = hpsi_mock.precond()[ig];
    }
#else
    DIAGOTEST::hmatrix_local = DIAGOTEST::hmatrix;
    DIAGOTEST::npw_local[0] = DIAGOTEST::npw;
    psi_local = psi;
    precondition_local = new double[DIAGOTEST::npw];
    for (int ig = 0; ig < DIAGOTEST::npw; ++ig)
    {
        precondition_local[ig] = hpsi_mock.precond()[ig];
    }
#endif

    psi_local.fix_k(0);
    using T = std::complex<double>;
    const int dim = DIAGOTEST::npw;
    const std::vector<T>& h_mat = DIAGOTEST::hmatrix_local;
    auto hpsi_func = [h_mat, dim](T* psi_in, T* hpsi_out, const int ld_psi, const int nvec) {
        const T one(1.0);
        const T zero(0.0);
        ModuleBase::gemm_op<T, base_device::DEVICE_CPU>()(
            'N', 'N',
            dim, nvec, dim,
            &one,
            h_mat.data(), dim,
            psi_in, ld_psi,
            &zero,
            hpsi_out, ld_psi);
    };

    // S = I (identity), so spsi is just a copy of psi_in
    auto spsi_func = [](T* psi_in, T* spsi_out, const int ld_psi, const int nvec) {
        std::copy(psi_in, psi_in + static_cast<size_t>(ld_psi) * nvec, spsi_out);
    };

    const int ld_psi = psi_local.get_current_ngk();
    const int david_ndim = 4;
    const int david_maxiter = 200;

#ifdef __MPI
    hsolver::diag_comm_info diag_comm(MPI_COMM_WORLD, myrank, nproc);
#else
    hsolver::diag_comm_info diag_comm(myrank, nproc);
#endif

    hsolver::DiagoDavid<T> david(precondition_local, nband, npw, david_ndim, diag_comm);

    std::vector<double> eigen(nband, 0.0);
    std::vector<double> ethr_band(nband, ethr);

    auto t_start = std::chrono::high_resolution_clock::now();
    int niter = david.diag(hpsi_func, spsi_func, ld_psi, psi_local.get_pointer(),
                           eigen.data(), ethr_band, david_maxiter);
    auto t_end = std::chrono::high_resolution_clock::now();
    double elapsed_ms = std::chrono::duration<double, std::milli>(t_end - t_start).count();

    for (int ib = 0; ib < nband; ++ib)
    {
        double err = std::abs(eigen[ib] - e_lapack[ib]);
        if (err > max_error)
        {
            max_error = err;
        }
    }

    if (myrank == 0)
    {
        std::cout << npw << "," << nband << "," << sparsity << ","
                  << nproc << "," << omp_threads << "," << niter << ","
                  << elapsed_ms << "," << max_error << std::endl;
    }

    delete[] DIAGOTEST::npw_local;
    delete[] precondition_local;

    MPI_Finalize();
    return 0;
}
