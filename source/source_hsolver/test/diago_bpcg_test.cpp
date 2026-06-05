#include "source_base/inverse_matrix.h"
#include "source_base/module_external/lapack_connector.h"
#include "source_pw/module_pwdft/structure_factor.h"
#include "source_psi/psi.h"
#include "source_hamilt/hamilt.h"
#include "source_pw/module_pwdft/hamilt_pw.h"
#include "../diago_iter_assist.h"
#include "../diago_bpcg.h"
#include "../kernels/bpcg_kernel_op.h"
#include "diago_mock.h"
#include "mpi.h"
#include "source_base/global_variable.h"
#include "source_base/parallel_comm.h"
#include "source_basis/module_pw/test/test_tool.h"

#include <gtest/gtest.h>
#include <complex>
#include <random>
#ifdef _OPENMP
#include <omp.h>
#endif

/************************************************
 *  unit test of functions in Diago_BPCG
 ***********************************************/

/**
 * Class Diago_BPCG is an approach for eigenvalue problems
 * This unittest test the function Diago_BPCG::diag() for FPTYPE=double, Device=cpu
 * with different examples.
 *  - the Hermite matrices (npw=500,1000) produced using random numbers and with sparsity of 0%, 60%, 80%
 *  - the Hamiltonian matrix read from "data-H", produced by using out_hs in INPUT of a LCAO calculation
 *  - a 2x2 Hermite matrix for learning and checking
 *
 * Note:
 * The test is passed when the eignvalues are closed to these calculated by LAPACK.
 * It is used together with a header file diago_mock.h.
 * The default Hermite matrix generated here is real symmetric, one can add an imaginary part
 * by changing two commented out lines in diago_mock.h.
 *
 */

// call lapack in order to compare to bpcg
void lapackEigen(int &npw, std::vector<std::complex<double>> &hm, double *e, bool outtime = false)
{
    clock_t start, end;
    start = clock();
    int lwork = 2 * npw;
    std::complex<double> *work2 = new std::complex<double>[lwork];
    double *rwork = new double[3 * npw - 2];
    int info = 0;
    char tmp_c1 = 'V', tmp_c2 = 'U';
    zheev_(&tmp_c1, &tmp_c2, &npw, hm.data(), &npw, e, work2, &lwork, rwork, &info);
    end = clock();
    if (outtime) {
        std::cout << "Lapack Run time: " << (double)(end - start) / CLOCKS_PER_SEC << " S" << std::endl;
}
    delete[] rwork;
    delete[] work2;
}

class DiagoBPCGPrepare
{
  public:
    DiagoBPCGPrepare(int nband, int npw, int sparsity, bool reorder, double eps, int maxiter, double threshold)
        : nband(nband), npw(npw), sparsity(sparsity), reorder(reorder), eps(eps), maxiter(maxiter),
          threshold(threshold)
    {
#ifdef __MPI	
		MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
        MPI_Comm_rank(MPI_COMM_WORLD, &mypnum);
#endif	
    }

    int nband, npw, sparsity, maxiter, notconv;
    // eps is the convergence threshold within cg_diago
    double eps, avg_iter;
    bool reorder;
    double threshold;
    int nprocs=1, mypnum=0;
    // threshold is the comparison standard between bpcg and lapack

    void CompareEigen(double *precondition, bool check_vectors = false)
    {
        // calculate eigenvalues by LAPACK;
        double *e_lapack = new double[npw];
        auto ev = DIAGOTEST::hmatrix;
        if (mypnum == 0) {
            lapackEigen(npw, ev, e_lapack, false);
        }
    #ifdef __MPI
        MPI_Bcast(e_lapack, npw, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    #endif
        // initial guess of psi by perturbing lapack psi
        ModuleBase::ComplexMatrix psiguess(nband, npw);
        std::default_random_engine p(1);
        std::uniform_int_distribution<unsigned> u(1, 10);
        for (int i = 0; i < nband; i++)
        {
            for (int j = 0; j < npw; j++)
            {
		        double rand = static_cast<double>(u(p))/10.;
                // psiguess(i,j) = ev(j,i)*(1+rand);
                psiguess(i, j) = ev[j * npw + i] * rand;
            }
        }
        // run bpcg
	//======================================================================
        double *en = new double[npw];
        int ik = 1;
        psi::Psi<std::complex<double>> psi;
	    psi.resize(ik,nband,npw);
	    //psi.fix_k(0);
        for (int i = 0; i < nband; i++)
        {
            for (int j = 0; j < npw; j++)
            {
	            psi(i,j)=psiguess(i,j);
	        }
	    }	

        psi::Psi<std::complex<double>> psi_local;
        double* precondition_local;
        DIAGOTEST::npw_local = new int[nprocs];
#ifdef __MPI				
	    DIAGOTEST::cal_division(DIAGOTEST::npw);
        DIAGOTEST::divide_hpsi(psi, psi_local, DIAGOTEST::hmatrix, DIAGOTEST::hmatrix_local); //will distribute psi and Hmatrix to each process
	    precondition_local = new double[DIAGOTEST::npw_local[mypnum]];
	    DIAGOTEST::divide_psi<double>(precondition,precondition_local);	
#else
	    DIAGOTEST::hmatrix_local = DIAGOTEST::hmatrix;
	    DIAGOTEST::npw_local[0] = DIAGOTEST::npw;
	    psi_local = psi;
	    precondition_local = new double[DIAGOTEST::npw];
	    for(int i=0;i<DIAGOTEST::npw;i++) precondition_local[i] = precondition[i];
#endif
        hsolver::DiagoBPCG<std::complex<double>> bpcg(precondition_local);
        psi_local.fix_k(0);
        double start, end;
        start = MPI_Wtime();
        using T = std::complex<double>;
        const int dim = DIAGOTEST::npw;
        const std::vector<T> &h_mat = DIAGOTEST::hmatrix_local;
        auto hpsi_func = [h_mat, dim](T *psi_in, T *hpsi_out,
                                const int ld_psi, const int nvec) {
            auto one = std::make_unique<T>(1.0);
            auto zero = std::make_unique<T>(0.0);
            const T *one_ = one.get();
            const T *zero_ = zero.get();

            base_device::DEVICE_CPU *ctx = {};
            // hpsi_out(dim * nvec) = h_mat(dim * dim) * psi_in(dim * nvec)
            ModuleBase::gemm_op<T, base_device::DEVICE_CPU>()(
                'N', 'N',
                dim, nvec, dim,
                one_,
                h_mat.data(), dim,
                psi_in, ld_psi,
                zero_,
                hpsi_out, ld_psi);
        };
        const int ndim = psi_local.get_current_ngk();
        bpcg.init_iter(nband, nband, npw, ndim);
        std::vector<double> ethr_band(nband, 1e-5);
        // One diag() call has a relatively small internal iteration cap; do a few passes
        // to reach LAPACK-close eigenvalues for random dense problems.
        for (int pass = 0; pass < 4; ++pass) {
            bpcg.diag(hpsi_func, psi_local.get_pointer(), en, ethr_band);
        }
        end = MPI_Wtime();
        //if(mypnum == 0) printf("diago time:%7.3f\n",end-start);
        delete [] DIAGOTEST::npw_local;
	    delete [] precondition_local;
	    //======================================================================
        for (int i = 0; i < nband; i++)
        {
            EXPECT_NEAR(en[i], e_lapack[i], threshold);
        }

        if (check_vectors && nprocs == 1)
        {
            std::vector<std::complex<double>> hpsi_check(nband * npw);
            hpsi_func(psi_local.get_pointer(), hpsi_check.data(), npw, nband);

            for (int ib = 0; ib < nband; ++ib)
            {
                double norm = 0.0;
                double residual_norm = 0.0;
                for (int ig = 0; ig < npw; ++ig)
                {
                    const std::complex<double> psi_value = psi_local(ib, ig);
                    const std::complex<double> residual = hpsi_check[ib * npw + ig] - en[ib] * psi_value;
                    norm += std::norm(psi_value);
                    residual_norm += std::norm(residual);
                }
                EXPECT_NEAR(norm, 1.0, 1e-10);
                EXPECT_NEAR(std::sqrt(residual_norm), 0.0, 1e-8);
            }

            for (int ib = 0; ib < nband; ++ib)
            {
                for (int jb = ib + 1; jb < nband; ++jb)
                {
                    std::complex<double> overlap = 0.0;
                    for (int ig = 0; ig < npw; ++ig)
                    {
                        overlap += std::conj(psi_local(ib, ig)) * psi_local(jb, ig);
                    }
                    EXPECT_NEAR(std::abs(overlap), 0.0, 1e-10);
                }
            }
        }

        delete[] en;
        delete[] e_lapack;
    }
};

class DiagoBPCGTest : public ::testing::TestWithParam<DiagoBPCGPrepare>
{
};

TEST_P(DiagoBPCGTest, RandomHamilt)
{
    DiagoBPCGPrepare dcp = GetParam();
    //std::cout << "npw=" << dcp.npw << ", nband=" << dcp.nband << ", sparsity="
    //		  << dcp.sparsity << ", eps=" << dcp.eps << std::endl;
    hsolver::DiagoIterAssist<std::complex<double>>::PW_DIAG_NMAX = dcp.maxiter;
    hsolver::DiagoIterAssist<std::complex<double>>::PW_DIAG_THR = dcp.eps;
    hsolver::DiagoIterAssist<std::complex<double>>::SCF_ITER = 1;
    //std::cout<<"maxiter "<<hsolver::DiagoIterAssist<std::complex<double>>::PW_DIAG_NMAX<<std::endl;
    //std::cout<<"eps "<<hsolver::DiagoIterAssist<std::complex<double>>::PW_DIAG_THR<<std::endl;
    HPsi<std::complex<double>> hpsi(dcp.nband, dcp.npw, dcp.sparsity);
    DIAGOTEST::hmatrix = hpsi.hamilt();

    DIAGOTEST::npw = dcp.npw;
    // ModuleBase::ComplexMatrix psi = hpsi.psi();
    dcp.CompareEigen(hpsi.precond());
}

INSTANTIATE_TEST_SUITE_P(VerifyCG,
                         DiagoBPCGTest,
                         ::testing::Values(
                             // nband, npw, sparsity, reorder, eps, maxiter, threshold
                             DiagoBPCGPrepare(6, 120, 0, true, 1e-5, 200, 5e-2)
                            //  DiagoBPCGPrepare(20, 500, 6, true, 1e-5, 300, 5e-2)
                            //  DiagoBPCGPrepare(20, 1000, 8, true, 1e-5, 300, 5e-2),
                            //  DiagoBPCGPrepare(40, 1000, 8, true, 1e-6, 300, 5e-2)
                            )); 
                            //DiagoBPCGPrepare(40, 2000, 8, true, 1e-5, 500, 1e-2))); 
			    // the last one is passed but time-consumming.

// check that the mock class HPsi work well
// in generating a Hermite matrix
TEST(DiagoBPCGTest, Hamilt)
{
    int dim = 2;
    int nbnd = 2;
    HPsi<std::complex<double>> hpsi(nbnd, dim);
    std::vector<std::complex<double>> hm = hpsi.hamilt();
    EXPECT_EQ(DIAGOTEST::h_nr, 2);
    EXPECT_EQ(DIAGOTEST::h_nc, 2);
    EXPECT_EQ(hm[0].imag(), 0.0);
    EXPECT_EQ(hm[DIAGOTEST::h_nc + 1].imag(), 0.0);
    EXPECT_EQ(conj(hm[DIAGOTEST::h_nc]).real(), hm[1].real());
    EXPECT_EQ(conj(hm[DIAGOTEST::h_nc]).imag(), hm[1].imag());
}

// bpcg for a 2x2 matrix (analytic eigenvalues: (7±sqrt(5))/2)
TEST(DiagoBPCGTest, TwoByTwo)
{
    const int dim = 2;
    const int nband = 2;
    std::vector<std::complex<double>> hm(dim * dim);
    hm[0] = {4.0, 0.0};
    hm[1] = {1.0, 0.0};
    hm[2] = {1.0, 0.0};
    hm[3] = {3.0, 0.0};

    DiagoBPCGPrepare dcp(nband, dim, 0, true, 1e-8, 80, 1e-8);
    hsolver::DiagoIterAssist<std::complex<double>>::PW_DIAG_NMAX = dcp.maxiter;
    hsolver::DiagoIterAssist<std::complex<double>>::PW_DIAG_THR = dcp.eps;
    hsolver::DiagoIterAssist<std::complex<double>>::SCF_ITER = 1;

    // simple positive precondition
    double precond[dim] = {1.0, 1.0};
    DIAGOTEST::hmatrix = hm;
    DIAGOTEST::npw = dim;
    dcp.CompareEigen(precond, true);
}

TEST(BpcgKernelOpTest, ApplyEigenvaluesUsesLeadingDimension)
{
    using T = std::complex<double>;
    const int nbase = 4101;
    const int nbase_x = nbase + 3;
    const int notconv = 2;
    const T untouched = {-9.0, 4.0};

    std::vector<T> vectors(nbase_x * notconv);
    std::vector<T> result(nbase_x * notconv, untouched);
    const double eigenvalues[notconv] = {2.0, -0.5};

    for (int m = 0; m < notconv; ++m)
    {
        for (int i = 0; i < nbase_x; ++i)
        {
            vectors[m * nbase_x + i] = T(0.25 * (i + 1), -0.1 * (m + 1));
        }
    }

    hsolver::apply_eigenvalues_op<T, base_device::DEVICE_CPU>()(
        nbase, nbase_x, notconv, result.data(), vectors.data(), eigenvalues);

    for (int m = 0; m < notconv; ++m)
    {
        for (int i = 0; i < nbase; ++i)
        {
            EXPECT_EQ(result[m * nbase_x + i], eigenvalues[m] * vectors[m * nbase_x + i]);
        }
        for (int i = nbase; i < nbase_x; ++i)
        {
            EXPECT_EQ(result[m * nbase_x + i], untouched);
        }
    }
}

TEST(BpcgKernelOpTest, PreconditionUsesBandOffsetAndFormula)
{
    using T = std::complex<double>;
    const int dim = 4;
    const int nbase = 2;
    const int notconv = 2;
    std::vector<T> psi_iter((nbase + notconv) * dim);
    const std::vector<T> original = {
        {1.0, 0.0}, {2.0, 0.0}, {3.0, 0.0}, {4.0, 0.0},
        {5.0, 0.0}, {6.0, 0.0}, {7.0, 0.0}, {8.0, 0.0},
        {1.0, 2.0}, {2.0, 3.0}, {3.0, 4.0}, {4.0, 5.0},
        {2.0, -1.0}, {3.0, -2.0}, {4.0, -3.0}, {5.0, -4.0}};
    psi_iter = original;

    const double precondition[dim] = {1.0, 2.5, 4.0, 7.0};
    const double eigenvalues[notconv] = {0.5, 3.0};

    hsolver::precondition_op<T, base_device::DEVICE_CPU>()(
        dim, psi_iter.data(), nbase, notconv, precondition, eigenvalues);

    for (int i = 0; i < nbase * dim; ++i)
    {
        EXPECT_EQ(psi_iter[i], original[i]);
    }

    for (int m = 0; m < notconv; ++m)
    {
        for (int i = 0; i < dim; ++i)
        {
            const double x = std::abs(precondition[i] - eigenvalues[m]);
            const double denom = 0.5 * (1.0 + x + std::sqrt(1.0 + (x - 1.0) * (x - 1.0)));
            const int idx = (nbase + m) * dim + i;
            EXPECT_NEAR(psi_iter[idx].real(), (original[idx] / denom).real(), 1e-14);
            EXPECT_NEAR(psi_iter[idx].imag(), (original[idx] / denom).imag(), 1e-14);
        }
    }
}

TEST(BpcgKernelOpTest, RefreshProjectedMatricesOnlyTouchesDiagonal)
{
    using T = std::complex<double>;
    const int n = 3;
    const int ldh = 5;
    const T one = {1.0, 0.0};
    const T h_sentinel = {-1.0, 0.5};
    const T s_sentinel = {-2.0, 0.5};
    const T v_sentinel = {-3.0, 0.5};
    const double eigenvalues[n] = {0.25, 1.5, 3.75};

    std::vector<T> hcc(ldh * ldh, h_sentinel);
    std::vector<T> scc(ldh * ldh, s_sentinel);
    std::vector<T> vcc(ldh * ldh, v_sentinel);

    hsolver::refresh_hcc_scc_vcc_op<T, base_device::DEVICE_CPU>()(
        n, hcc.data(), scc.data(), vcc.data(), ldh, eigenvalues, one);

    for (int col = 0; col < ldh; ++col)
    {
        for (int row = 0; row < ldh; ++row)
        {
            const int idx = col * ldh + row;
            if (row == col && row < n)
            {
                EXPECT_EQ(hcc[idx], T(eigenvalues[row], 0.0));
                EXPECT_EQ(scc[idx], one);
                EXPECT_EQ(vcc[idx], one);
            }
            else
            {
                EXPECT_EQ(hcc[idx], h_sentinel);
                EXPECT_EQ(scc[idx], s_sentinel);
                EXPECT_EQ(vcc[idx], v_sentinel);
            }
        }
    }
}

TEST(BpcgKernelOpTest, ApplyEigenvaluesMatchesSingleThreadResult)
{
#ifndef _OPENMP
    GTEST_SKIP() << "OpenMP is not enabled in this build";
#else
    using T = std::complex<double>;
    const int nbase = 5000;
    const int nbase_x = nbase + 7;
    const int notconv = 3;
    std::vector<T> vectors(nbase_x * notconv);
    std::vector<T> result_single(nbase_x * notconv);
    std::vector<T> result_multi(nbase_x * notconv);
    const double eigenvalues[notconv] = {1.25, -2.0, 0.125};

    for (int m = 0; m < notconv; ++m)
    {
        for (int i = 0; i < nbase_x; ++i)
        {
            vectors[m * nbase_x + i] = T(0.01 * (i % 97) + m, -0.02 * (i % 31));
        }
    }

    const int old_threads = omp_get_max_threads();
    omp_set_num_threads(1);
    hsolver::apply_eigenvalues_op<T, base_device::DEVICE_CPU>()(
        nbase, nbase_x, notconv, result_single.data(), vectors.data(), eigenvalues);

    omp_set_num_threads(4);
    hsolver::apply_eigenvalues_op<T, base_device::DEVICE_CPU>()(
        nbase, nbase_x, notconv, result_multi.data(), vectors.data(), eigenvalues);
    omp_set_num_threads(old_threads);

    for (size_t i = 0; i < result_single.size(); ++i)
    {
        EXPECT_EQ(result_multi[i], result_single[i]);
    }
#endif
}

// check that lapack work well
// for an eigenvalue problem
/*TEST(DiagoBPCGTest, ZHEEV)
{
    int dim = 100;
    int nbnd = 2;
    HPsi hpsi(nbnd, dim);
    std::vector<std::complex<double>> hm = hpsi.hamilt();
    std::vector<std::complex<double>> hm_backup = hm;
    ModuleBase::ComplexMatrix eig(dim, dim);
    double e[dim];
    // using zheev to do a direct test
    lapackEigen(dim, hm, e);
    eig = transpose(hm, true) * hm_backup * hm;
    // for (int i=0;i<dim;i++) std::cout<< " e[i] "<<e[i]<<std::endl;
    for (int i = 0; i < dim; i++)
    {
        EXPECT_NEAR(e[i], eig(i, i).real(), 1e-10);
    }
}*/


TEST(DiagoBPCGTest, readH)
{
    // read Hamilt matrix from file data-H
    std::vector<std::complex<double>> hm;
    std::ifstream ifs;
    std::string filename = "H-KPoints-Si2.dat";
    ifs.open(filename);
    // open file and check status
    if (!ifs.is_open())
    {
        std::cout << "Error opening file " << filename << std::endl;
        exit(1);
    }
    DIAGOTEST::readh(ifs, hm);
    ifs.close();
    int dim = DIAGOTEST::npw;
    int nband = 10; // not nband < dim, here dim = 26 in data-H
    // nband, npw, sub, sparsity, reorder, eps, maxiter, threshold
    DiagoBPCGPrepare dcp(nband, dim, 0, true, 1e-5, 500, 1e-1);
    hsolver::DiagoIterAssist<std::complex<double>>::PW_DIAG_NMAX = dcp.maxiter;
    hsolver::DiagoIterAssist<std::complex<double>>::PW_DIAG_THR = dcp.eps;
    hsolver::DiagoIterAssist<std::complex<double>>::SCF_ITER = 1;
    HPsi<std::complex<double>> hpsi;
    hpsi.create(nband, dim);
    // use the matrix read from file
    DIAGOTEST::hmatrix = hm;
    DIAGOTEST::npw = dim;
    dcp.CompareEigen(hpsi.precond());
}

int main(int argc, char **argv)
{
	int nproc = 1, myrank = 0;

#ifdef __MPI
	int nproc_in_pool, kpar=1, mypool, rank_in_pool;
    setupmpi(argc,argv,nproc, myrank);
    divide_pools(nproc, myrank, nproc_in_pool, kpar, mypool, rank_in_pool);
    // In unit tests we don't do band-parallel splitting; keep BP_WORLD as the full pool communicator.
    MPI_Comm_dup(POOL_WORLD, &BP_WORLD);
    GlobalV::NPROC_IN_POOL = nproc;
#else
	MPI_Init(&argc, &argv);	
#endif

    testing::InitGoogleTest(&argc, argv);
    ::testing::TestEventListeners &listeners = ::testing::UnitTest::GetInstance()->listeners();
    if (myrank != 0) { delete listeners.Release(listeners.default_result_printer());
}

    int result = RUN_ALL_TESTS();
    if (myrank == 0 && result != 0)
    {
        std::cout << "ERROR:some tests are not passed" << std::endl;
        return result;
	}

    MPI_Finalize();
	return 0;
}
