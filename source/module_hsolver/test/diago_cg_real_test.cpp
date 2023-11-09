#include "module_base/inverse_matrix.h"
#include "module_base/lapack_connector.h"
#include "module_hamilt_pw/hamilt_pwdft/structure_factor.h"
#include "module_psi/psi.h"
#include "module_hamilt_general/hamilt.h"
#include "module_hamilt_pw/hamilt_pwdft/hamilt_pw.h"
#include "../diago_cg.h"
#include "../diago_iter_assist.h"
#include "diago_mock.h"
#include "mpi.h"
#include "module_basis/module_pw/test/test_tool.h"
#include <complex>

#include "gtest/gtest.h"
#include <random>

/************************************************
 *  unit test of functions in Diago_CG
 ***********************************************/

 /**
  * Class Diago_CG is an approach for eigenvalue problems
  * This unittest test the function Diago_CG::diag() for FPTYPE=double, Device=cpu
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

  // call lapack in order to compare to cg
void lapackEigen(int& npw, std::vector<double>& hm, double* e, bool outtime = false)
{
    int info = 0;
    auto tmp = hm;
    clock_t start, end;
    start = clock();
    char tmp_c1 = 'V', tmp_c2 = 'U';

    double work_tmp;
    constexpr int minus_one = -1;
    dsyev_(&tmp_c1, &tmp_c2, &npw, tmp.data(), &npw, e, &work_tmp, &minus_one, &info);		// get best lwork

    const int lwork = work_tmp;
    double* work2 = new double[lwork];
    dsyev_(&tmp_c1, &tmp_c2, &npw, tmp.data(), &npw, e, work2, &lwork, &info);
    end = clock();
    if (info) std::cout << "ERROR: Lapack solver, info=" << info << std::endl;
    if (outtime) std::cout << "Lapack Run time: " << (double)(end - start) / CLOCKS_PER_SEC << " S" << std::endl;
    delete[] work2;
}

class DiagoCGPrepare
{
public:
    DiagoCGPrepare(int nband, int npw, int sparsity, bool reorder, double eps, int maxiter, double threshold)
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
    int nprocs = 1, mypnum = 0;
    // threshold is the comparison standard between cg and lapack

    void CompareEigen(double* precondition)
    {
        // calculate eigenvalues by LAPACK;
        double* e_lapack = new double[npw];
        auto ev = DIAGOTEST::hmatrix_d;
        if (mypnum == 0)  lapackEigen(npw, ev, e_lapack, false);
        // initial guess of psi by perturbing lapack psi
        ModuleBase::matrix psiguess(nband, npw);
        std::default_random_engine p(1);
        std::uniform_int_distribution<unsigned> u(1, 10);
        for (int i = 0; i < nband; i++)
        {
            for (int j = 0; j < npw; j++)
            {
                double rand = static_cast<double>(u(p)) / 10.;
                // psiguess(i,j) = ev(j,i)*(1+rand);
                psiguess(i, j) = ev[j * DIAGOTEST::h_nc + i] * rand;
            }
        }
        // run cg
    //======================================================================
        double* en = new double[npw];
        int ik = 1;
        hamilt::Hamilt<double>* ha;
        ha = new hamilt::HamiltPW<double>(nullptr, nullptr, nullptr);
        int* ngk = new int[1];
        psi::Psi<double> psi;
        psi.resize(ik, nband, npw);
        //psi.fix_k(0);
        for (int i = 0; i < nband; i++)
        {
            for (int j = 0; j < npw; j++)
            {
                psi(i, j) = psiguess(i, j);
            }
        }

        psi::Psi<double> psi_local;
        double* precondition_local;
        DIAGOTEST::npw_local = new int[nprocs];
#ifdef __MPI				
        DIAGOTEST::cal_division(DIAGOTEST::npw);
        DIAGOTEST::divide_hpsi(psi, psi_local, DIAGOTEST::hmatrix_d, DIAGOTEST::hmatrix_local_d); //will distribute psi and Hmatrix to each process
        precondition_local = new double[DIAGOTEST::npw_local[mypnum]];
        DIAGOTEST::divide_psi<double>(precondition, precondition_local);
#else
        DIAGOTEST::hmatrix_local_d = DIAGOTEST::hmatrix_d;
        DIAGOTEST::npw_local[0] = DIAGOTEST::npw;
        psi_local = psi;
        precondition_local = new double[DIAGOTEST::npw];
        for (int i = 0;i < DIAGOTEST::npw;i++) precondition_local[i] = precondition[i];
#endif
        hsolver::DiagoCG<double> cg(precondition_local);
        psi_local.fix_k(0);
        double start, end;
        start = MPI_Wtime();
        cg.diag(ha, psi_local, en);
        end = MPI_Wtime();
        //if(mypnum == 0) printf("diago time:%7.3f\n",end-start);
        delete[] DIAGOTEST::npw_local;
        delete[] precondition_local;
        //======================================================================
        for (int i = 0; i < nband; i++)
        {
            EXPECT_NEAR(en[i], e_lapack[i], threshold);
        }

        delete[] en;
        delete[] e_lapack;
        delete ha;
    }
};

class DiagoCGTest : public ::testing::TestWithParam<DiagoCGPrepare>
{
};

TEST_P(DiagoCGTest, RandomHamilt)
{
    DiagoCGPrepare dcp = GetParam();
    hsolver::DiagoIterAssist<double>::PW_DIAG_NMAX = dcp.maxiter;
    hsolver::DiagoIterAssist<double>::PW_DIAG_THR = dcp.eps;
    HPsi<double> hpsi(dcp.nband, dcp.npw, dcp.sparsity);
    DIAGOTEST::hmatrix_d = hpsi.hamilt();

    DIAGOTEST::npw = dcp.npw;
    // ModuleBase::matrix psi = hpsi.psi();
    dcp.CompareEigen(hpsi.precond());
}

INSTANTIATE_TEST_SUITE_P(VerifyCG,
    DiagoCGTest,
    ::testing::Values(
        // nband, npw, sparsity, reorder, eps, maxiter, threshold
        DiagoCGPrepare(10, 500, 0, true, 1e-5, 300, 1e-3),
        DiagoCGPrepare(20, 500, 6, true, 1e-5, 300, 1e-3),
        DiagoCGPrepare(20, 1000, 8, true, 1e-5, 300, 1e-3),
        DiagoCGPrepare(40, 1000, 8, true, 1e-6, 300, 1e-3)));
//DiagoCGPrepare(40, 2000, 8, true, 1e-5, 500, 1e-2))); 
// the last one is passed but time-consumming.

// check that the mock class HPsi work well
// in generating a Hermite matrix
TEST(DiagoCGTest, Hamilt)
{
    int dim = 2;
    int nbnd = 2;
    HPsi<double> hpsi(nbnd, dim);
    std::vector<double> hm = hpsi.hamilt();
    EXPECT_EQ(DIAGOTEST::h_nr, 2);
    EXPECT_EQ(DIAGOTEST::h_nc, 2);
    EXPECT_EQ(hm[DIAGOTEST::h_nc], hm[1]);
}

// cg for a 2x2 matrix
#ifdef __MPI
#else
TEST(DiagoCGTest, TwoByTwo)
{
    int dim = 2;
    int nband = 2;
    ModuleBase::matrix hm(2, 2);
    hm(0, 0) = 4.0;
    hm(0, 1) = 1.0;
    hm(1, 0) = 1.0;
    hm(1, 1) = 3.0;
    // nband, npw, sub, sparsity, reorder, eps, maxiter, threshold
    DiagoCGPrepare dcp(nband, dim, 0, true, 1e-4, 50, 1e-10);
    hsolver::DiagoIterAssist<double>::PW_DIAG_NMAX = dcp.maxiter;
    hsolver::DiagoIterAssist<double>::PW_DIAG_THR = dcp.eps;
    HPsi<double> hpsi;
    hpsi.create(nband, dim);
    DIAGOTEST::hmatrix_d = hm;
    DIAGOTEST::npw = dim;
    dcp.CompareEigen(hpsi.precond());
}
#endif

TEST(DiagoCGTest, readH)
{
    // read Hamilt matrix from file data-H
    std::vector<double> hm;
    std::ifstream ifs;
    ifs.open("H-GammaOnly-Si64.dat");
    DIAGOTEST::readh(ifs, hm);
    ifs.close();
    int dim = DIAGOTEST::npw;
    int nband = 10; // not nband < dim, here dim = 26 in data-H
    // nband, npw, sub, sparsity, reorder, eps, maxiter, threshold
    DiagoCGPrepare dcp(nband, dim, 0, true, 1e-5, 500, 1e-3);
    hsolver::DiagoIterAssist<double>::PW_DIAG_NMAX = dcp.maxiter;
    hsolver::DiagoIterAssist<double>::PW_DIAG_THR = dcp.eps;
    HPsi<double> hpsi;
    hpsi.create(nband, dim);
    DIAGOTEST::hmatrix_d = hpsi.hamilt();
    DIAGOTEST::npw = dim;
    dcp.CompareEigen(hpsi.precond());
}

int main(int argc, char** argv)
{
    int nproc = 1, myrank = 0;

#ifdef __MPI
    int nproc_in_pool, kpar = 1, mypool, rank_in_pool;
    setupmpi(argc, argv, nproc, myrank);
    divide_pools(nproc, myrank, nproc_in_pool, kpar, mypool, rank_in_pool);
    GlobalV::NPROC_IN_POOL = nproc;
#else
    MPI_Init(&argc, &argv);
#endif

    testing::InitGoogleTest(&argc, argv);
    ::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();
    if (myrank != 0) delete listeners.Release(listeners.default_result_printer());

    int result = RUN_ALL_TESTS();
    if (myrank == 0 && result != 0)
    {
        std::cout << "ERROR:some tests are not passed" << std::endl;
        return result;
    }

    MPI_Finalize();
    return 0;
}
