#include"module_hsolver/diago_david.h"
#include"module_hsolver/diago_iter_assist.h"
#include"module_hamilt_pw/hamilt_pwdft/hamilt_pw.h"
#include"diago_mock.h"
#include "module_psi/psi.h"
#include"gtest/gtest.h"
#include "module_base/inverse_matrix.h"
#include "module_base/lapack_connector.h"
#include "module_basis/module_pw/test/test_tool.h"
#include"mpi.h"

#define CONVTHRESHOLD 1e-0
#define DETAILINFO false


/************************************************
*  unit test of class Diago_David
***********************************************/

/**
 * Class Diago_David is used to solve the eigenvalues
 * This unittest test the function Diago_David::diag() for FPTYPE=double and Device=cpu
 * with different examples.
 * 	- the hamilt matrix (npw=100,500,1000) produced by random with sparsity of 90%
 *  - the hamilt matrix (npw=100,500,1000) produced by random with sparsity of 50%
 *  - the hamilt matrix (npw=100,500,1000) produced by random with sparsity of 0%
 *  - the hamilt matrix read from "data-H"
 *
 * The test is passed when the eignvalues are closed to these calculated by LAPACK.
 *
 */

 //use lapack to calcualte eigenvalue of matrix hm
 //NOTE: after finish this function, hm stores the eigen vectors.
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
class DiagoDavPrepare
{
public:
    DiagoDavPrepare(int nband, int npw, int sparsity, int order, double eps, int maxiter) :
        nband(nband), npw(npw), sparsity(sparsity), order(order), eps(eps), maxiter(maxiter)
    {
#ifdef __MPI	
        MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
        MPI_Comm_rank(MPI_COMM_WORLD, &mypnum);
#endif					
    }

    int nband, npw, sparsity, order, maxiter, notconv;
    double eps, avg_iter;
    int nprocs = 1, mypnum = 0;

    void CompareEigen(psi::Psi<double>& phi, double* precondition)
    {
        //calculate eigenvalues by LAPACK;
        double* e_lapack = new double[npw];
        double* ev;
        if (mypnum == 0) lapackEigen(npw, DIAGOTEST::hmatrix_d, e_lapack, DETAILINFO);

        //do Diago_David::diag()
        double* en = new double[npw];
        hamilt::Hamilt<double>* phm;
        phm = new hamilt::HamiltPW<double>(nullptr, nullptr, nullptr);
        hsolver::DiagoDavid<double> dav(precondition);
        hsolver::DiagoDavid<double>::PW_DIAG_NDIM = order;
        hsolver::DiagoIterAssist<double>::PW_DIAG_NMAX = maxiter;
        hsolver::DiagoIterAssist<double>::PW_DIAG_THR = eps;
        GlobalV::NPROC_IN_POOL = nprocs;
        phi.fix_k(0);

        double use_time = 0.0;
#ifdef __MPI		
        double start = 0.0, end = 0.0;
        start = MPI_Wtime();
#else
        clock_t start, end;
        start = clock();
#endif	

        dav.diag(phm, phi, en);

#ifdef __MPI		
        end = MPI_Wtime();
        use_time = end - start;
#else		
        end = clock();
        use_time = (double)(end - start);
#endif		

        if (mypnum == 0)
        {
            if (DETAILINFO) std::cout << "diag Run time: " << use_time << std::endl;
            for (int i = 0;i < nband;i++)
            {
                EXPECT_NEAR(en[i], e_lapack[i], CONVTHRESHOLD);
            }
        }
        delete[] en;
        delete phm;
        delete[] e_lapack;
    }
};

class DiagoDavTest : public ::testing::TestWithParam<DiagoDavPrepare> {};

TEST_P(DiagoDavTest, RandomHamilt)
{
    DiagoDavPrepare ddp = GetParam();
    if (DETAILINFO && ddp.mypnum == 0) std::cout << "npw=" << ddp.npw << ", nband=" << ddp.nband << ", sparsity="
        << ddp.sparsity << ", eps=" << ddp.eps << std::endl;

    HPsi<double> hpsi(ddp.nband, ddp.npw, ddp.sparsity);
    DIAGOTEST::hmatrix_d = hpsi.hamilt();
    DIAGOTEST::npw = ddp.npw;
    DIAGOTEST::npw_local = new int[ddp.nprocs];
    psi::Psi<double> psi = hpsi.psi();
    psi::Psi<double> psi_local;
    double* precondition_local;

#ifdef __MPI				
    DIAGOTEST::cal_division(DIAGOTEST::npw);
    DIAGOTEST::divide_hpsi(psi, psi_local, DIAGOTEST::hmatrix_d, DIAGOTEST::hmatrix_local_d);
    precondition_local = new double[DIAGOTEST::npw_local[ddp.mypnum]];
    DIAGOTEST::divide_psi<double>(hpsi.precond(), precondition_local);
#else
    DIAGOTEST::hmatrix_local_d = DIAGOTEST::hmatrix_d;
    DIAGOTEST::npw_local[0] = DIAGOTEST::npw;
    psi_local = psi;
    precondition_local = new double[DIAGOTEST::npw];
    for (int i = 0;i < DIAGOTEST::npw;i++) precondition_local[i] = (hpsi.precond())[i];
#endif

    ddp.CompareEigen(psi_local, precondition_local);
    delete[] DIAGOTEST::npw_local;
    delete[] precondition_local;
}


INSTANTIATE_TEST_SUITE_P(VerifyDiag, DiagoDavTest, ::testing::Values(
    //DiagoDavPrepare(int nband, int npw, int sparsity, int order,double eps,int maxiter)
    DiagoDavPrepare(10, 100, 0, 4, 1e-5, 500),
    DiagoDavPrepare(20, 500, 7, 4, 1e-5, 500)
    //DiagoDavPrepare(50,1000,8,4,1e-5,500)
    //DiagoDavPrepare(20,2000,8,4,1e-5,500)
));

TEST(DiagoDavRealSystemTest, dataH)
{
    std::vector<double> hmatrix;
    std::ifstream ifs("H-GammaOnly-Si64.dat");
    DIAGOTEST::readh(ifs, hmatrix);
    ifs.close();
    DIAGOTEST::hmatrix_d = hmatrix;
    int nband = std::max(DIAGOTEST::npw / 20, 1);

    DiagoDavPrepare ddp(nband, DIAGOTEST::npw, 0, 2, 1e-5, 500);

    HPsi<double> hpsi(nband, DIAGOTEST::npw);
    psi::Psi<double> psi = hpsi.psi();
    DIAGOTEST::npw_local = new int[ddp.nprocs];
    psi::Psi<double> psi_local;
    double* precondition_local;

#ifdef __MPI				
    DIAGOTEST::cal_division(DIAGOTEST::npw);
    DIAGOTEST::divide_hpsi(psi, psi_local, DIAGOTEST::hmatrix_d, DIAGOTEST::hmatrix_local_d);
    precondition_local = new double[DIAGOTEST::npw_local[ddp.mypnum]];
    DIAGOTEST::divide_psi<double>(hpsi.precond(), precondition_local);
#else
    DIAGOTEST::hmatrix_local_d = DIAGOTEST::hmatrix_d;
    DIAGOTEST::npw_local[0] = DIAGOTEST::npw;
    psi_local = psi;
    precondition_local = new double[DIAGOTEST::npw];
    for (int i = 0;i < DIAGOTEST::npw;i++) precondition_local[i] = (hpsi.precond())[i];
#endif

    ddp.CompareEigen(psi_local, precondition_local);

    delete[] DIAGOTEST::npw_local;
    delete[] precondition_local;
}

int main(int argc, char** argv)
{
    int nproc = 1, myrank = 0;

#ifdef __MPI
    int nproc_in_pool, kpar = 1, mypool, rank_in_pool;
    setupmpi(argc, argv, nproc, myrank);
    divide_pools(nproc, myrank, nproc_in_pool, kpar, mypool, rank_in_pool);
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
