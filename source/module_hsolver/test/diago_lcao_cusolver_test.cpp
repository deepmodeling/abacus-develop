#include <vector>

#include "gtest/gtest.h"
#include "module_hsolver/diago_blas.h"
#include "module_hsolver/test/diago_elpa_utils.h"
#include "mpi.h"
#include "string.h"
#ifdef __ELPA
#include "module_hsolver/diago_elpa.h"
#endif
#ifdef __CUSOLVER_LCAO
#include "module_hsolver/diago_cusolver.h"
#endif

#define PASSTHRESHOLD 1e-10
#define DETAILINFO false
#define PRINT_HS false
#define REPEATRUN 1

/************************************************
 *  unit test of LCAO diagonalization
 ***********************************************/

/**
 * Tested function:
 *  - hsolver::DiagoElpa::diag (for ELPA)
 *  - hsolver::DiagoBlas::diag (for Scalapack)
 *
 * The 2d block cyclic distribution of H/S matrix is done by
 * self-realized functions in module_hsolver/test/diago_elpa_utils.h
 */

template <typename T>
class HamiltTEST : public hamilt::Hamilt<T>
{
  public:
    int desc[9];
    int nrow, ncol;
    std::vector<T> h_local;
    std::vector<T> s_local;

    void matrix(hamilt::MatrixBlock<T>& hk_in, hamilt::MatrixBlock<T>& sk_in)
    {
        hk_in = hamilt::MatrixBlock<T>{this->h_local.data(), (size_t)this->nrow, (size_t)this->ncol, this->desc};
        sk_in = hamilt::MatrixBlock<T>{this->s_local.data(), (size_t)this->nrow, (size_t)this->ncol, this->desc};
    }

    void constructHamilt(const int iter, const hamilt::MatrixBlock<double> rho)
    {
    }
    void updateHk(const int ik)
    {
    }
};

template <class T>
class DiagoPrepare
{
  public:
    DiagoPrepare(int nlocal,
                 int nbands,
                 int nb2d,
                 int sparsity,
                 std::string ks_solver,
                 std::string hfname,
                 std::string sfname)
        : nlocal(nlocal),
          nbands(nbands),
          nb2d(nb2d),
          sparsity(sparsity),
          ks_solver(ks_solver),
          hfname(hfname),
          sfname(sfname)
    {
        MPI_Comm_size(MPI_COMM_WORLD, &dsize);
        MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

        if (ks_solver == "scalapack_gvx")
            dh = new hsolver::DiagoBlas<T>;
#ifdef __CUSOLVER_LCAO
        else if (ks_solver == "cusolver")
            dh = new hsolver::DiagoCusolver<T>;
#endif
        else
        {
            if (myrank == 0)
                std::cout << "ERROR: undefined ks_solver: " << ks_solver << std::endl;
            exit(1);
        }
    }

    int dsize, myrank;
    int nlocal, nbands, nb2d, sparsity;
    double hsolver_time = 0.0, lapack_time = 0.0;
    std::string ks_solver, sfname, hfname;
    HamiltTEST<T> hmtest;
    std::vector<T> h;
    std::vector<T> s;
    std::vector<T> h_local;
    std::vector<T> s_local;
    hsolver::DiagH<T>* dh = 0;
    psi::Psi<T> psi;
    std::vector<double> e_solver;
    std::vector<double> e_lapack;
    std::vector<double> abc;
    int icontxt;

    bool read_HS()
    {
        bool readhfile = false;
        bool readsfile = false;
        if (this->myrank == 0)
        {
            int hdim, sdim;
            readhfile = LCAO_DIAGO_TEST::read_hs<std::vector<T>>(hfname, this->h);
            readsfile = LCAO_DIAGO_TEST::read_hs<std::vector<T>>(sfname, this->s);
            hdim = sqrt(this->h.size());
            sdim = sqrt(this->s.size());
            if (hdim != sdim)
            {
                printf("Error: dimensions of H and S are not equal, %d, %d\n", hdim, sdim);
                readhfile = readsfile = false;
            }
            nlocal = hdim;
        }
        MPI_Bcast(&nlocal, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&readhfile, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
        MPI_Bcast(&readsfile, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
        nbands = nlocal / 2;
        if (readhfile && readsfile)
            return true;
        return false;
    }

    bool produce_HS()
    {
        bool ok = this->read_HS();

        e_solver.resize(nlocal, 0.0);
        e_lapack.resize(nlocal, 1.0);

        return ok;
    }

    void print_hs()
    {
        if (!PRINT_HS)
            return;
        if (myrank == 0)
        {
            std::ofstream fp("hmatrix.dat");
            LCAO_DIAGO_TEST::print_matrix(fp, this->h.data(), nlocal, nlocal, true);
            fp.close();
            fp.open("smatrix.dat");
            LCAO_DIAGO_TEST::print_matrix(fp, this->s.data(), nlocal, nlocal, true);
            fp.close();
        }
    }

    void pb2d()
    {
        int nprows, npcols, myprow, mypcol;
        LCAO_DIAGO_TEST::process_2d(nprows, npcols, myprow, mypcol, icontxt);

        hmtest.nrow = LCAO_DIAGO_TEST::na_rc(nlocal,
                                             nb2d,
                                             nprows,
                                             myprow); // the number of row of the new_matrix in each process
        hmtest.ncol = LCAO_DIAGO_TEST::na_rc(nlocal,
                                             nb2d,
                                             npcols,
                                             mypcol); // the number of column of the new_matrix in each process
        int ISRC = 0, info;
        descinit_(hmtest.desc, &nlocal, &nlocal, &nb2d, &nb2d, &ISRC, &ISRC, &icontxt, &(hmtest.nrow), &info);
        if (info != 0)
        {
            printf("Invalid blacs-distribution. Abort!\n");
            exit(1);
        }
    }

    void distribute_data()
    {

        int local_size = hmtest.nrow * hmtest.ncol;
        this->h_local.resize(local_size);
        this->s_local.resize(local_size);

        LCAO_DIAGO_TEST::distribute_data<T>(this->h.data(),
                                            this->h_local.data(),
                                            nlocal,
                                            nb2d,
                                            hmtest.nrow,
                                            hmtest.ncol,
                                            icontxt);
        LCAO_DIAGO_TEST::distribute_data<T>(this->s.data(),
                                            this->s_local.data(),
                                            nlocal,
                                            nb2d,
                                            hmtest.nrow,
                                            hmtest.ncol,
                                            icontxt);
        psi.resize(1, hmtest.ncol, hmtest.nrow);
    }

    void set_env()
    {
        GlobalV::NLOCAL = nlocal;
        GlobalV::NBANDS = nbands;
        GlobalV::DSIZE = dsize;
    }

    void diago()
    {
        this->pb2d();
        this->distribute_data();
        this->print_hs();
        this->set_env();

        double starttime = 0.0, endtime = 0.0;
        MPI_Barrier(MPI_COMM_WORLD);
        starttime = MPI_Wtime();
        for (int i = 0; i < REPEATRUN; i++)
        {
            hmtest.h_local = this->h_local;
            hmtest.s_local = this->s_local;
            dh->diag(&hmtest, psi, e_solver.data());
        }
        endtime = MPI_Wtime();
        hsolver_time = (endtime - starttime) / REPEATRUN;
        delete dh;
    }

    void diago_lapack()
    {
        double starttime = 0.0, endtime = 0.0;
        starttime = MPI_Wtime();
        for (int i = 0; i < REPEATRUN; i++)
            LCAO_DIAGO_TEST::lapack_diago(this->h.data(), this->s.data(), this->e_lapack.data(), nlocal);
        endtime = MPI_Wtime();
        lapack_time = (endtime - starttime) / REPEATRUN;
    }

    bool compare_eigen(std::stringstream& out_info)
    {
        double maxerror = 0.0;
        int iindex = 0;
        bool pass = true;
        for (int i = 0; i < nbands; i++)
        {
            // EXPECT_NEAR(e_lapack[i], e_solver[i], PASSTHRESHOLD);
            double error = std::abs(e_lapack[i] - e_solver[i]);
            if (error > maxerror)
            {
                maxerror = error;
                iindex = i;
            }
            if (error > PASSTHRESHOLD)
                pass = false;
        }

        std::cout << "H/S matrix are read from " << hfname << ", " << sfname << std::endl;
        std::cout << "solver=" << ks_solver << ", NLOCAL=" << nlocal << ", nbands=" << nbands << ", nb2d=" << nb2d;
        std::cout << std::endl;
        out_info << "solver time: " << hsolver_time << "s, LAPACK time(1 core):" << lapack_time << "s" << std::endl;
        out_info << "Maximum difference between ks_hsolver and LAPACK is " << maxerror << " (" << iindex
                 << "-th eigenvalue), the pass threshold is " << PASSTHRESHOLD << std::endl;

        if (DETAILINFO)
        {
            std::cout << out_info.str();
            out_info.str("");
            out_info.clear();
        }
        return pass;
    }
};

class DiagoGammaOnlyTest : public ::testing::TestWithParam<DiagoPrepare<double>>
{
};

TEST_P(DiagoGammaOnlyTest, LCAO)
{
    std::stringstream out_info;
    DiagoPrepare<double> dp = GetParam();
    ASSERT_TRUE(dp.produce_HS());
    dp.diago();

    if (dp.myrank == 0)
    {
        dp.diago_lapack();
        bool pass = dp.compare_eigen(out_info);
        EXPECT_TRUE(pass) << out_info.str();
    }

    MPI_Barrier(MPI_COMM_WORLD);
}

INSTANTIATE_TEST_SUITE_P(
    DiagoTest,
    DiagoGammaOnlyTest,
    ::testing::Values( // int nlocal, int nbands, int nb2d, int sparsity, std::string ks_solver_in, std::string hfname,
                       // std::string sfname DiagoPrepare<double>(0, 0, 1, 0, "genelpa", "H-GammaOnly-Si2.dat",
                       // "S-GammaOnly-Si2.dat")
#ifdef __CUSOLVER_LCAO
        DiagoPrepare<double>(0, 0, 32, 0, "cusolver", "H-GammaOnly-Si64.dat", "S-GammaOnly-Si64.dat"),
#endif
        DiagoPrepare<double>(0, 0, 1, 0, "scalapack_gvx", "H-GammaOnly-Si2.dat", "S-GammaOnly-Si2.dat"),
        DiagoPrepare<double>(0, 0, 32, 0, "scalapack_gvx", "H-GammaOnly-Si64.dat", "S-GammaOnly-Si64.dat")));

class DiagoKPointsTest : public ::testing::TestWithParam<DiagoPrepare<std::complex<double>>>
{
};
TEST_P(DiagoKPointsTest, LCAO)
{
    std::stringstream out_info;
    DiagoPrepare<std::complex<double>> dp = GetParam();
    ASSERT_TRUE(dp.produce_HS());
    dp.diago();

    if (dp.myrank == 0)
    {
        dp.diago_lapack();
        bool pass = dp.compare_eigen(out_info);
        EXPECT_TRUE(pass) << out_info.str();
    }
    MPI_Barrier(MPI_COMM_WORLD);
}
INSTANTIATE_TEST_SUITE_P(
    DiagoTest,
    DiagoKPointsTest,
    ::testing::Values( // int nlocal, int nbands, int nb2d, int sparsity, std::string ks_solver_in, std::string hfname,
                       // std::string sfname DiagoPrepare<std::complex<double>>(800, 400, 32, 7, "genelpa", "", ""),
#ifdef __CUSOLVER_LCAO
        DiagoPrepare<std::complex<double>>(0, 0, 1, 0, "cusolver", "H-KPoints-Si2.dat", "S-KPoints-Si2.dat"),
#endif
        // DiagoPrepare<std::complex<double>>(0, 0, 32, 0, "genelpa", "H-KPoints-Si64.dat", "S-KPoints-Si64.dat"),
        DiagoPrepare<std::complex<double>>(0, 0, 1, 0, "scalapack_gvx", "H-KPoints-Si2.dat", "S-KPoints-Si2.dat"),
        DiagoPrepare<std::complex<double>>(0, 0, 32, 0, "scalapack_gvx", "H-KPoints-Si64.dat", "S-KPoints-Si64.dat")));

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int mypnum, dsize;
    MPI_Comm_size(MPI_COMM_WORLD, &dsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mypnum);

    testing::InitGoogleTest(&argc, argv);
    // Parallel_Global::split_diag_world(dsize);
    ::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();
    if (mypnum != 0)
    {
        delete listeners.Release(listeners.default_result_printer());
    }
    int result = RUN_ALL_TESTS();

    if (mypnum == 0 && result != 0)
    {
        std::cout << "ERROR:some tests are not passed" << std::endl;
        return result;
    }
    else
    {
        MPI_Finalize();
        return 0;
    }
}