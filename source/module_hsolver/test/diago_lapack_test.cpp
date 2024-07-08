// Author: Zhang Xiaoyang
// A modified version of diago_lcao_test.cpp
// Remove some useless functions and dependencies. Serialized the full code
// and refactored some function.

#include "gtest/gtest.h"
#include "string.h"
#include <vector>

#include "module_hsolver/diago_lapack.h"
#include "module_base/lapack_connector.h"

#define PASSTHRESHOLD 1e-10
#define DETAILINFO false
#define PRINT_HS false
#define REPEATRUN 1

// A hamilt class used for test. It will be removed in the future.

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
    if (outtime)
        std::cout << "Lapack Run time: " << (double)(end - start) / CLOCKS_PER_SEC << " S" << std::endl;
    delete[] rwork;
    delete[] work2;
}


void lapack_diago(double *hmatrix, double *smatrix, double *e, int &nFull)
{
    const int itype = 1; // solve A*X=(lambda)*B*X
    const char jobz = 'V'; // 'N':only calc eigenvalue, 'V': eigenvalues and eigenvectors
    const char uplo = 'U'; // Upper triangles
    int lwork = (nFull + 2) * nFull, info = 0;
    double *ev = new double[nFull * nFull];

    double *a = new double[nFull * nFull];
    double *b = new double[nFull * nFull];
    for (int i = 0; i < nFull * nFull; i++)
    {
        a[i] = hmatrix[i];
        b[i] = smatrix[i];
    }

    dsygv_(&itype, &jobz, &uplo, &nFull, a, &nFull, b, &nFull, e, ev, &lwork, &info);
    if (info != 0)
    {
        std::cout << "ERROR: solvered by LAPACK error, info=" << info << std::endl;
        exit(1);
    }

    delete[] a;
    delete[] b;
    delete[] ev;
}

void lapack_diago(std::complex<double> *hmatrix, std::complex<double> *smatrix, double *e, int &nFull)
{
    const int itype = 1; // solve A*X=(lambda)*B*X
    const char jobz = 'V'; // 'N':only calc eigenvalue, 'V': eigenvalues and eigenvectors
    const char uplo = 'U'; // Upper triangles
    int lwork = (nFull + 1) * nFull, info = 0;
    double *rwork = new double[3 * nFull - 2];
    std::complex<double> *ev = new std::complex<double>[nFull * nFull];

    std::complex<double> *a = new std::complex<double>[nFull * nFull];
    std::complex<double> *b = new std::complex<double>[nFull * nFull];
    for (int i = 0; i < nFull * nFull; i++)
    {
        a[i] = hmatrix[i];
        b[i] = smatrix[i];
    }

    zhegv_(&itype, &jobz, &uplo, &nFull, a, &nFull, b, &nFull, e, ev, &lwork, rwork, &info);
    if (info != 0)
    {
        std::cout << "ERROR: solvered by LAPACK error, info=" << info << std::endl;
        exit(1);
    }

    delete[] a;
    delete[] b;
    delete[] ev;
    delete[] rwork;
}

// The serialized version of functions from diago_elpa_utils

template <class T> bool read_hs(std::string fname, T &matrix)
{
    int ndim;
    std::ifstream inf(fname);
    if(! inf.is_open())
    {
        std::cout << "Error: open file " << fname << " failed, skip!" << std::endl;
        return false;
    }
    inf >> ndim;
    matrix.resize(ndim*ndim);
    for (int i = 0; i < ndim; i++)
    {
        for (int j = i; j < ndim; j++)
        {
            inf >> matrix[i * ndim + j];
            if (i != j)
            {
                matrix[j * ndim + i] = matrix[i * ndim + j];
            }
        }
    }
    inf.close();
    return true;
}

template <class T> inline void print_matrix(std::ofstream &fp, T *matrix, int &nrow, int &ncol, bool row_first)
{
    int coef_row = row_first ? ncol : 1;
    int coef_col = row_first ? 1 : nrow;
    for (int i = 0; i < nrow; i++)
    {
        for (int j = 0; j < ncol; j++)
        {
            fp << std::setw(15) << matrix[i * coef_row + j * coef_col] << " ";
        }
        fp << std::endl;
    }
}


template <class T>
class DiagoLapackPrepare
{
    public:
    DiagoLapackPrepare(int nlocal,
                       int nbands,
                       int nb2d,
                       int sparsity,
                       std::string hfname,
                       std::string sfname)
        : nlocal(nlocal), nbands(nbands), nb2d(nb2d), sparsity(sparsity), hfname(hfname),
          sfname(sfname)
    {
        dh = new hsolver::DiagoLapack<T>;
    }

    int nlocal, nbands, nb2d, sparsity;
    std::string sfname, hfname;
    std::vector<T> h;
    std::vector<T> s;
    HamiltTEST<T> hmtest;
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
        int hdim, sdim;
        readhfile = read_hs<std::vector<T>>(hfname, this->h);
        readsfile = read_hs<std::vector<T>>(sfname, this->s);
        hdim = sqrt(this->h.size());
        sdim = sqrt(this->s.size());
        if (hdim != sdim)
        {
            printf("Error: dimensions of H and S are not equal, %d, %d\n", hdim, sdim);
            readhfile = readsfile = false;
        }
        nlocal = hdim;
        nbands = nlocal / 2;
        if (readhfile && readsfile)
            std::cout << "READ FINISH" << std::endl;
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
        std::ofstream fp("hmatrix.dat");
        print_matrix(fp, this->h.data(), nlocal, nlocal, true);
        fp.close();
        fp.open("smatrix.dat");
        print_matrix(fp, this->s.data(), nlocal, nlocal, true);
        fp.close();
    }

    void pb2d()
    {
        hmtest.h_local = this->h;
        hmtest.s_local = this->s;

        hmtest.nrow = nlocal;
        hmtest.ncol = nlocal;
    }

    void set_env()
    {
        GlobalV::NLOCAL = nlocal;
        GlobalV::NBANDS = nbands;
    }

    void diago()
    {
        this->pb2d();
        this->print_hs();
        this->set_env();

        for (int i = 0; i < REPEATRUN; i++)
        {
            dh->diag(&hmtest, psi, e_solver.data());
        }
        delete dh;
    }

    void diago_lapack()
    {
        for (int i = 0; i < REPEATRUN; i++)
            lapack_diago(this->h.data(), this->s.data(), this->e_lapack.data(), nlocal);
    }

    bool compare_eigen(std::stringstream& out_info)
    {
        double maxerror = 0.0;
        int iindex = 0;
        bool pass = true;
        for (int i = 0; i < nbands; i++)
        {
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
        std::cout << ", NLOCAL=" << nlocal << ", nbands=" << nbands << ", nb2d=" << nb2d;
        std::cout << std::endl;
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

class DiagoLapackGammaOnlyTest : public ::testing::TestWithParam<DiagoLapackPrepare<double>>
{
};

TEST_P(DiagoLapackGammaOnlyTest, LCAO)
{
    std::stringstream out_info;
    DiagoLapackPrepare<double> dp = GetParam();
    ASSERT_TRUE(dp.produce_HS());
    dp.diago();

    dp.diago_lapack();
    bool pass = dp.compare_eigen(out_info);
    EXPECT_TRUE(pass) << out_info.str();
}

INSTANTIATE_TEST_SUITE_P(
    DiagoLapackTest,
    DiagoLapackGammaOnlyTest,
    ::testing::Values( // int nlocal, int nbands, int nb2d, int sparsity, std::string ks_solver_in, std::string hfname,
                       // std::string sfname DiagoLapackPrepare<double>(0, 0, 1, 0, "genelpa", "H-GammaOnly-Si2.dat",
                       // "S-GammaOnly-Si2.dat")
        DiagoLapackPrepare<double>(0, 0, 1, 0, "H-GammaOnly-Si2.dat", "S-GammaOnly-Si2.dat"),
        DiagoLapackPrepare<double>(0, 0, 32, 0, "H-GammaOnly-Si64.dat", "S-GammaOnly-Si64.dat")));


class DiagoLapackKPointsTest : public ::testing::TestWithParam<DiagoLapackPrepare<std::complex<double>>>
{
};
TEST_P(DiagoLapackKPointsTest, LCAO)
{
    std::stringstream out_info;
    DiagoLapackPrepare<std::complex<double>> dp = GetParam();
    ASSERT_TRUE(dp.produce_HS());
    dp.diago();

    dp.diago_lapack();
    bool pass = dp.compare_eigen(out_info);
    EXPECT_TRUE(pass) << out_info.str();
}

INSTANTIATE_TEST_SUITE_P(
    DiagoLapackTest,
    DiagoLapackKPointsTest,
    ::testing::Values( // int nlocal, int nbands, int nb2d, int sparsity, std::string ks_solver_in, std::string hfname,
                       // std::string sfname DiagoLapackPrepare<std::complex<double>>(800, 400, 32, 7, "genelpa", "", ""),
        DiagoLapackPrepare<std::complex<double>>(0, 0, 1, 0, "H-KPoints-Si2.dat", "S-KPoints-Si2.dat"),
        DiagoLapackPrepare<std::complex<double>>(0, 0, 32, 0, "H-KPoints-Si64.dat", "S-KPoints-Si64.dat")));

int main(int argc, char** argv)
{
    testing::InitGoogleTest(&argc, argv);
    ::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();
    delete listeners.Release(listeners.default_result_printer());
    int result = RUN_ALL_TESTS();
    if (result != 0)
    {
        std::cout << "ERROR:some tests are not passed" << std::endl;
        return result;
    }
    else
    {
        return 0;
    }
}