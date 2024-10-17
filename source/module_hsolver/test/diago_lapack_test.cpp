// Author: Zhang Xiaoyang
// A modified version of diago_lcao_test.cpp
#define private public
#include "module_parameter/parameter.h"
#undef private
// Remove some useless functions and dependencies. Serialized the full code
// and refactored some function.

#include "gtest/gtest.h"
#include <cstring>
#include <vector>

#include "module_hsolver/diago_lapack.h"
#include "module_base/lapack_connector.h"

#define PASSTHRESHOLD 1e-5
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

bool read_solution(std::string fname, std::vector<double> &result)
{
    int ndim;
    std::ifstream inf(fname);
    if(! inf.is_open())
    {
        std::cout << "Error: open file " << fname << " failed, skip!" << std::endl;
        return false;
    }
    inf >> ndim;
    for (int i = 0; i < ndim; i++)
    {
        inf >> result[i];
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
                       std::string sfname,
                       std::string solutionfname)
        : nlocal(nlocal), nbands(nbands), nb2d(nb2d), sparsity(sparsity), hfname(hfname),
          sfname(sfname), solutionfname(solutionfname)
    {
        // dh = new hsolver::DiagoLapack<T>;
    }

    int nlocal, nbands, nb2d, sparsity;
    std::string sfname, hfname, solutionfname;
    std::vector<T> h;
    std::vector<T> s;
    HamiltTEST<T> hmtest;
    // hsolver::DiagH<T>* dh = nullptr;
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
        if (readhfile && readsfile) {
            return true;
        }  
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
        if (!PRINT_HS) {
            return;
        }
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
        PARAM.sys.nlocal = nlocal;
        PARAM.input.nbands = nbands;
    }

    void diago()
    {
        this->pb2d();
        this->print_hs();
        this->set_env();

        for (int i = 0; i < REPEATRUN; i++)
        {
            hsolver::DiagoLapack<T> dh;
            dh.diag(&hmtest, psi, e_solver.data());
            // dh->diag(&hmtest, psi, e_solver.data());
        }
        // delete dh;
    }

    void read_SOLUTION()
    {
        read_solution(this->solutionfname, this->e_lapack);
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
            if (error > PASSTHRESHOLD) {
                pass = false;
            }
        }
        std::cout << "H/S matrix are read from " << hfname << ", " << sfname << std::endl;
        std::cout << ", NLOCAL=" << nlocal << ", nbands=" << nbands << ", nb2d=" << nb2d;
        std::cout << std::endl;
        std::cout << "Solution is read from " << solutionfname << std::endl << std::endl;
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

    dp.read_SOLUTION();
    bool pass = dp.compare_eigen(out_info);
    EXPECT_TRUE(pass) << out_info.str();
}

INSTANTIATE_TEST_SUITE_P(
    DiagoLapackTest,
    DiagoLapackGammaOnlyTest,
    ::testing::Values(
        DiagoLapackPrepare<double>(0, 0, 1, 0, "H-GammaOnly-Si2.dat", "S-GammaOnly-Si2.dat", "GammaOnly-Si2-Solution.dat"),
        DiagoLapackPrepare<double>(0, 0, 32, 0, "H-GammaOnly-Si64.dat", "S-GammaOnly-Si64.dat", "GammaOnly-Si64-Solution.dat")));


class DiagoLapackKPointsTest : public ::testing::TestWithParam<DiagoLapackPrepare<std::complex<double>>>
{
};
TEST_P(DiagoLapackKPointsTest, LCAO)
{
    std::stringstream out_info;
    DiagoLapackPrepare<std::complex<double>> dp = GetParam();
    ASSERT_TRUE(dp.produce_HS());
    dp.diago();

    dp.read_SOLUTION();
    bool pass = dp.compare_eigen(out_info);
    EXPECT_TRUE(pass) << out_info.str();
}

INSTANTIATE_TEST_SUITE_P(
    DiagoLapackTest,
    DiagoLapackKPointsTest,
    ::testing::Values(
        DiagoLapackPrepare<std::complex<double>>(0, 0, 1, 0, "H-KPoints-Si2.dat", "S-KPoints-Si2.dat", "KPoints-Si2-Solution.dat"),
        DiagoLapackPrepare<std::complex<double>>(0, 0, 32, 0, "H-KPoints-Si64.dat", "S-KPoints-Si64.dat", "KPoints-Si64-Solution.dat")));

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