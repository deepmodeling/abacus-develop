
#include "diago_elpa_utils.h"
#include "mpi.h"

#include "gtest/gtest.h"
#include <type_traits>

#define DETAILINFO       false
#define PASSTHRESHOLD    1e-12
#define PRINT_NEW_MATRIX false // print out the H matrix after doing block-cycle distribution
#define PRINT_HS         false

/************************************************
 *  unit test of ELPA
 ***********************************************/

/**
 * Test the diagonalization by ELPA with parallel
 * H/S matrix is complex double or double.
 */

template <class T> class ElpaPrepare
{
  public:
    ElpaPrepare(int n, int nblk, int sparsity) : n(n), nblk(nblk), sparsity(sparsity)
    {
        MPI_Comm_rank(MPI_COMM_WORLD, &mypnum);
    }

    ElpaPrepare()
    {
        MPI_Comm_rank(MPI_COMM_WORLD, &mypnum);
    }

    ~ElpaPrepare()
    {
        delete[] hmatrix;
        delete[] smatrix;
        delete[] e_lapack;
        delete[] e_elpa;
    }

    int n = 10; // dimension of H/S matrix
    int nblk = 1; // dimension of sub-matrix
    int sparsity = 100; // sparsity of H matrix;
    T *hmatrix = nullptr;
    T *smatrix = nullptr;
    double *e_lapack = nullptr; // eigenvalues solved by LAPACK
    double *e_elpa = nullptr;
    double lapack_time = 0.0, elpa_time = 0.0;
    int mypnum = 0;
    std::string sfname, hfname;

    bool read_HS()
    {
        bool readhfile = false;
        bool readsfile = false;
        if (mypnum == 0)
        {
            std::vector<T> htmp, stmp;
            readhfile = LCAO_DIAGO_TEST::read_hs<std::vector<T>>(hfname, htmp);
            readsfile = LCAO_DIAGO_TEST::read_hs<std::vector<T>>(sfname, stmp);
            if (htmp.size() != stmp.size())
            {
                std::cout << "Error: dimensions of H and S are not equal: "
                    << htmp.size() << "," << stmp.size() << std::endl;
                readhfile = readsfile = false;
            }

            hmatrix = new T[htmp.size()];
            smatrix = new T[stmp.size()];

            for (int i = 0; i < htmp.size(); i++)
            {
                hmatrix[i] = htmp[i];
                smatrix[i] = stmp[i];
            }

            n = sqrt(stmp.size());
            nblk = 1;
        }
        MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&nblk, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&readhfile, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
        MPI_Bcast(&readsfile, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
        e_lapack = new double[n];
        e_elpa = new double[n];
        for (int i = 0; i < n; i++)
        {
            e_lapack[i] = 0.0;
            e_elpa[i] = 1.0;
        }
        if (readhfile && readsfile)
            return true;
        else
            return false;
    }

    void new_matrix()
    {
        hmatrix = new T[n * n];
        smatrix = new T[n * n];
        e_lapack = new double[n];
        e_elpa = new double[n];
        for (int i = 0; i < n; i++)
        {
            e_lapack[i] = 0.0;
            e_elpa[i] = 1.0;
        }
    }

    void print_hs()
    {
        if (!PRINT_HS)
        {
            return;
        }
        if (mypnum == 0)
        {
            std::ofstream fp("hmatrix.dat");
            LCAO_DIAGO_TEST::print_matrix(fp, hmatrix, n, n, true);
            fp.close();
            fp.open("smatrix.dat");
            LCAO_DIAGO_TEST::print_matrix(fp, smatrix, n, n, true);
            fp.close();
        }
    }

    void lapack_ev()
    {
        double starttime, endtime;
        starttime = MPI_Wtime();
        LCAO_DIAGO_TEST::lapack_diago(hmatrix, smatrix, e_lapack, n);
        endtime = MPI_Wtime();
        lapack_time = endtime - starttime;
    }
};

class ElpaDoubleTest : public ::testing::TestWithParam<ElpaPrepare<double>>
{
};
class ElpaComplexDoubleTest : public ::testing::TestWithParam<ElpaPrepare<std::complex<double>>>
{
};

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    testing::InitGoogleTest(&argc, argv);

    int mypnum;
    MPI_Comm_rank(MPI_COMM_WORLD, &mypnum);
    ::testing::TestEventListeners &listeners = ::testing::UnitTest::GetInstance()->listeners();
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
