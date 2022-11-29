#include "module_base/blacs_connector.h"
#include "module_base/lapack_connector.h"
#include "module_base/scalapack_connector.h"
#include "mpi.h"

#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <random>
#include <type_traits>
#include <vector>

#ifdef __ELPA
extern "C"
{
#include "module_hsolver/my_elpa.h"
}
#endif

namespace LCAO_DIAGO_TEST
{
void process_2d(int &nprows, int &npcols, int &myprow, int &mypcol, int &icontxt) // map the processes to 2D matrix
{
    int nprocs, mypnum;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &mypnum);

    for (npcols = int(sqrt(double(nprocs))); npcols >= 1; --npcols)
    {
        if (nprocs % npcols == 0)
            break;
    }
    nprows = nprocs / npcols;

    int *usermap = new int[nprocs];
    for (int i = 0; i < nprows; i++)
    {
        for (int j = 0; j < npcols; j++)
            usermap[j * nprows + i] = i * npcols + j;
    }

    int what = 0, ic = 0;
    Cblacs_get(ic, what, &icontxt);
    Cblacs_pinfo(&mypnum, &nprocs);
    Cblacs_gridmap(&icontxt, usermap, nprows, nprows, npcols);
    Cblacs_gridinfo(icontxt, &nprows, &npcols, &myprow, &mypcol);
    delete[] usermap;
}

inline int na_rc(int n, int nblk, int nprows, int myprow) // calculate the row/column number in each process
{
    int i, na = 0;
    for (i = 0; i < n; i++)
        if ((i % (nblk * nprows)) / nblk == myprow)
            na++;
    return na;
}

// distribute hmatrix to each process
// hmatrix is row-first, and new_matrix is column-first
template <class T>
void distribute_data(T *hmatrix, T *new_matrix, int &nFull, int &nblk, int &na_rows, int &na_cols, int &icontxt)
{
    int mypnum;
    MPI_Comm_rank(MPI_COMM_WORLD, &mypnum);
    int nprows, npcols, myprow, mypcol;
    Cblacs_gridinfo(icontxt, &nprows, &npcols, &myprow, &mypcol);

    int SENDPROW = 0, SENDPCOL = 0, tag = 0;
    T *tmp = new T[nblk];

    // do iteration for matrix, distribute old_matrix to each process, pass a block each time
    for (int row = 0; row < nFull; row++)
    {
        int recv_prow = (row / nblk) % nprows; // the row number of recive process
        int nm_row = ((row / nblk) / nprows) * nblk + row % nblk; // row number of block in new_matrix
        for (int col = 0; col < nFull; col += nblk)
        {
            int recv_pcol = (col / nblk) % npcols; // the column number of recive process
            int nm_col = ((col / nblk) / npcols) * nblk + col % nblk;
            int pass_length = std::min(nFull - col, nblk); // nFull may not be devided by nblk;

            if (myprow == SENDPROW && mypcol == SENDPCOL)
            {
                if (myprow == recv_prow && mypcol == recv_pcol)
                {
                    for (int i = 0; i < pass_length; i++)
                    {
                        new_matrix[(nm_col + i) * na_rows + nm_row] = hmatrix[nFull * row + col + i];
                    }
                }
                else
                {
                    if (std::is_same<T, double>::value)
                        MPI_Send(&(hmatrix[nFull * row + col]),
                                 pass_length,
                                 MPI_DOUBLE,
                                 recv_prow * npcols + recv_pcol,
                                 tag,
                                 MPI_COMM_WORLD);
                    else
                        MPI_Send(&(hmatrix[nFull * row + col]),
                                 pass_length,
                                 MPI_DOUBLE_COMPLEX,
                                 recv_prow * npcols + recv_pcol,
                                 tag,
                                 MPI_COMM_WORLD);
                }
            }
            else if (myprow == recv_prow && mypcol == recv_pcol)
            {
                if (std::is_same<T, double>::value)
                    MPI_Recv(tmp, pass_length, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                else
                    MPI_Recv(tmp, pass_length, MPI_DOUBLE_COMPLEX, 0, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                for (int i = 0; i < pass_length; i++)
                {
                    new_matrix[(nm_col + i) * na_rows + nm_row] = tmp[i];
                }
            }
        }
    }
    delete[] tmp;
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

} // namespace LCAO_DIAGO_TEST
