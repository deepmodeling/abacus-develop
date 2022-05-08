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
extern "C"
{
#include "../my_elpa.h"
}

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

int elpa_sethandle(elpa_t &handle,
                   int nFull,
                   int nev,
                   int local_nrows,
                   int local_ncols,
                   int nblk,
                   int process_row,
                   int process_col,
                   int icontxt)
{
    int error;
    elpa_init(20210430);
    handle = elpa_allocate(&error);
    if (error != 0)
        printf("ERROR: elpa_allocate error=%d!\n", error);

    elpa_set_integer(handle, "na", nFull, &error);
    if (error != 0)
        printf("ERROR: set na error=%d!\n", error);

    elpa_set_integer(handle, "nev", nev, &error);
    if (error != 0)
        printf("ERROR: set nev error=%d!\n", error);

    elpa_set_integer(handle, "local_nrows", local_nrows, &error);
    if (error != 0)
        printf("ERROR: set local_nrows error=%d!\n", error);

    elpa_set_integer(handle, "local_ncols", local_ncols, &error);
    if (error != 0)
        printf("ERROR: set local_ncols error=%d!\n", error);

    elpa_set_integer(handle, "nblk", nblk, &error);
    if (error != 0)
        printf("ERROR: set nblk error=%d!\n", error);

    elpa_set_integer(handle, "mpi_comm_parent", MPI_Comm_c2f(MPI_COMM_WORLD), &error);
    if (error != 0)
        printf("ERROR: set mpi_comm_parent error=%d!\n", error);

    elpa_set_integer(handle, "process_row", process_row, &error);
    if (error != 0)
        printf("ERROR: set process_row error=%d!\n", error);

    elpa_set_integer(handle, "process_col", process_col, &error);
    if (error != 0)
        printf("ERROR: set process_col error=%d!\n", error);

    elpa_set_integer(handle, "blacs_context", (int)icontxt, &error);
    if (error != 0)
        printf("ERROR: set blacs_context error=%d!\n", error);

    elpa_set_integer(handle, "cannon_for_generalized", 0, &error);
    if (error != 0)
        printf("ERROR: set cannon_for_generalized error=%d!\n", error);

    /* Setup */
    elpa_setup(handle); /* Set tunables */
    return 0;
}

void elpa_ev_solver(elpa_t &handle,
                    double *new_hmatrix,
                    double *new_smatrix,
                    double *e_elpa,
                    double *q_elpa,
                    int &is_already_decomposed,
                    int &error)
{
    // elpa_generalized_eigenvectors_d(handle, new_hmatrix, new_smatrix, e_elpa, q_elpa, is_already_decomposed, &error);
    elpa_generalized_eigenvalues_d(handle, new_hmatrix, new_smatrix, e_elpa, is_already_decomposed, &error);
}

void elpa_ev_solver(elpa_t &handle,
                    std::complex<double> *new_hmatrix,
                    std::complex<double> *new_smatrix,
                    double *e_elpa,
                    std::complex<double> *q_elpa,
                    int &is_already_decomposed,
                    int &error)
{
    elpa_generalized_eigenvectors_dc(handle,
                                     reinterpret_cast<double _Complex *>(new_hmatrix),
                                     reinterpret_cast<double _Complex *>(new_smatrix),
                                     e_elpa,
                                     reinterpret_cast<double _Complex *>(q_elpa),
                                     is_already_decomposed,
                                     &error);
}

template <class T>
void elpa_diago(T *hmatrix,
                T *smatrix,
                double *e_elpa,
                int &nFull,
                int &nblk,
                bool print_new_hmatrix = false,
                bool detail_info = false) // hmatrix is a 2D matrix, row-first.
{
    int nprocs, mypnum;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &mypnum);

    // map the processes to 2D matrix
    int nprows, npcols, myprow, mypcol, icontxt;
    process_2d(nprows, npcols, myprow, mypcol, icontxt);
    if (mypnum == 0 && detail_info)
    {
        printf("nprocs=%d, nprows=%d, npcols=%d\n", nprocs, nprows, npcols);
    }

    int na_rows = na_rc(nFull, nblk, nprows, myprow); // the number of row of the new_matrix in each process
    int na_cols = na_rc(nFull, nblk, npcols, mypcol); // the number of column of the new_matrix in each process

    int desc[9];
    int ISRC = 0;
    int info;
    descinit_(desc, &nFull, &nFull, &nblk, &nblk, &ISRC, &ISRC, &icontxt, &na_rows, &info);
    if (info != 0)
    {
        printf("Invalid blacs-distribution. Abort!\n");
        exit(1);
    }

    T *new_hmatrix = new T[na_rows * na_cols]; // this should be column-first
    T *new_smatrix = new T[na_rows * na_cols]; // this should be column-first
    distribute_data(hmatrix, new_hmatrix, nFull, nblk, na_rows, na_cols, icontxt);
    distribute_data(smatrix, new_smatrix, nFull, nblk, na_rows, na_cols, icontxt);

    MPI_Barrier(MPI_COMM_WORLD);

    if (print_new_hmatrix)
    {
        std::string file("./file_");
        file += std::to_string(mypnum) + std::string(".txt");
        std::ofstream fp(file);
        print_matrix(fp, new_hmatrix, na_rows, na_cols, false);
        fp.close();
    }

    elpa_t handle;
    int error;
    elpa_sethandle(handle, nFull, nFull, na_rows, na_cols, nblk, myprow, mypcol, icontxt);
    // setup solver
    // ELPA_SOLVER_2STAGE is not work when complied by GNU,
    // elpa_set_integer(handle, "solver", ELPA_SOLVER_2STAGE, &error);
    // elpa_set_integer(handle, "real_kernel", ELPA_2STAGE_REAL_AVX_BLOCK2, &error);

    double starttime, endtime;
    int is_already_decomposed = 0;
    T *q_elpa = new T[na_rows * na_cols];

    starttime = MPI_Wtime();
    elpa_ev_solver(handle, new_hmatrix, new_smatrix, e_elpa, q_elpa, is_already_decomposed, error);
    endtime = MPI_Wtime();

    if (mypnum == 0 && error != 0)
    {
        printf("ERROR: elpa_generalized_eigenvectors_d error=%d", error);
    }
    if (mypnum == 0 && detail_info)
    {
        printf("ELPA solving time is: %1f s\n", endtime - starttime);
    }

    elpa_deallocate(handle, &error);
    elpa_uninit(&error);

    delete[] new_hmatrix;
    delete[] new_smatrix;
    delete[] q_elpa;
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

double conj(double &a)
{
    return a;
}

// produce the random H/S matrix
template <class T> void random_hs(T *hmatrix, T *smatrix, int &n, int &sparsity)
{
    int min = 0;
    int max = 9999;
    std::default_random_engine e(100);
    std::uniform_int_distribution<unsigned> u(min, max);
    if (sparsity < 0)
        sparsity = 0;
    if (sparsity > 10)
        sparsity = 10;

    T *psi = new T[n * n];
    double *norm = new double[n];

    for (int i = 0; i < n; i++)
    {
        for (int j = i; j < n; j++)
        {
            double mincoef = 0.0;
            if (int(u(e) % 10) > int(sparsity - 1))
                mincoef = 1.0;
            double realp = pow(-1.0,u(e)%2) * static_cast<double>(u(e))/max;

            hmatrix[i * n + j] = realp;
            if (std::is_same<T, std::complex<double>>::value)
            {
                reinterpret_cast<double *>(&hmatrix[i * n + j])[1] = 0.0;
            }
            if (i != j)
            {
                hmatrix[j * n + i] = hmatrix[i * n + j] = mincoef * hmatrix[i * n + j];
            }
        }

        norm[i] = 0.0;
        for (int j = 0; j < n; j++)
        {
            psi[i * n + j] = pow(-1.0, u(e) % 2) * static_cast<double>(u(e)) / max;
            if (std::is_same<T, std::complex<double>>::value)
            {
                reinterpret_cast<double *>(&psi[i * n + j])[1] = pow(-1.0, u(e) % 2) * static_cast<double>(u(e)) / max;
            }
            norm[i] += std::norm(psi[i * n + j]);
        }
        norm[i] = sqrt(norm[i]);
    }

    for (int i = 0; i < n; i++)
    {
        for (int j = i; j < n; j++)
        {
            smatrix[i * n + j] = 0.0;
            if (std::is_same<T, std::complex<double>>::value)
            {
                reinterpret_cast<double *>(&hmatrix[i * n + j])[1] = 0.0;
            }

            for (int k = 0; k < n; k++)
                smatrix[i * n + j] += conj(psi[i * n + k]) * psi[j * n + k];
            smatrix[i * n + j] = smatrix[i * n + j] / norm[i] / norm[j];
            smatrix[j * n + i] = conj(smatrix[i * n + j]);
        }
    }
    delete[] psi;
    delete[] norm;
}

} // namespace LCAO_DIAGO_TEST
