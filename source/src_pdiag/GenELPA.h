//------------------------------------>8======================================
//  Copyright (c) 2016, Yu Shen (shenyu@ustc.edu.cn)
//  All rights reserved.
//
//  Redistribution and use in source and binary forms, with or without
//  modification, are permitted provided that the following conditions are met:
//      * Redistributions of source code must retain the above copyright
//        notice, this list of conditions and the following disclaimer.
//      * Redistributions in binary form must reproduce the above copyright
//        notice, this list of conditions and the following disclaimer in the
//        documentation and/or other materials provided with the distribution.
//      * Neither the name of the <organization> nor the
//        names of its contributors may be used to endorse or promote products
//        derived from this software without specific prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//  DISCLAIMED. IN NO EVENT Yu Shen SHALL BE LIABLE FOR ANY
//  DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
//  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
//  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
//  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
//  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//====================================8<----------------------------------------
#ifdef __MPI
#include <complex>
using namespace std;
const int  CHOLESKY_CRITICAL_SIZE=1000;
// Sometime ELPA fails on small system while pdotrf from scalapack is lower efficient on big system,
// therefore we combined them together, and this value is used to decide whether we use elpa or scalapack

int pdSolveGenEigen1(int nev, int nFull, int narows, int nacols, int *desc,
                     double *a, double *b, double *ev, double *q, double *work,
                     MPI_Comm mpi_comm_world, int my_blacs_ctxt,
                     int &method,
                     bool wantEigenVector, bool wantDebug);

int pdSolveGenEigen2(int nev, int nFull, int narows, int nacols, int *desc,
                     double *a, double *b, double *ev, double *q, double *work,
                     MPI_Comm mpi_comm_world, int my_blacs_ctxt,
                     int &method, int THIS_REAL_ELPA_KERNEL_API, int useQR,
                     bool wantEigenVector, bool wantDebug);

int pdDecomposeRightMatrix1(int nFull, int narows, int nacols, int *desc,
                           double *b, double *ev, double *q, double *work,
                           MPI_Comm mpi_comm_world, int mpi_comm_rows, int mpi_comm_cols,
                           int &method);


int pdDecomposeRightMatrix2(int nFull, int narows, int nacols, int *desc,
                           double *b, double *ev, double *q, double *work,
                           MPI_Comm mpi_comm_world, int mpi_comm_rows, int mpi_comm_cols,
                           int &method, int THIS_REAL_ELPA_KERNEL_API, int useQR);

int pdCheloskyDecomposeRightMatrix(int nFull, int narows, int nacols, int *desc, double *b,
                                   MPI_Comm mpi_comm_world, int mpi_comm_rows, int mpi_comm_cols,
                                   int &method,
                                   bool wantDebug);

int pdDiagonalizeRightMatrix1(int nFull, int narows, int nacols, int *desc,
                              double *b, double *ev, double *q, double *work,
                              int mpi_comm_rows, int mpi_comm_cols);

int pdDiagonalizeRightMatrix2(int nFull, int narows, int nacols, int *desc,
                              double *b, double *ev, double *q, double *work,
                              MPI_Comm mpi_comm_world, int mpi_comm_rows, int mpi_comm_cols,
                              int THIS_REAL_ELPA_KERNEL_API, int useQR);

int pdSolveEigen1(int nev, int nFull, int narows, int nacols, int *desc,
                 double *a, double *b, double *ev, double *q, double *work,
                 int mpi_comm_rows, int mpi_comm_cols,
                 int method,
                 bool wantEigenVector, bool wantDebug);

int pdSolveEigen2(int nev, int nFull, int narows, int nacols, int *desc,
                 double *a, double *b, double *ev, double *q, double *work,
                 MPI_Comm mpi_comm_world, int mpi_comm_rows, int mpi_comm_cols,
                 int method, int THIS_REAL_ELPA_KERNEL_API, int useQR,
                 bool wantEigenVector, bool wantDebug);

int pzSolveGenEigen1(int nev, int nFull, int narows, int nacols, int *desc,
                     std::complex<double> *a, std::complex<double> *b, double *ev, std::complex<double> *q, std::complex<double> *work,
                     MPI_Comm mpi_comm_world, int my_blacs_ctxt,
                     int &method,
                     bool wantEigenVector, bool wantDebug);

int pzSolveGenEigen2(int nev, int nFull, int narows, int nacols, int *desc,
                     std::complex<double> *a, std::complex<double> *b, double *ev, std::complex<double> *q, std::complex<double> *work,
                     MPI_Comm mpi_comm_world, int my_blacs_ctxt,
                     int &method, int THIS_REAL_ELPA_KERNEL_API,
                     bool wantEigenVector, bool wantDebug);

int pzDecomposeRightMatrix1(int nFull, int narows, int nacols, int *desc,
                           std::complex<double> *b, double *ev, std::complex<double> *q, std::complex<double> *work,
                           MPI_Comm mpi_comm_world, int mpi_comm_rows, int mpi_comm_cols,
                           int &method);

int pzDecomposeRightMatrix2(int nFull, int narows, int nacols, int *desc,
                           std::complex<double> *b, double *ev, std::complex<double> *q, std::complex<double> *work,
                           MPI_Comm mpi_comm_world, int mpi_comm_rows, int mpi_comm_cols,
                           int &method, int THIS_REAL_ELPA_KERNEL_API);

int pzCheloskyDecomposeRightMatrix(int nFull, int narows, int nacols, int *desc, std::complex<double> *b,
                                   MPI_Comm mpi_comm_world, int mpi_comm_rows, int mpi_comm_cols,
                                   int &method,
                                   bool wantDebug);

int pzDiagonalizeRightMatrix1(int nFull, int narows, int nacols, int *desc,
                              std::complex<double> *b, double *ev, std::complex<double> *q, std::complex<double> *work,
                              int mpi_comm_rows, int mpi_comm_cols);

int pzDiagonalizeRightMatrix2(int nFull, int narows, int nacols, int *desc,
                              std::complex<double> *b, double *ev, std::complex<double> *q, std::complex<double> *work,
                              MPI_Comm mpi_comm_world, int mpi_comm_rows, int mpi_comm_cols,
                              int THIS_REAL_ELPA_KERNEL_API);

int pzSolveEigen1(int nev, int nFull, int narows, int nacols, int *desc,
                 std::complex<double> *a, std::complex<double> *b, double *ev, std::complex<double> *q, std::complex<double> *work,
                 int mpi_comm_rows, int mpi_comm_cols,
                 int method,
                 bool wantEigenVector, bool wantDebug);

int pzSolveEigen2(int nev, int nFull, int narows, int nacols, int *desc,
                 std::complex<double> *a, std::complex<double> *b, double *ev, std::complex<double> *q, std::complex<double> *work,
                 MPI_Comm mpi_comm_world, int mpi_comm_rows, int mpi_comm_cols,
                 int method, int THIS_REAL_ELPA_KERNEL_API,
                 bool wantEigenVector, bool wantDebug);

#endif
