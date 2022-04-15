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

//------------------------------------>8======================================
//  DEBUG step information
//  Macro Steps
//  Step        Function                        Description
//  A           p?DecomposeRightMatrix[1|2]     decomposing matrix B in the right side
//  A1          p?CheloskyDecomposeRightMatrix  using method 1 to decompose B
//  A2          p?CheloskyDecomposeRightMatrix  using method 2 to decompose B
//  A3          p?DiagonalizeRightMatrix[1|2]   using method 3 to decompose B
//  B           p?SolveEigen[1|2]               solving the standard eigenvalue equation
//  Micro steps:
//  Method  Step    Function                                    Description
//  1   1       p?potrf                                     cholesky decomposing B to U
//  1   2   for-cycle                                   clear_lower_part
//  1       3       elpa_invert_trm_[real|std::complex]              invert U to U^-1
//  1       4       p?symm                                      calculating A*U-1
//  1       5       p?gemm                                      calculating U-^*A*U-1
//  1       6       elpa_solve_evp_[real|std::complex]_[1|2]stage    solve eigenvalue
//  1       7       p?trmm                                      calculating eigenvector
//
//  2       1       elpa_cholesky_[real|std::complex]                cholesky decomposing B to U
//  2   2     for-cycle                     clear_lower_part
//  2       3       elpa_invert_trm_[real|std::complex]              invert U to U^-1
//  2       4       p?symm                                      calculating A*U-1
//  2       5       p?gemm                                      calculating U-^*A*U-1
//  2       6       elpa_solve_evp_[real|std::complex]_[1|2]stage    solve eigenvalue
//  2       7       p?trmm                                      calculating eigenvector
//
//  3       1       elpa_solve_evp_[real|std::complex]_[1|2]stage    calculating eigenvalue and eigenvector of B
//  3   2       for-cycle                       calculating q*ev
//  3       3       p?gemm                                      calculating B^-1/2
//  3       4       p?gemm                                      calculating B*A^T
//  3       5       p?gemm                                      calculating B*(B*A^T)^T
//  3       6       elpa_solve_evp_[real|std::complex]_[1|2]stage    solve eigenvalue
//  3       7       p?gemm                                      calculating eigenvector
// ?: d for double precision, z for double std::complex precision
//====================================8<----------------------------------------
#include <complex>
#include <cmath>
#include <cfloat>
#ifdef __MPI
#include "mpi.h"
#endif
extern "C"
{
//    #include "pblas.h"
    #include "../module_base/blacs_connector.h"
//    #include "scalapack.h"
    #include "my_elpa.h"
	#include "../module_base/scalapack_connector.h"
}
#include "GenELPA.h"

#ifdef _DEBUG
    #include <iostream>
    #include <sstream>
    #include <ctime>
    #include <cstring>
#endif

using namespace std;

#ifdef _DEBUG
    inline void timer(int myid, const char function[], int method, const char step[], double &t0)
    {
        //MPI_Barrier(MPI_COMM_WORLD);
        clock_t cur;
        cur=clock();
        double t1;
        std::stringstream outlog;
        if(t0<=0)  // t0 < 0 means this is the init call before the function
        {
            t0=double(cur);
            outlog.str("");
            outlog<<"DEBUG: Process "<<myid<<" Call "<<function<<std::endl;
            std::cout<<outlog.str();
        }
        else {
            t1=double(cur);
            outlog.str("");
            outlog<<"DEBUG: Process "<<myid<<" Method "<<method<<" Step "
                  <<step<<" "<<function<<" time: "<<(t1-t0)/CLOCKS_PER_SEC<<" s"<<std::endl;
            std::cout<<outlog.str();
        }
    }
#endif

int pdSolveGenEigen1(int nev, int nFull, int narows, int nacols, int *desc,
                     double *a, double *b, double *ev, double *q, double *work,
                     MPI_Comm mpi_comm_world, int my_blacs_ctxt,
                     int &method,
                     bool wantEigenVector, bool wantDebug)
{
    int info;
    int mpi_comm_rows, mpi_comm_cols;
    int nprows, npcols, myprow, mypcol;	// nblk;

    #ifdef _DEBUG
        int myid;
        time_t t0, t1;
        double t;
        std::stringstream outlog;
        MPI_Comm_rank(MPI_COMM_WORLD, &myid);
        outlog.str("");
        outlog<<"DEBUG: Process "<<myid<<" Enter pdSolveGenEigen1"<<std::endl;
        MPI_Barrier(MPI_COMM_WORLD);
        std::cout<<outlog.str();
    #endif

    info=0;
    Cblacs_gridinfo(my_blacs_ctxt, &nprows, &npcols, &myprow, &mypcol);
    info=get_elpa_communicators(mpi_comm_world, myprow, mypcol, &mpi_comm_rows, &mpi_comm_cols);

    if(info != 0) return info;

    #ifdef _DEBUG
        t=-1;
        timer(myid, "pdDecomposeRightMatrix1", method, "A", t);
    #endif
    info=pdDecomposeRightMatrix1(nFull, narows, nacols, desc,
                                b, ev, q, work,
                                mpi_comm_world, mpi_comm_rows, mpi_comm_cols, method);
    #ifdef _DEBUG
        timer(myid, "pdDecomposeRightMatrix1", method, "A", t);
    #endif

    if(info != 0) return info;

    #ifdef _DEBUG
        t=-1;
        timer(myid, "pdSolveEigen1", method, "B", t);
    #endif
    info=pdSolveEigen1(nev, nFull, narows, nacols, desc,
                       a, b, ev, q, work,
                       mpi_comm_rows, mpi_comm_cols, method,
                       wantEigenVector, wantDebug);
    #ifdef _DEBUG
        timer(myid, "pdSolveEigen1", method, "B", t);
    #endif
    return info;
}

int pdSolveGenEigen2(int nev, int nFull, int narows, int nacols, int *desc,
                     double *a, double *b, double *ev, double *q, double *work,
                     MPI_Comm mpi_comm_world, int my_blacs_ctxt,
                     int &method, int THIS_REAL_ELPA_KERNEL_API, int useQR,
                     bool wantEigenVector, bool wantDebug)
{
    int info;
    int mpi_comm_rows, mpi_comm_cols;
    int nprows, npcols, myprow, mypcol;	// nblk;

    #ifdef _DEBUG
        int myid;
        double t;
        std::stringstream outlog;
        MPI_Comm_rank(MPI_COMM_WORLD, &myid);
        outlog.str("");
        outlog<<"DEBUG: Process "<<myid<<" Enter pdSolveGenEigen2"<<std::endl;
        MPI_Barrier(MPI_COMM_WORLD);
        std::cout<<outlog.str();
    #endif

    info=0;
    Cblacs_gridinfo(my_blacs_ctxt, &nprows, &npcols, &myprow, &mypcol);
    info=get_elpa_communicators(mpi_comm_world, myprow, mypcol, &mpi_comm_rows, &mpi_comm_cols);

    if(info != 0) return info;
    #ifdef _DEBUG
        t=-1;
        timer(myid, "pdDecomposeRightMatrix2", method, "A", t);
    #endif
    info=pdDecomposeRightMatrix2(nFull, narows, nacols, desc,
                                b, ev, q, work,
                                mpi_comm_world, mpi_comm_rows, mpi_comm_cols, method,
                                THIS_REAL_ELPA_KERNEL_API, useQR);
    #ifdef _DEBUG
        timer(myid, "pdDecomposeRightMatrix2", method, "A", t);
    #endif

    if(info != 0) return info;

    #ifdef _DEBUG
        t=-1;
        timer(myid, "pdSolveEigen2", method, "B", t);
    #endif
    info=pdSolveEigen2(nev, nFull, narows, nacols, desc,
                       a, b, ev, q, work,
                       mpi_comm_world, mpi_comm_rows, mpi_comm_cols, method,
                       THIS_REAL_ELPA_KERNEL_API, useQR,
                       wantEigenVector, wantDebug);
    #ifdef _DEBUG
        timer(myid, "pdSolveEigen2", method, "B", t);
    #endif
    return info;
}

int pdDecomposeRightMatrix1(int nFull, int narows, int nacols, int *desc,
                           double *b, double *ev, double *q, double *work,
                           MPI_Comm mpi_comm_world, int mpi_comm_rows, int mpi_comm_cols,
                           int &method)
{
    int info=0; // for elpa functions, 1 is for success, 0 is for failure
    int wantDebug=false;
    #ifdef _DEBUG
        int myid;
        double t=-1;
        std::stringstream outlog;
        MPI_Comm_rank(MPI_COMM_WORLD, &myid);
        outlog.str("");
        outlog<<"DEBUG: Process "<<myid<<" Enter pdDecomposeRightMatrix1"<<std::endl;
        MPI_Barrier(MPI_COMM_WORLD);
        std::cout<<outlog.str();
    #endif
    switch(method)
    {
        case 0:
            #ifdef _DEBUG
                t=-1;
                timer(myid, "pdCheloskyDecomposeRightMatrix", method, "A", t);
            #endif
            info=pdCheloskyDecomposeRightMatrix(nFull, narows, nacols, desc, b,
                                            mpi_comm_world, mpi_comm_rows, mpi_comm_cols, method, wantDebug);
            #ifdef _DEBUG
                timer(myid, "pdCheloskyDecomposeRightMatrix", method, "A", t);
            #endif
            if(info!=0)
            {
                method=3;
                #ifdef _DEBUG
                    t=-1;
                    timer(myid, "pdDiagonalizeRightMatrix1", method, "A3", t);
                #endif
                info=pdDiagonalizeRightMatrix1(nFull, narows, nacols, desc,
                                            b, ev, q, work,
                                            mpi_comm_rows, mpi_comm_cols);
                #ifdef _DEBUG
                    timer(myid, "pdDiagonalizeRightMatrix1", method, "A3", t);
                #endif
            }
            break;
        case 1:
        case 2:
            #ifdef _DEBUG
                t=-1;
                timer(myid, "pdCheloskyDecomposeRightMatrix", method, "A", t);
            #endif
            info=pdCheloskyDecomposeRightMatrix(nFull, narows, nacols, desc, b,
                                            mpi_comm_world, mpi_comm_rows, mpi_comm_cols, method, wantDebug);
            #ifdef _DEBUG
                timer(myid, "pdCheloskyDecomposeRightMatrix", method, "A", t);
            #endif
            break;
        case 3:
            #ifdef _DEBUG
                t=-1;
                timer(myid, "pdDiagonalizeRightMatrix1", method, "A3", t);
            #endif
            info=pdDiagonalizeRightMatrix1(nFull, narows, nacols, desc,
                                        b, ev, q, work,
                                        mpi_comm_rows, mpi_comm_cols);
            #ifdef _DEBUG
                timer(myid, "pdDiagonalizeRightMatrix1", method, "A3", t);
            #endif
            break;
        default:
            info=10;
    }
    return info;
}

int pdDecomposeRightMatrix2(int nFull, int narows, int nacols, int *desc,
                           double *b, double *ev, double *q, double *work,
                           MPI_Comm mpi_comm_world, int mpi_comm_rows, int mpi_comm_cols,
                           int &method, int THIS_REAL_ELPA_KERNEL_API, int useQR)
{
    int info=0; // for elpa functions, 1 is for success, 0 is for failure
    int wantDebug=false;
    //int allinfo=0;

    #ifdef _DEBUG
        int myid;
        double t;
        std::stringstream outlog;
        MPI_Comm_rank(MPI_COMM_WORLD, &myid);
        outlog.str("");
        outlog<<"DEBUG: Process "<<myid<<" Enter pdDecomposeRightMatrix2"<<std::endl;
        MPI_Barrier(MPI_COMM_WORLD);
        std::cout<<outlog.str();
    #endif
    switch(method)
    {
        case 0:
            #ifdef _DEBUG
                t=-1;
                timer(myid, "pdCheloskyDecomposeRightMatrix", method, "A", t);
            #endif
            info=pdCheloskyDecomposeRightMatrix(nFull, narows, nacols, desc, b,
                                            mpi_comm_world, mpi_comm_rows, mpi_comm_cols, method, wantDebug);
            #ifdef _DEBUG
                timer(myid, "pdCheloskyDecomposeRightMatrix", method, "A", t);
            #endif
            if(info!=0)
            {
                method=3;
                #ifdef _DEBUG
                    t=-1;
                    timer(myid, "pdDiagonalizeRightMatrix2", method, "A3", t);
                #endif
                info=pdDiagonalizeRightMatrix2(nFull, narows, nacols, desc,
                                            b, ev, q, work,
                                            mpi_comm_world, mpi_comm_rows, mpi_comm_cols,
                                            THIS_REAL_ELPA_KERNEL_API, useQR);
                #ifdef _DEBUG
                    timer(myid, "pdDiagonalizeRightMatrix2", method, "A3", t);
                #endif
            }
            break;
        case 1:
        case 2:
            #ifdef _DEBUG
                t=-1;
                timer(myid, "pdCheloskyDecomposeRightMatrix", method, "A", t);
            #endif
            info=pdCheloskyDecomposeRightMatrix(nFull, narows, nacols, desc, b,
                                            mpi_comm_world, mpi_comm_rows, mpi_comm_cols, method, wantDebug);
            #ifdef _DEBUG
                timer(myid, "pdCheloskyDecomposeRightMatrix", method, "A", t);
            #endif
            break;
        case 3:
            #ifdef _DEBUG
                t=-1;
                timer(myid, "pdDiagonalizeRightMatrix2", method, "A3", t);
            #endif
            info=pdDiagonalizeRightMatrix2(nFull, narows, nacols, desc,
                                        b, ev, q, work,
                                        mpi_comm_world, mpi_comm_rows, mpi_comm_cols,
                                        THIS_REAL_ELPA_KERNEL_API, useQR);
            #ifdef _DEBUG
                timer(myid, "pdDiagonalizeRightMatrix2", method, "A3", t);
            #endif
            break;
        default:
            info=10;
    }
    return info;
}

int pdCheloskyDecomposeRightMatrix(int nFull, int narows, int nacols, int *desc, double *b,
                                   MPI_Comm mpi_comm_world, int mpi_comm_rows, int mpi_comm_cols,
                                   int &method,
                                   bool wantDebug)
{
    int info, allinfo;
    int my_blacs_ctxt;
    int nprows, npcols, myprow, mypcol, nblk;
    char uplo;	//transa, transb, side, diag;
    int isrc=1, jsrc=1;
    //double alpha, beta;
    int lda;

    #ifdef _DEBUG
        int myid;
        double t;
        std::stringstream outlog;
        MPI_Comm_rank(MPI_COMM_WORLD, &myid);
        outlog.str("");
        outlog<<"DEBUG: Process "<<myid<<" Enter pdCheloskyDecomposeRightMatrix"<<std::endl;
        MPI_Barrier(MPI_COMM_WORLD);
        std::cout<<outlog.str();
    #endif
    //info=1; // for elpa functions, 1 is for success, 0 is for failure
    my_blacs_ctxt=desc[1];
    nblk=desc[4];
    lda=desc[8];
    Cblacs_gridinfo(my_blacs_ctxt, &nprows, &npcols, &myprow, &mypcol);

    // do cholesky decomposing of overlap matrix B: B=U^T*U, put U to b
    if(method==1 || (method==0 && nFull<CHOLESKY_CRITICAL_SIZE))
    {
        method=1;
        uplo='U';
        #ifdef _DEBUG
            t=-1;
            timer(myid, "pdpotrf_", method, "1", t);
        #endif
        pdpotrf_(&uplo, &nFull, b, &isrc, &jsrc, desc, &info);
        #ifdef _DEBUG
            timer(myid, "pdpotrf_", method, "1", t);
        #endif
        MPI_Allreduce(&info, &allinfo, 1, MPI_INT, MPI_SUM, mpi_comm_world);
        if(allinfo != 0)
        {
            if(method==1)
            {
                info=10;
                return info;
            }
            else
            {
                method=2;
                #ifdef _DEBUG
                    t=-1;
                    timer(myid, "elpa_cholesky_real", method, "1", t);
                #endif
                info=elpa_cholesky_real(nFull, b, narows, nblk, nacols, mpi_comm_rows, mpi_comm_cols, wantDebug);
                #ifdef _DEBUG
                    timer(myid, "elpa_cholesky_real", method, "1", t);
                #endif
                MPI_Allreduce(&info, &allinfo, 1, MPI_INT, MPI_MIN, mpi_comm_world);
                if(allinfo != 1)
                {
                    info=10;
                    return info;
                }
            }
        }
    }
    else if(method==2 || (method==0 && nFull>CHOLESKY_CRITICAL_SIZE))
    {
        method=2;
        #ifdef _DEBUG
            t=-1;
            timer(myid, "elpa_cholesky_real", method, "1", t);
        #endif
        wantDebug=true;
        info=elpa_cholesky_real(nFull, b, narows, nblk, nacols, mpi_comm_rows, mpi_comm_cols, wantDebug);
        #ifdef _DEBUG
            timer(myid, "elpa_cholesky_real", method, "1", t);
        #endif
        MPI_Allreduce(&info, &allinfo, 1, MPI_INT, MPI_MIN, mpi_comm_world);
        if(allinfo != 1)
        {
            if(method==2)
            {
                info=10;
                return info;
            }
            else
            {
                method=1;
                uplo='U';
                #ifdef _DEBUG
                    t=-1;
                    timer(myid, "pdpotrf_", method, "1", t);
                #endif
                pdpotrf_(&uplo, &nFull, b, &isrc, &jsrc, desc, &info);
                #ifdef _DEBUG
                    timer(myid, "pdpotrf_", method, "1", t);
                #endif
                MPI_Allreduce(&info, &allinfo, 1, MPI_INT, MPI_SUM, mpi_comm_world);
                if(allinfo != 0)
                {
                    info=10;
                    return info;
                }
            }
        }
    }
    else
    {
        info=20;
        return info;
    }

    //clear_lower_part of b
    #ifdef _DEBUG
        t=-1;
        timer(myid, "clear_lower_part", method, "2", t);
    #endif
    for(int j=0; j<nacols; ++j)
    {
        int jGlobal=Local_Orbital_wfc::globalIndex(j, nblk, npcols, mypcol);
        for(int i=0; i<narows; ++i)
        {
            int iGlobal=Local_Orbital_wfc::globalIndex(i, nblk, nprows, myprow);
            if(iGlobal>jGlobal) b[i+j*narows]=0;
        }
    }
    #ifdef _DEBUG
        timer(myid, "clear_lower_part", method, "2", t);
    #endif
    // invert U to U^-1, put to b
    #ifdef _DEBUG
        t=-1;
        timer(myid, "elpa_invert_trm_real", method, "3", t);
    #endif
    info=elpa_invert_trm_real(nFull, b, lda, nblk, nacols, mpi_comm_rows, mpi_comm_cols, wantDebug);
    #ifdef _DEBUG
        timer(myid, "elpa_invert_trm_real", method, "3", t);
    #endif
    if(info != 1)
    {
        info=20;
        return info;
    }

    return 0;
}


int pdDiagonalizeRightMatrix1(int nFull, int narows, int nacols, int *desc,
                              double *b, double *ev, double *q, double *work,
                              int mpi_comm_rows, int mpi_comm_cols)
{
    int info;
    int my_blacs_ctxt;
    int nprows, npcols, myprow, mypcol, nblk;
    char transa, transb;
    int isrc=1, jsrc=1;
    double alpha, beta;
    int lda;

    #ifdef _DEBUG
        int myid;
        double t;
        int method=3;
        std::stringstream outlog;
        MPI_Comm_rank(MPI_COMM_WORLD, &myid);
        outlog.str("");
        outlog<<"DEBUG: Process "<<myid<<" Enter pdDiagonalizeRightMatrix1"<<std::endl;
        MPI_Barrier(MPI_COMM_WORLD);
        std::cout<<outlog.str();
    #endif
    info=1; // for elpa functions, 1 is for success, 0 is for failure
    my_blacs_ctxt=desc[1];
    nblk=desc[4];
    lda=desc[8];
    Cblacs_gridinfo(my_blacs_ctxt, &nprows, &npcols, &myprow, &mypcol);

    // calculate the eigenvectors of overlap matrix b, put to q
    // calculate the eigenvalues of overlap matrix b, put to ev
    #ifdef _DEBUG
        t=-1;
        timer(myid, "elpa_solve_evp_real_1stage", method, "1", t);
    #endif
    info=elpa_solve_evp_real_1stage(nFull, nFull, b, lda, ev, q, lda, nblk,
                                    nacols, mpi_comm_rows, mpi_comm_cols);
    #ifdef _DEBUG
        timer(myid, "elpa_solve_evp_real_1stage", method, "1", t);
    #endif
    if(info != 1)
    {
        info=100;
        return info;
    }

    // calculate B^{-1/2}_{i,j}=\sum_k q_{i,k}*ev_k^{-1/2}*q_{j,k} and put to b, which will be b^-1/2
    // calculate q*ev and put to work
    #ifdef _DEBUG
        t=-1;
        timer(myid, "q*ev", method, "2", t);
    #endif
    for(int i=0; i<nacols; ++i)
    {
        int eidx=Local_Orbital_wfc::globalIndex(i, nblk, npcols, mypcol);
        double ev_sqrt=ev[eidx]>DBL_MIN?1.0/sqrt(ev[eidx]):0;
        for(int j=0; j<narows; ++j)
            work[i*lda+j]=q[i*lda+j]*ev_sqrt;
    }
    #ifdef _DEBUG
        timer(myid, "q*ev", method, "2", t);
    #endif
    // calculate qevq=qev*q^T, put to b, which is B^{-1/2}
    transa='N';
    transb='T';
    alpha=1.0;
    beta=0.0;
    #ifdef _DEBUG
        t=-1;
        timer(myid, "pdgemm_", method, "3", t);
    #endif
    pdgemm_(&transa, &transb, &nFull, &nFull, &nFull,
            &alpha, work, &isrc, &jsrc, desc,
                    q,    &isrc, &jsrc, desc,
            &beta,  b,    &isrc, &jsrc, desc);
    #ifdef _DEBUG
        timer(myid, "pdgemm_", method, "3", t);
    #endif
    return 0;
}

int pdDiagonalizeRightMatrix2(int nFull, int narows, int nacols, int *desc,
                              double *b, double *ev, double *q, double *work,
                              MPI_Comm mpi_comm_world, int mpi_comm_rows, int mpi_comm_cols,
                              int THIS_REAL_ELPA_KERNEL_API, int useQR)
{
    int info;
    int my_blacs_ctxt;
    int nprows, npcols, myprow, mypcol, nblk;
    char transa, transb;
    int isrc=1, jsrc=1;
    double alpha, beta;
    int lda;

    #ifdef _DEBUG
        int myid;
        double t;
        int method=3;
        std::stringstream outlog;
        MPI_Comm_rank(MPI_COMM_WORLD, &myid);
        outlog.str("");
        outlog<<"DEBUG: Process "<<myid<<" Enter pdDiagonalizeRightMatrix2"<<std::endl;
        MPI_Barrier(MPI_COMM_WORLD);
        std::cout<<outlog.str();
    #endif
    info=1; // for elpa functions, 1 is for success, 0 is for failure
    my_blacs_ctxt=desc[1];
    nblk=desc[4];
    lda=desc[8];
    Cblacs_gridinfo(my_blacs_ctxt, &nprows, &npcols, &myprow, &mypcol);

    // calculate the eigenvectors of overlap matrix b, put to q
    // calculate the eigenvalues of overlap matrix b, put to ev
    #ifdef _DEBUG
        t=-1;
        timer(myid, "elpa_solve_evp_real_2stage", method, "1", t);
    #endif
    info=elpa_solve_evp_real_2stage(nFull, nFull, b, lda, ev, q, lda, nblk, nacols,
                                    mpi_comm_rows, mpi_comm_cols, mpi_comm_world,
                                    THIS_REAL_ELPA_KERNEL_API, useQR);
    #ifdef _DEBUG
        timer(myid, "elpa_solve_evp_real_2stage", method, "1", t);
    #endif
    if(info != 1)
    {
        info=100;
        return info;
    }

    // calculate B^{-1/2}_{i,j}=\sum_k q_{i,k}*ev_k^{-1/2}*q_{j,k} and put to b, which will be b^-1/2
    // calculate q*ev and put to work
    #ifdef _DEBUG
        t=-1;
        timer(myid, "q*ev", method, "2", t);
    #endif
    for(int i=0; i<nacols; ++i)
    {
        int eidx=Local_Orbital_wfc::globalIndex(i, nblk, npcols, mypcol);
        //double ev_sqrt=1.0/sqrt(ev[eidx]);
        double ev_sqrt=ev[eidx]>DBL_MIN?1.0/sqrt(ev[eidx]):0;
        for(int j=0; j<narows; ++j)
            work[i*lda+j]=q[i*lda+j]*ev_sqrt;
    }
    #ifdef _DEBUG
        timer(myid, "q*ev", method, "2", t);
    #endif
    // calculate qevq=qev*q^T, put to b, which is B^{-1/2}
    transa='N';
    transb='T';
    alpha=1.0;
    beta=0.0;
    #ifdef _DEBUG
        t=-1;
        timer(myid, "pdgemm_", method, "3", t);
    #endif
    pdgemm_(&transa, &transb, &nFull, &nFull, &nFull,
            &alpha, work, &isrc, &jsrc, desc,
                    q,    &isrc, &jsrc, desc,
            &beta,  b,    &isrc, &jsrc, desc);
    #ifdef _DEBUG
        timer(myid, "pdgemm_", method, "3", t);
    #endif
    return 0;
}

int pdSolveEigen1(int nev, int nFull, int narows, int nacols, int *desc,
                 double *a, double *b, double *ev, double *q, double *work,
                 int mpi_comm_rows, int mpi_comm_cols, int method,
                 bool wantEigenVector, bool wantDebug)
{
    int info;
    int my_blacs_ctxt;
    int nprows, npcols, myprow, mypcol, nblk;
    char transa, transb, side, uplo, diag;
    int isrc=1, jsrc=1;
    double alpha, beta;
    int lda;

    #ifdef _DEBUG
        int myid;
        double t;
        std::stringstream outlog;
        MPI_Comm_rank(MPI_COMM_WORLD, &myid);
        outlog.str("");
        outlog<<"DEBUG: Process "<<myid<<" Enter pdSolveEigen1"<<std::endl;
        MPI_Barrier(MPI_COMM_WORLD);
        std::cout<<outlog.str();
    #endif

    info=1; // for elpa functions, 1 is for success, 0 is for failure
    my_blacs_ctxt=desc[1];
    nblk=desc[4];
    lda=desc[8];
    Cblacs_gridinfo(my_blacs_ctxt, &nprows, &npcols, &myprow, &mypcol);

    if(method==1 || method==2) // use cholesky decomposing result of overlap matrix (b)
    {
        // calculate A*U^-1, put to work
        transa='T';
        transb='N';
        alpha=1.0;
        beta=0.0;
        #ifdef _DEBUG
            t=-1;
            timer(myid, "pdgemm_", method, "4", t);
        #endif
        pdgemm_(&transa, &transb, &nFull, &nFull, &nFull,
                &alpha, a, &isrc, &isrc, desc,
                        b, &isrc, &isrc, desc,
                &beta,  work, &isrc, &isrc, desc);
        #ifdef _DEBUG
            timer(myid, "pdgemm_", method, "4", t);
        #endif
        // calculate U^-T^(A*U^-1), put to a
        transa='T';
        transb='N';
        alpha=1.0;
        beta=0.0;
        #ifdef _DEBUG
            t=-1;
            timer(myid, "pdgemm_", method, "5", t);
        #endif
        pdgemm_(&transa, &transb, &nFull, &nFull, &nFull,
                &alpha, b, &isrc, &isrc, desc,
                        work, &isrc, &isrc, desc,
                &beta,  a, &isrc, &isrc, desc);
        #ifdef _DEBUG
            timer(myid, "pdgemm_", method, "5", t);
        #endif
        // calculate the eigenvalues and eigenvectors, put to ev and q
        #ifdef _DEBUG
            t=-1;
            timer(myid, "elpa_solve_evp_real_1stage", method, "6", t);
        #endif
        info=elpa_solve_evp_real_1stage(nFull, nev, a, lda, ev, q, lda, nblk,
                                        nacols, mpi_comm_rows, mpi_comm_cols);
        #ifdef _DEBUG
            timer(myid, "elpa_solve_evp_real_1stage", method, "6", t);
        #endif
        if(info != 1)
        {
            info=100;
            return info;
        }
        if(wantEigenVector)
        {
            // transform the eigenvectors to original general equation, let U^-1*q, and put to q
            side='L';
            uplo='U';
            transa='N';
            diag='N';
            #ifdef _DEBUG
                t=-1;
                timer(myid, "pdtrmm_", method, "7", t);
            #endif
            pdtrmm_(&side, &uplo, &transa,  &diag, &nFull, &nFull,
                    &alpha, b, &isrc, &jsrc, desc,
                            q, &isrc, &jsrc, desc);
            #ifdef _DEBUG
                timer(myid, "pdtrmm_", method, "7", t);
            #endif
        }
    }
    else if(method==3)  //use diagonalized result of overlap matrix (b)
    {
        // calculate b*a^T and put to work
        transa='N';
        transb='T';
        alpha=1.0;
        beta=0.0;
        #ifdef _DEBUG
            t=-1;
            timer(myid, "pdgemm_", method, "4", t);
        #endif
        pdgemm_(&transa, &transb, &nFull, &nFull, &nFull,
                &alpha, b,       &isrc, &jsrc, desc,
                        a,    &isrc, &jsrc, desc,
                &beta,  work, &isrc, &jsrc, desc);
        #ifdef _DEBUG
            timer(myid, "pdgemm_", method, "4", t);
        #endif
        // calculate b*work^T and put to a -- origian A*x=v*B*x was transform to a*x'=v*x'
        #ifdef _DEBUG
            t=-1;
            timer(myid, "pdgemm_", method, "5", t);
        #endif
        pdgemm_(&transa, &transb, &nFull, &nFull, &nFull,
                &alpha, b,       &isrc, &jsrc, desc,
                        work, &isrc, &jsrc, desc,
                &beta,  a,       &isrc, &jsrc, desc);
        #ifdef _DEBUG
            timer(myid, "pdgemm_", method, "5", t);
        #endif
        // calculate eigenvectors of a*x'=v*x'
        #ifdef _DEBUG
            t=-1;
            timer(myid, "elpa_solve_evp_real_1stage", method, "6", t);
        #endif
        info=elpa_solve_evp_real_1stage(nFull, nFull, a, lda, ev, work, lda, nblk,
                                        nacols, mpi_comm_rows, mpi_comm_cols);
        #ifdef _DEBUG
            timer(myid, "elpa_solve_evp_real_1stage", method, "6", t);
        #endif
        if(info != 1)
        {
            info=200;
            return info;
        }
        if(wantEigenVector)
        {
            // transform the eigenvectors to original general equation, let b^T*q, and put to q
            transa='T';
            transb='N';
            #ifdef _DEBUG
                t=-1;
                timer(myid, "pdgemm_", method, "4", t);
            #endif
            pdgemm_(&transa, &transb, &nFull, &nFull, &nFull,
                    &alpha, b,       &isrc, &jsrc, desc,
                            work, &isrc, &jsrc, desc,
                    &beta,  q,       &isrc, &jsrc, desc);
            #ifdef _DEBUG
                timer(myid, "pdgemm_", method, "4", t);
            #endif
        }
    }
    else
    {
        info=30;
        return info;
    }

    return 0;
}

int pdSolveEigen2(int nev, int nFull, int narows, int nacols, int *desc,
                 double *a, double *b, double *ev, double *q, double *work,
                 MPI_Comm mpi_comm_world, int mpi_comm_rows, int mpi_comm_cols, int method,
                 int THIS_REAL_ELPA_KERNEL_API, int useQR,
                 bool wantEigenVector, bool wantDebug)
{
    int info;
    int my_blacs_ctxt;
    int nprows, npcols, myprow, mypcol, nblk;
    char transa, transb, side, uplo, diag;
    int isrc=1, jsrc=1;
    double alpha, beta;
    int lda;

    #ifdef _DEBUG
        int myid;
        double t;
        std::stringstream outlog;
        MPI_Comm_rank(MPI_COMM_WORLD, &myid);
        outlog.str("");
        outlog<<"DEBUG: Process "<<myid<<" Enter pdSolveEigen2"<<std::endl;
        MPI_Barrier(MPI_COMM_WORLD);
        std::cout<<outlog.str();
    #endif

    info=1; // for elpa functions, 1 is for success, 0 is for failure
    my_blacs_ctxt=desc[1];
    nblk=desc[4];
    lda=desc[8];
    Cblacs_gridinfo(my_blacs_ctxt, &nprows, &npcols, &myprow, &mypcol);

    if(method==1||method==2) // use chmelosky decomposing result of overlap matrix (b)
    {
        // calculate A*U^-1, put to work
        transa='T';
        transb='N';
        alpha=1.0;
        beta=0.0;
        #ifdef _DEBUG
            t=-1;
            timer(myid, "pdgemm_", method, "4", t);
        #endif
        pdgemm_(&transa, &transb, &nFull, &nFull, &nFull,
                &alpha, a, &isrc, &isrc, desc,
                        b, &isrc, &isrc, desc,
                &beta,  work, &isrc, &isrc, desc);
        #ifdef _DEBUG
            timer(myid, "pdgemm_", method, "4", t);
        #endif
        // calculate U^-T^(A*U^-1), put to a
        transa='T';
        transb='N';
        alpha=1.0;
        beta=0.0;
        #ifdef _DEBUG
            t=-1;
            timer(myid, "pdgemm_", method, "5", t);
        #endif
        pdgemm_(&transa, &transb, &nFull, &nFull, &nFull,
                &alpha, b, &isrc, &isrc, desc,
                        work, &isrc, &isrc, desc,
                &beta,  a, &isrc, &isrc, desc);
        #ifdef _DEBUG
            timer(myid, "pdgemm_", method, "5", t);
        #endif
        // calculate the eigenvalues and eigenvectors, put to ev and q
        #ifdef _DEBUG
            t=-1;
            timer(myid, "elpa_solve_evp_real_2stage", method, "6", t);
        #endif
        info=elpa_solve_evp_real_2stage(nFull, nev, a, lda, ev, q, lda, nblk, nacols,
                                        mpi_comm_rows, mpi_comm_cols, mpi_comm_world,
                                        THIS_REAL_ELPA_KERNEL_API, useQR);
        #ifdef _DEBUG
            timer(myid, "elpa_solve_evp_real_2stage", method, "6", t);
        #endif
        if(info != 1)
        {
            info=100;
            return info;
        }
        if(wantEigenVector)
        {
            // transform the eigenvectors to original general equation, let U^-1*q, and put to q
            side='L';
            uplo='U';
            transa='N';
            diag='N';
            #ifdef _DEBUG
                t=-1;
                timer(myid, "pdtrmm_", method, "7", t);
            #endif
            pdtrmm_(&side, &uplo, &transa,  &diag, &nFull, &nFull,
                    &alpha, b, &isrc, &jsrc, desc,
                            q, &isrc, &jsrc, desc);
            #ifdef _DEBUG
                timer(myid, "pdtrmm_", method, "7", t);
            #endif
        }
    }
    else if(method==3)  //use diagonalized result of overlap matrix (b)
    {
        // calculate b*a^T and put to work
        transa='N';
        transb='T';
        alpha=1.0;
        beta=0.0;
        #ifdef _DEBUG
            t=-1;
            timer(myid, "pdgemm_", method, "4", t);
        #endif
        pdgemm_(&transa, &transb, &nFull, &nFull, &nFull,
                &alpha, b,       &isrc, &jsrc, desc,
                        a,    &isrc, &jsrc, desc,
                &beta,  work, &isrc, &jsrc, desc);
        #ifdef _DEBUG
            timer(myid, "pdgemm_", method, "4", t);
        #endif
        // calculate b*work^T and put to a -- origian A*x=v*B*x was transform to a*x'=v*x'
        #ifdef _DEBUG
            t=-1;
            timer(myid, "pdgemm_", method, "5", t);
        #endif
        pdgemm_(&transa, &transb, &nFull, &nFull, &nFull,
                &alpha, b,       &isrc, &jsrc, desc,
                        work, &isrc, &jsrc, desc,
                &beta,  a,       &isrc, &jsrc, desc);
        #ifdef _DEBUG
            timer(myid, "pdgemm_", method, "5", t);
        #endif
        // calculate eigenvectors of a*x'=v*x'
        #ifdef _DEBUG
            t=-1;
            timer(myid, "elpa_solve_evp_real_2stage", method, "6", t);
        #endif
        info=elpa_solve_evp_real_2stage(nFull, nev, a, lda, ev, work, lda, nblk, nacols,
                                        mpi_comm_rows, mpi_comm_cols, mpi_comm_world,
                                        THIS_REAL_ELPA_KERNEL_API, useQR);
        #ifdef _DEBUG
            timer(myid, "elpa_solve_evp_real_2stage", method, "6", t);
        #endif
        if(info != 1)
        {
            info=200;
            return info;
        }
        if(wantEigenVector)
        {
            // transform the eigenvectors to original general equation, let b^T*q, and put to q
            transa='T';
            transb='N';
            #ifdef _DEBUG
                t=-1;
                timer(myid, "pdgemm_", method, "7", t);
            #endif
            pdgemm_(&transa, &transb, &nFull, &nFull, &nFull,
                    &alpha, b,       &isrc, &jsrc, desc,
                            work, &isrc, &jsrc, desc,
                    &beta,  q,       &isrc, &jsrc, desc);
            #ifdef _DEBUG
                timer(myid, "pdgemm_", method, "7", t);
            #endif
        }
    }
    else
    {
        info=30;
        return info;
    }
    return 0;
}


int pzSolveGenEigen1(int nev, int nFull, int narows, int nacols, int *desc,
                     std::complex<double> *a, std::complex<double> *b, double *ev, std::complex<double> *q, std::complex<double> *work,
                     MPI_Comm mpi_comm_world, int my_blacs_ctxt,
                     int &method,
                     bool wantEigenVector, bool wantDebug)
{
    int info;
    int mpi_comm_rows, mpi_comm_cols;
    int nprows, npcols, myprow, mypcol;	// nblk;

    #ifdef _DEBUG
        int myid;
        time_t t0, t1;
        double t;
        std::stringstream outlog;
        MPI_Comm_rank(MPI_COMM_WORLD, &myid);
        outlog.str("");
        outlog<<"DEBUG: Process "<<myid<<" Enter pzSolveGenEigen1"<<std::endl;
        MPI_Barrier(MPI_COMM_WORLD);
        std::cout<<outlog.str();
    #endif

    info=0;
    Cblacs_gridinfo(my_blacs_ctxt, &nprows, &npcols, &myprow, &mypcol);
    info=get_elpa_communicators(mpi_comm_world, myprow, mypcol, &mpi_comm_rows, &mpi_comm_cols);

    if(info != 0) return info;

    #ifdef _DEBUG
        t=-1;
        timer(myid, "pzDecomposeRightMatrix1", method, "A", t);
    #endif
    info=pzDecomposeRightMatrix1(nFull, narows, nacols, desc,
                                b, ev, q, work,
                                mpi_comm_world, mpi_comm_rows, mpi_comm_cols, method);
    #ifdef _DEBUG
        timer(myid, "pzDecomposeRightMatrix1", method, "A", t);
    #endif

    if(info != 0) return info;

    #ifdef _DEBUG
        t=-1;
        timer(myid, "pzSolveEigen1", method, "B", t);
    #endif
    info=pzSolveEigen1(nev, nFull, narows, nacols, desc,
                       a, b, ev, q, work,
                       mpi_comm_rows, mpi_comm_cols, method,
                       wantEigenVector, wantDebug);
    #ifdef _DEBUG
        timer(myid, "pzSolveEigen1", method, "B", t);
    #endif
    return info;
}

int pzSolveGenEigen2(int nev, int nFull, int narows, int nacols, int *desc,
                     std::complex<double> *a, std::complex<double> *b, double *ev, std::complex<double> *q, std::complex<double> *work,
                     MPI_Comm mpi_comm_world, int my_blacs_ctxt,
                     int &method, int THIS_REAL_ELPA_KERNEL_API,
                     bool wantEigenVector, bool wantDebug)
{
    int info;
    int mpi_comm_rows, mpi_comm_cols;
    int nprows, npcols, myprow, mypcol;	// nblk;

    #ifdef _DEBUG
        int myid;
        double t;
        std::stringstream outlog;
        MPI_Comm_rank(MPI_COMM_WORLD, &myid);
        outlog.str("");
        outlog<<"DEBUG: Process "<<myid<<" Enter pdSolveGenEigen2"<<std::endl;
        MPI_Barrier(MPI_COMM_WORLD);
        std::cout<<outlog.str();
    #endif

    info=0;
    Cblacs_gridinfo(my_blacs_ctxt, &nprows, &npcols, &myprow, &mypcol);
    info=get_elpa_communicators(mpi_comm_world, myprow, mypcol, &mpi_comm_rows, &mpi_comm_cols);

    if(info != 0) return info;
    #ifdef _DEBUG
        t=-1;
        timer(myid, "pzDecomposeRightMatrix2", method, "A", t);
    #endif
    info=pzDecomposeRightMatrix2(nFull, narows, nacols, desc,
                                b, ev, q, work,
                                mpi_comm_world, mpi_comm_rows, mpi_comm_cols, method,
                                THIS_REAL_ELPA_KERNEL_API);
    #ifdef _DEBUG
        timer(myid, "pzDecomposeRightMatrix2", method, "A", t);
    #endif

    if(info != 0) return info;

    #ifdef _DEBUG
        t=-1;
        timer(myid, "pzSolveEigen2", method, "B", t);
    #endif
    info=pzSolveEigen2(nev, nFull, narows, nacols, desc,
                       a, b, ev, q, work,
                       mpi_comm_world, mpi_comm_rows, mpi_comm_cols, method,
                       THIS_REAL_ELPA_KERNEL_API,
                       wantEigenVector, wantDebug);
    #ifdef _DEBUG
        timer(myid, "pzSolveEigen2", method, "B", t);
    #endif
    return info;
}

int pzDecomposeRightMatrix1(int nFull, int narows, int nacols, int *desc,
                           std::complex<double> *b, double *ev, std::complex<double> *q, std::complex<double> *work,
                           MPI_Comm mpi_comm_world, int mpi_comm_rows, int mpi_comm_cols,
                           int &method)
{
    int info=0; // for elpa functions, 1 is for success, 0 is for failure
    int wantDebug=false;
    #ifdef _DEBUG
        int myid;
        double t=-1;
        std::stringstream outlog;
        MPI_Comm_rank(MPI_COMM_WORLD, &myid);
        outlog.str("");
        outlog<<"DEBUG: Process "<<myid<<" Enter pzDecomposeRightMatrix1"<<std::endl;
        MPI_Barrier(MPI_COMM_WORLD);
        std::cout<<outlog.str();
    #endif
    switch(method)
    {
        case 0:
            #ifdef _DEBUG
                t=-1;
                timer(myid, "pzCheloskyDecomposeRightMatrix", method, "A", t);
            #endif
            info=pzCheloskyDecomposeRightMatrix(nFull, narows, nacols, desc, b,
                                            mpi_comm_world, mpi_comm_rows, mpi_comm_cols, method, wantDebug);
            #ifdef _DEBUG
                timer(myid, "pzCheloskyDecomposeRightMatrix", method, "A", t);
            #endif
            if(info!=0)
            {
                method=3;
                #ifdef _DEBUG
                    t=-1;
                    timer(myid, "pzDiagonalizeRightMatrix1", method, "A3", t);
                #endif
                info=pzDiagonalizeRightMatrix1(nFull, narows, nacols, desc,
                                            b, ev, q, work,
                                            mpi_comm_rows, mpi_comm_cols);
                #ifdef _DEBUG
                    timer(myid, "pzDiagonalizeRightMatrix1", method, "A3", t);
                #endif
            }
            break;
        case 1:
        case 2:
            #ifdef _DEBUG
                t=-1;
                timer(myid, "pzCheloskyDecomposeRightMatrix", method, "A", t);
            #endif
            info=pzCheloskyDecomposeRightMatrix(nFull, narows, nacols, desc, b,
                                                mpi_comm_world, mpi_comm_rows, mpi_comm_cols, method, wantDebug);
            #ifdef _DEBUG
                timer(myid, "pzCheloskyDecomposeRightMatrix", method, "A", t);
            #endif
            break;
        case 3:
            #ifdef _DEBUG
                t=-1;
                timer(myid, "pzDiagonalizeRightMatrix1", method, "A3", t);
            #endif
            info=pzDiagonalizeRightMatrix1(nFull, narows, nacols, desc,
                                        b, ev, q, work,
                                        mpi_comm_rows, mpi_comm_cols);
            #ifdef _DEBUG
                timer(myid, "pzDiagonalizeRightMatrix1", method, "A3", t);
            #endif
            break;
        default:
            info=10;
    }
    return info;
}

int pzDecomposeRightMatrix2(int nFull, int narows, int nacols, int *desc,
                           std::complex<double> *b, double *ev, std::complex<double> *q, std::complex<double> *work,
                           MPI_Comm mpi_comm_world, int mpi_comm_rows, int mpi_comm_cols,
                           int &method, int THIS_REAL_ELPA_KERNEL_API)
{
    int info=0; // for elpa functions, 1 is for success, 0 is for failure
    int wantDebug=false;
    //int allinfo=0;

    #ifdef _DEBUG
        int myid;
        double t;
        std::stringstream outlog;
        MPI_Comm_rank(MPI_COMM_WORLD, &myid);
        outlog.str("");
        outlog<<"DEBUG: Process "<<myid<<" Enter pzDecomposeRightMatrix2"<<std::endl;
        MPI_Barrier(MPI_COMM_WORLD);
        std::cout<<outlog.str();
    #endif
    switch(method)
    {
        case 0:
            #ifdef _DEBUG
                t=-1;
                timer(myid, "pzCheloskyDecomposeRightMatrix", method, "A", t);
            #endif
            info=pzCheloskyDecomposeRightMatrix(nFull, narows, nacols, desc, b,
                                            		mpi_comm_world, mpi_comm_rows, mpi_comm_cols,
                                            		method, wantDebug);
            #ifdef _DEBUG
                timer(myid, "pzCheloskyDecomposeRightMatrix", method, "A", t);
            #endif
            if(info!=0)
            {
                method=3;
                #ifdef _DEBUG
                    t=-1;
                    timer(myid, "pzDiagonalizeRightMatrix2", method, "A3", t);
                #endif
                info=pzDiagonalizeRightMatrix2(nFull, narows, nacols, desc,
                                            b, ev, q, work,
                                            mpi_comm_world, mpi_comm_rows, mpi_comm_cols,
                                            THIS_REAL_ELPA_KERNEL_API);
                #ifdef _DEBUG
                    timer(myid, "pzDiagonalizeRightMatrix2", method, "A3", t);
                #endif
            }
            break;
        case 1:
        case 2:
            #ifdef _DEBUG
                t=-1;
                timer(myid, "pzCheloskyDecomposeRightMatrix", method, "A", t);
            #endif
            info=pzCheloskyDecomposeRightMatrix(nFull, narows, nacols, desc, b,
                                            	  mpi_comm_world, mpi_comm_rows, mpi_comm_cols,
                                            	  method, wantDebug);
            #ifdef _DEBUG
                timer(myid, "pzCheloskyDecomposeRightMatrix", method, "A", t);
            #endif
            break;
        case 3:
            #ifdef _DEBUG
                t=-1;
                timer(myid, "pzDiagonalizeRightMatrix2", method, "A3", t);
            #endif
            info=pzDiagonalizeRightMatrix2(nFull, narows, nacols, desc,
                                           b, ev, q, work,
                                           mpi_comm_world, mpi_comm_rows, mpi_comm_cols,
                                           THIS_REAL_ELPA_KERNEL_API);
//						int pzDiagonalizeRightMatrix2(int nFull, int narows, int nacols, int *desc,
//						                              std::complex<doubles> *b, double *ev, std::complex<double> *q, std::complex<double> *work,
//						                              MPI_Comm mpi_comm_world, int mpi_comm_rows, int mpi_comm_cols,
//						                              int THIS_REAL_ELPA_KERNEL_API);
            #ifdef _DEBUG
                timer(myid, "pzDiagonalizeRightMatrix2", method, "A3", t);
            #endif
            break;
        default:
            info=10;
    }
    return info;
}

int pzCheloskyDecomposeRightMatrix(int nFull, int narows, int nacols, int *desc, std::complex<double> *b,
                                   MPI_Comm mpi_comm_world, int mpi_comm_rows, int mpi_comm_cols,
                                   int &method,
                                   bool wantDebug)
{
    int info, allinfo;
    int my_blacs_ctxt;
    int nprows, npcols, myprow, mypcol, nblk;
    char uplo;
    int isrc=1, jsrc=1;
    int lda;

    #ifdef _DEBUG
        int myid;
        double t;
        std::stringstream outlog;
        MPI_Comm_rank(MPI_COMM_WORLD, &myid);
        outlog.str("");
        outlog<<"DEBUG: Process "<<myid<<" Enter pzCheloskyDecomposeRightMatrix"<<std::endl;
        MPI_Barrier(MPI_COMM_WORLD);
        std::cout<<outlog.str();
    #endif
    //info=1; // for elpa functions, 1 is for success, 0 is for failure
    my_blacs_ctxt=desc[1];
    nblk=desc[4];
    lda=desc[8];
    Cblacs_gridinfo(my_blacs_ctxt, &nprows, &npcols, &myprow, &mypcol);

    // do cholesky decomposing of overlap matrix B: B=U^T*U, put U to b
    if(method==1 || (method==0 && nFull<CHOLESKY_CRITICAL_SIZE))
    {
        method=1;
        uplo='U';
        #ifdef _DEBUG
            t=-1;
            timer(myid, "pzpotrf_", method, "1", t);
        #endif
        pzpotrf_(&uplo, &nFull, b, &isrc, &jsrc, desc, &info);
        #ifdef _DEBUG
            timer(myid, "pzpotrf_", method, "1", t);
        #endif
        MPI_Allreduce(&info, &allinfo, 1, MPI_INT, MPI_SUM, mpi_comm_world);
        if(allinfo != 0)
        {
            if(method==1)
            {
                info=10;
                return info;
            }
            else
            {
                method=2;
                #ifdef _DEBUG
                    t=-1;
                    timer(myid, "elpa_cholesky_complex", method, "1", t);
                #endif
                info=elpa_cholesky_complex(nFull, b, narows, nblk, nacols, mpi_comm_rows, mpi_comm_cols, wantDebug);
                #ifdef _DEBUG
                    timer(myid, "elpa_cholesky_complex", method, "1", t);
                #endif
                MPI_Allreduce(&info, &allinfo, 1, MPI_INT, MPI_MIN, mpi_comm_world);
                if(allinfo != 1)
                {
                    info=10;
                    return info;
                }
            }
        }
    }
    else if(method==2 || (method==0 && nFull>CHOLESKY_CRITICAL_SIZE))
    {
        method=2;
        #ifdef _DEBUG
            t=-1;
            timer(myid, "elpa_cholesky_complex", method, "1", t);
        #endif
        wantDebug=true;
        info=elpa_cholesky_complex(nFull, b, narows, nblk, nacols, mpi_comm_rows, mpi_comm_cols, wantDebug);
        #ifdef _DEBUG
            timer(myid, "elpa_cholesky_complex", method, "1", t);
        #endif
        MPI_Allreduce(&info, &allinfo, 1, MPI_INT, MPI_MIN, mpi_comm_world);
        if(allinfo != 1)
        {
            if(method==2)
            {
                info=10;
                return info;
            }
            else
            {
                method=1;
                uplo='U';
                #ifdef _DEBUG
                    t=-1;
                    timer(myid, "pzpotrf_", method, "1", t);
                #endif
                pzpotrf_(&uplo, &nFull, b, &isrc, &jsrc, desc, &info);
                #ifdef _DEBUG
                    timer(myid, "pzpotrf_", method, "1", t);
                #endif
                MPI_Allreduce(&info, &allinfo, 1, MPI_INT, MPI_SUM, mpi_comm_world);
                if(allinfo != 0)
                {
                    info=10;
                    return info;
                }
            }
        }
    }
    else
    {
        info=20;
        return info;
    }

    //clear_lower_part of b
    #ifdef _DEBUG
        t=-1;
        timer(myid, "clear_lower_part", method, "2", t);
    #endif
    for(int j=0; j<nacols; ++j)
    {
        int jGlobal=Local_Orbital_wfc::globalIndex(j, nblk, npcols, mypcol);
        for(int i=0; i<narows; ++i)
        {
            int iGlobal=Local_Orbital_wfc::globalIndex(i, nblk, nprows, myprow);
            if(iGlobal>jGlobal) b[i+j*narows]=0;
        }
    }
    #ifdef _DEBUG
        timer(myid, "clear_lower_part", method, "2", t);
    #endif
    // invert U to U^-1, put to b
    #ifdef _DEBUG
        t=-1;
        timer(myid, "elpa_invert_trm_complex", method, "3", t);
    #endif
    info=elpa_invert_trm_complex(nFull, b, lda, nblk, nacols, mpi_comm_rows, mpi_comm_cols, wantDebug);
    #ifdef _DEBUG
        timer(myid, "elpa_invert_trm_complex", method, "3", t);
    #endif
    if(info != 1)
    {
        info=20;
        return info;
    }

    return 0;
}


int pzDiagonalizeRightMatrix1(int nFull, int narows, int nacols, int *desc,
                              std::complex<double> *b, double *ev, std::complex<double> *q, std::complex<double> *work,
                              int mpi_comm_rows, int mpi_comm_cols)
{
    int info;
    int my_blacs_ctxt;
    int nprows, npcols, myprow, mypcol, nblk;
    char transa, transb;
    int isrc=1, jsrc=1;
    double alpha, beta;
    int lda;

    #ifdef _DEBUG
        int myid;
        double t;
        int method=3;
        std::stringstream outlog;
        MPI_Comm_rank(MPI_COMM_WORLD, &myid);
        outlog.str("");
        outlog<<"DEBUG: Process "<<myid<<" Enter pzDiagonalizeRightMatrix1"<<std::endl;
        MPI_Barrier(MPI_COMM_WORLD);
        std::cout<<outlog.str();
    #endif
    info=1; // for elpa functions, 1 is for success, 0 is for failure
    my_blacs_ctxt=desc[1];
    nblk=desc[4];
    lda=desc[8];
    Cblacs_gridinfo(my_blacs_ctxt, &nprows, &npcols, &myprow, &mypcol);

    // calculate the eigenvectors of overlap matrix b, put to q
    // calculate the eigenvalues of overlap matrix b, put to ev
    #ifdef _DEBUG
        t=-1;
        timer(myid, "elpa_solve_evp_complex_1stage", method, "1", t);
    #endif
    info=elpa_solve_evp_complex_1stage(nFull, nFull, b, lda, ev, q, lda, nblk,
                               nacols, mpi_comm_rows, mpi_comm_cols);
    #ifdef _DEBUG
        timer(myid, "elpa_solve_evp_complex_1stage", method, "1", t);
    #endif
    if(info != 1)
    {
        info=100;
        return info;
    }

    // calculate B^{-1/2}_{i,j}=\sum_k q_{i,k}*ev_k^{-1/2}*q_{j,k} and put to b, which will be b^-1/2
    // calculate q*ev and put to work
    #ifdef _DEBUG
        t=-1;
        timer(myid, "q*ev", method, "2", t);
    #endif
    for(int i=0; i<nacols; ++i)
    {
        int eidx=Local_Orbital_wfc::globalIndex(i, nblk, npcols, mypcol);
        double ev_sqrt=ev[eidx]>DBL_MIN?1.0/sqrt(ev[eidx]):0;
        for(int j=0; j<narows; ++j)
            work[i*lda+j]=q[i*lda+j]*ev_sqrt;
    }
    #ifdef _DEBUG
        timer(myid, "q*ev", method, "2", t);
    #endif
    // calculate qevq=qev*q^T, put to b, which is B^{-1/2}
    transa='N';
    transb='C';
    alpha=1.0;
    beta=0.0;
    #ifdef _DEBUG
        t=-1;
        timer(myid, "pzgemm_", method, "3", t);
    #endif
    pzgemm_(&transa, &transb, &nFull, &nFull, &nFull,
            &alpha, work, &isrc, &jsrc, desc,
                    q, &isrc, &jsrc, desc,
            &beta,  b, &isrc, &jsrc, desc);
    #ifdef _DEBUG
        timer(myid, "pzgemm_", method, "3", t);
    #endif
    return 0;
}

int pzDiagonalizeRightMatrix2(int nFull, int narows, int nacols, int *desc,
                              std::complex<double> *b, double *ev, std::complex<double> *q, std::complex<double> *work,
                              MPI_Comm mpi_comm_world, int mpi_comm_rows, int mpi_comm_cols,
                              int THIS_REAL_ELPA_KERNEL_API)
{
    int info;
    int my_blacs_ctxt;
    int nprows, npcols, myprow, mypcol, nblk;
    char transa, transb;
    int isrc=1, jsrc=1;
    double alpha, beta;
    int lda;

    #ifdef _DEBUG
        int myid;
        double t;
        int method=3;
        std::stringstream outlog;
        MPI_Comm_rank(MPI_COMM_WORLD, &myid);
        outlog.str("");
        outlog<<"DEBUG: Process "<<myid<<" Enter pdDiagonalizeRightMatrix2"<<std::endl;
        MPI_Barrier(MPI_COMM_WORLD);
        std::cout<<outlog.str();
    #endif
    info=1; // for elpa functions, 1 is for success, 0 is for failure
    my_blacs_ctxt=desc[1];
    nblk=desc[4];
    lda=desc[8];
    Cblacs_gridinfo(my_blacs_ctxt, &nprows, &npcols, &myprow, &mypcol);

    // calculate the eigenvectors of overlap matrix b, put to q
    // calculate the eigenvalues of overlap matrix b, put to ev
    #ifdef _DEBUG
        t=-1;
        timer(myid, "elpa_solve_evp_complex_2stage", method, "1", t);
    #endif
    info=elpa_solve_evp_complex_2stage(nFull, nFull, b, lda, ev, q, lda, nblk, nacols,
                                    mpi_comm_rows, mpi_comm_cols, mpi_comm_world,
                                    THIS_REAL_ELPA_KERNEL_API);
    #ifdef _DEBUG
        timer(myid, "elpa_solve_evp_complex_2stage", method, "1", t);
    #endif
    if(info != 1)
    {
        info=100;
        return info;
    }

    // calculate B^{-1/2}_{i,j}=\sum_k q_{i,k}*ev_k^{-1/2}*q_{j,k} and put to b, which will be b^-1/2
    // calculate q*ev and put to work
    #ifdef _DEBUG
        t=-1;
        timer(myid, "q*ev", method, "2", t);
    #endif
    for(int i=0; i<nacols; ++i)
    {
        int eidx=Local_Orbital_wfc::globalIndex(i, nblk, npcols, mypcol);
        //double ev_sqrt=1.0/sqrt(ev[eidx]);
        double ev_sqrt=ev[eidx]>DBL_MIN?1.0/sqrt(ev[eidx]):0;
        for(int j=0; j<narows; ++j)
            work[i*lda+j]=q[i*lda+j]*ev_sqrt;
    }
    #ifdef _DEBUG
        timer(myid, "q*ev", method, "2", t);
    #endif
    // calculate qevq=qev*q^T, put to b, which is B^{-1/2}
    transa='N';
    transb='C';
    alpha=1.0;
    beta=0.0;
    #ifdef _DEBUG
        t=-1;
        timer(myid, "pzgemm_", method, "3", t);
    #endif
    pzgemm_(&transa, &transb, &nFull, &nFull, &nFull,
            &alpha, work, &isrc, &jsrc, desc,
                    q, &isrc, &jsrc, desc,
            &beta,  b, &isrc, &jsrc, desc);
    #ifdef _DEBUG
        timer(myid, "pzgemm_", method, "3", t);
    #endif
    return 0;
}

int pzSolveEigen1(int nev, int nFull, int narows, int nacols, int *desc,
                 std::complex<double> *a, std::complex<double> *b, double *ev, std::complex<double> *q, std::complex<double> *work,
                 int mpi_comm_rows, int mpi_comm_cols, int method,
                 bool wantEigenVector, bool wantDebug)
{
    int info;
    int my_blacs_ctxt;
    int nprows, npcols, myprow, mypcol, nblk;
    char transa, transb, side, uplo, diag;
    int isrc=1, jsrc=1;
    double alpha, beta;
    int lda;

    #ifdef _DEBUG
        int myid;
        double t;
        std::stringstream outlog;
        MPI_Comm_rank(MPI_COMM_WORLD, &myid);
        outlog.str("");
        outlog<<"DEBUG: Process "<<myid<<" Enter pzSolveEigen1"<<std::endl;
        MPI_Barrier(MPI_COMM_WORLD);
        std::cout<<outlog.str();
    #endif

    info=1; // for elpa functions, 1 is for success, 0 is for failure
    my_blacs_ctxt=desc[1];
    nblk=desc[4];
    lda=desc[8];
    Cblacs_gridinfo(my_blacs_ctxt, &nprows, &npcols, &myprow, &mypcol);

    if(method==1 || method==2) // use cholesky decomposing result of overlap matrix (b)
    {
        // calculate A*U^-1, put to work
        transa='C';
        transb='N';
        alpha=1.0;
        beta=0.0;
        #ifdef _DEBUG
            t=-1;
            timer(myid, "pzgemm_", method, "4", t);
        #endif
        pzgemm_(&transa, &transb, &nFull, &nFull, &nFull,
                &alpha, a, &isrc, &isrc, desc,
                        b, &isrc, &isrc, desc,
                &beta,  work, &isrc, &isrc, desc);
        #ifdef _DEBUG
            timer(myid, "pzgemm_", method, "4", t);
        #endif
        // calculate U^-T^(A*U^-1), put to a
        transa='C';
        transb='N';
        alpha=1.0;
        beta=0.0;
        #ifdef _DEBUG
            t=-1;
            timer(myid, "pzgemm_", method, "5", t);
        #endif
        pzgemm_(&transa, &transb, &nFull, &nFull, &nFull,
                &alpha, b, &isrc, &isrc, desc,
                        work, &isrc, &isrc, desc,
                &beta,  a, &isrc, &isrc, desc);
        #ifdef _DEBUG
            timer(myid, "pzgemm_", method, "5", t);
        #endif
        // calculate the eigenvalues and eigenvectors, put to ev and q
        #ifdef _DEBUG
            t=-1;
            timer(myid, "elpa_solve_evp_complex_1stage", method, "6", t);
        #endif
        info=elpa_solve_evp_complex_1stage(nFull, nev, a, lda, ev, q, lda, nblk,
                                           nacols, mpi_comm_rows, mpi_comm_cols);
        #ifdef _DEBUG
            timer(myid, "elpa_solve_evp_complex_1stage", method, "6", t);
        #endif
        if(info != 1)
        {
            info=100;
            return info;
        }
        if(wantEigenVector)
        {
            // transform the eigenvectors to original general equation, let U^-1*q, and put to q
            side='L';
            uplo='U';
            transa='N';
            diag='N';
            #ifdef _DEBUG
                t=-1;
                timer(myid, "pztrmm_", method, "7", t);
            #endif
            pztrmm_(&side, &uplo, &transa,  &diag, &nFull, &nFull,
                    &alpha, b, &isrc, &jsrc, desc,
                            q, &isrc, &jsrc, desc);
            #ifdef _DEBUG
                timer(myid, "pztrmm_", method, "7", t);
            #endif
        }
    }
    else if(method==3)  //use diagonalized result of overlap matrix (b)
    {
        // calculate b*a^T and put to work
        transa='N';
        transb='C';
        alpha=1.0;
        beta=0.0;
        #ifdef _DEBUG
            t=-1;
            timer(myid, "pzgemm_", method, "4", t);
        #endif
        pzgemm_(&transa, &transb, &nFull, &nFull, &nFull,
                &alpha, b, &isrc, &jsrc, desc,
                        a, &isrc, &jsrc, desc,
                &beta,  work, &isrc, &jsrc, desc);
        #ifdef _DEBUG
            timer(myid, "pzgemm_", method, "4", t);
        #endif
        // calculate b*work^T and put to a -- origian A*x=v*B*x was transform to a*x'=v*x'
        #ifdef _DEBUG
            t=-1;
            timer(myid, "pzgemm_", method, "5", t);
        #endif
        pzgemm_(&transa, &transb, &nFull, &nFull, &nFull,
                &alpha, b, &isrc, &jsrc, desc,
                        work, &isrc, &jsrc, desc,
                &beta,  a, &isrc, &jsrc, desc);
        #ifdef _DEBUG
            timer(myid, "pzgemm_", method, "5", t);
        #endif
        // calculate eigenvectors of a*x'=v*x'
        #ifdef _DEBUG
            t=-1;
            timer(myid, "elpa_solve_evp_complex_1stage", method, "6", t);
        #endif
        info=elpa_solve_evp_complex_1stage(nFull, nFull, a, lda, ev, work, lda, nblk,
                                           nacols, mpi_comm_rows, mpi_comm_cols);
        #ifdef _DEBUG
            timer(myid, "elpa_solve_evp_complex_1stage", method, "6", t);
        #endif
        if(info != 1)
        {
            info=200;
            return info;
        }
        if(wantEigenVector)
        {
            // transform the eigenvectors to original general equation, let b^T*q, and put to q
            transa='C';
            transb='N';
            #ifdef _DEBUG
                t=-1;
                timer(myid, "pzgemm_", method, "4", t);
            #endif
            pzgemm_(&transa, &transb, &nFull, &nFull, &nFull,
                    &alpha, b, &isrc, &jsrc, desc,
                            work, &isrc, &jsrc, desc,
                    &beta,  q, &isrc, &jsrc, desc);
            #ifdef _DEBUG
                timer(myid, "pzgemm_", method, "4", t);
            #endif
        }
    }
    else
    {
        info=30;
        return info;
    }

    return 0;
}

int pzSolveEigen2(int nev, int nFull, int narows, int nacols, int *desc,
                 std::complex<double> *a, std::complex<double> *b, double *ev, std::complex<double> *q, std::complex<double> *work,
                 MPI_Comm mpi_comm_world, int mpi_comm_rows, int mpi_comm_cols, int method,
                 int THIS_REAL_ELPA_KERNEL_API,
                 bool wantEigenVector, bool wantDebug)
{
    int info;
    int my_blacs_ctxt;
    int nprows, npcols, myprow, mypcol, nblk;
    char transa, transb, side, uplo, diag;
    int isrc=1, jsrc=1;
    double alpha, beta;
    int lda;

    #ifdef _DEBUG
        int myid;
        double t;
        std::stringstream outlog;
        MPI_Comm_rank(MPI_COMM_WORLD, &myid);
        outlog.str("");
        outlog<<"DEBUG: Process "<<myid<<" Enter pzSolveEigen2"<<std::endl;
        MPI_Barrier(MPI_COMM_WORLD);
        std::cout<<outlog.str();
    #endif

    info=1; // for elpa functions, 1 is for success, 0 is for failure
    my_blacs_ctxt=desc[1];
    nblk=desc[4];
    lda=desc[8];
    Cblacs_gridinfo(my_blacs_ctxt, &nprows, &npcols, &myprow, &mypcol);

    if(method==1||method==2) // use chmelosky decomposing result of overlap matrix (b)
    {
        // calculate A*U^-1, put to work
        transa='C';
        transb='N';
        alpha=1.0;
        beta=0.0;
        #ifdef _DEBUG
            t=-1;
            timer(myid, "pzgemm_", method, "4", t);
        #endif
        pzgemm_(&transa, &transb, &nFull, &nFull, &nFull,
                &alpha, a, &isrc, &isrc, desc,
                        b, &isrc, &isrc, desc,
                &beta,  work, &isrc, &isrc, desc);
        #ifdef _DEBUG
            timer(myid, "pzgemm_", method, "4", t);
        #endif
        // calculate U^-T^(A*U^-1), put to a
        transa='C';
        transb='N';
        alpha=1.0;
        beta=0.0;
        #ifdef _DEBUG
            t=-1;
            timer(myid, "pzgemm_", method, "5", t);
        #endif
        pzgemm_(&transa, &transb, &nFull, &nFull, &nFull,
                &alpha, b, &isrc, &isrc, desc,
                        work, &isrc, &isrc, desc,
                &beta,  a, &isrc, &isrc, desc);
        #ifdef _DEBUG
            timer(myid, "pzgemm_", method, "5", t);
        #endif
        // calculate the eigenvalues and eigenvectors, put to ev and q
        #ifdef _DEBUG
            t=-1;
            timer(myid, "elpa_solve_evp_complex_2stage", method, "6", t);
        #endif
        info=elpa_solve_evp_complex_2stage(nFull, nev, a, lda, ev, q, lda, nblk, nacols,
                                        mpi_comm_rows, mpi_comm_cols, mpi_comm_world,
                                        THIS_REAL_ELPA_KERNEL_API);
        #ifdef _DEBUG
            timer(myid, "elpa_solve_evp_complex_2stage", method, "6", t);
        #endif
        if(info != 1)
        {
            info=100;
            return info;
        }
        if(wantEigenVector)
        {
            // transform the eigenvectors to original general equation, let U^-1*q, and put to q
            side='L';
            uplo='U';
            transa='N';
            diag='N';
            #ifdef _DEBUG
                t=-1;
                timer(myid, "pztrmm_", method, "7", t);
            #endif
            pztrmm_(&side, &uplo, &transa,  &diag, &nFull, &nFull,
                    &alpha, b, &isrc, &jsrc, desc,
                            q, &isrc, &jsrc, desc);
            #ifdef _DEBUG
                timer(myid, "pztrmm_", method, "7", t);
            #endif
        }
    }
    else if(method==3)  //use diagonalized result of overlap matrix (b)
    {
        // calculate b*a^T and put to work
        transa='N';
        transb='C';
        alpha=1.0;
        beta=0.0;
        #ifdef _DEBUG
            t=-1;
            timer(myid, "pzgemm_", method, "4", t);
        #endif
        pzgemm_(&transa, &transb, &nFull, &nFull, &nFull,
                &alpha, b, &isrc, &jsrc, desc,
                        a, &isrc, &jsrc, desc,
                &beta,  work, &isrc, &jsrc, desc);
        #ifdef _DEBUG
            timer(myid, "pzgemm_", method, "4", t);
        #endif
        // calculate b*work^T and put to a -- origian A*x=v*B*x was transform to a*x'=v*x'
        #ifdef _DEBUG
            t=-1;
            timer(myid, "pzgemm_", method, "5", t);
        #endif
        pzgemm_(&transa, &transb, &nFull, &nFull, &nFull,
                &alpha, b, &isrc, &jsrc, desc,
                        work, &isrc, &jsrc, desc,
                &beta,  a, &isrc, &jsrc, desc);
        #ifdef _DEBUG
            timer(myid, "pzgemm_", method, "5", t);
        #endif
        // calculate eigenvectors of a*x'=v*x'
        #ifdef _DEBUG
            t=-1;
            timer(myid, "elpa_solve_evp_complex_2stage", method, "6", t);
        #endif
        info=elpa_solve_evp_complex_2stage(nFull, nev, a, lda, ev, work, lda, nblk, nacols,
                                           mpi_comm_rows, mpi_comm_cols, mpi_comm_world,
                                           THIS_REAL_ELPA_KERNEL_API);
        #ifdef _DEBUG
            timer(myid, "elpa_solve_evp_complex_2stage", method, "6", t);
        #endif
        if(info != 1)
        {
            info=200;
            return info;
        }
        if(wantEigenVector)
        {
            // transform the eigenvectors to original general equation, let b^T*q, and put to q
            transa='C';
            transb='N';
            #ifdef _DEBUG
                t=-1;
                timer(myid, "pzgemm_", method, "7", t);
            #endif
            pzgemm_(&transa, &transb, &nFull, &nFull, &nFull,
                    &alpha, b, &isrc, &jsrc, desc,
                            work, &isrc, &jsrc, desc,
                    &beta,  q, &isrc, &jsrc, desc);
            #ifdef _DEBUG
                timer(myid, "pzgemm_", method, "7", t);
            #endif
        }
    }
    else
    {
        info=30;
        return info;
    }
    return 0;
}