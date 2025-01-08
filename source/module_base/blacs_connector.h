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
// blacs
    // Initialization
#ifndef BLACS_CONNECTOR_H
#define BLACS_CONNECTOR_H

#include <complex>

extern "C"
{
	void Cblacs_pinfo(int *myid, int *nprocs);
	void Cblacs_get(int icontxt, int what, int *val);
	void Cblacs_gridmap(int* icontxt, int *usermap, int ldumap, int nprow, int npcol);
		// Informational and Miscellaneous
	void Cblacs_gridinfo(int icontxt, int* nprow, int *npcol, int *myprow, int *mypcol);
    void Cblacs_gridinit(int* icontxt, char* layout, int nprow, int npcol);
    void Cblacs_gridexit(int icontxt);
    int Cblacs_pnum(int icontxt, int prow, int pcol);
    void Cblacs_pcoord(int icontxt, int pnum, int *prow, int *pcol);
	void Cblacs_exit(int icontxt);

    // broadcast (send/recv)
    void Cigebs2d(int ConTxt, char *scope, char *top, int m, int n, int *A, int lda);
    void Cigebr2d(int ConTxt, char *scope, char *top, int m, int n, int *A, int lda, int rsrc, int csrc);

    void Csgebs2d(int ConTxt, char *scope, char *top, int m, int n, float *A, int lda);
    void Csgebr2d(int ConTxt, char *scope, char *top, int m, int n, float *A, int lda, int rsrc, int csrc);

    void Cdgebs2d(int ConTxt, char *scope, char *top, int m, int n, double *A, int lda);
    void Cdgebr2d(int ConTxt, char *scope, char *top, int m, int n, double *A, int lda, int rsrc, int csrc);

    void Ccgebs2d(int ConTxt, char *scope, char *top, int m, int n, std::complex<float> *A, int lda);
    void Ccgebr2d(int ConTxt, char *scope, char *top, int m, int n, std::complex<float> *A, int lda, int rsrc, int csrc);

    void Czgebs2d(int ConTxt, char *scope, char *top, int m, int n, std::complex<double> *A, int lda);
    void Czgebr2d(int ConTxt, char *scope, char *top, int m, int n, std::complex<double> *A, int lda, int rsrc, int csrc);
}

// unified interface for broadcast
template <typename T>
void Cxgebs2d(int ConTxt, char *scope, char *top, int m, int n, T *A, int lda)
{
    static_assert(
        std::is_same<T, int>::value ||
        std::is_same<T, float>::value ||
        std::is_same<T, double>::value ||
        std::is_same<T,std::complex<float>>::value ||
        std::is_same<T,std::complex<double>>::value,
        "Type not supported");

	if (std::is_same<T, int>::value) {
        Cigebs2d(ConTxt, scope, top, m, n, reinterpret_cast<int*>(A), lda);
    }
	if (std::is_same<T, float>::value) {
        Csgebs2d(ConTxt, scope, top, m, n, reinterpret_cast<float*>(A), lda);
    }
	if (std::is_same<T, double>::value) {
        Cdgebs2d(ConTxt, scope, top, m, n, reinterpret_cast<double*>(A), lda);
    }
	if (std::is_same<T, std::complex<float>>::value) {
        Ccgebs2d(ConTxt, scope, top, m, n, reinterpret_cast<std::complex<float>*>(A), lda);
    }
	if (std::is_same<T, std::complex<double>>::value) {
        Czgebs2d(ConTxt, scope, top, m, n, reinterpret_cast<std::complex<double>*>(A), lda);
    }
}

template <typename T>
void Cxgebr2d(int ConTxt, char *scope, char *top, int m, int n, T *A, int lda, int rsrc, int csrc)
{
    static_assert(
        std::is_same<T, int>::value ||
        std::is_same<T, float>::value ||
        std::is_same<T, double>::value ||
        std::is_same<T,std::complex<float>>::value ||
        std::is_same<T,std::complex<double>>::value,
        "Type not supported");

	if (std::is_same<T, int>::value) {
        Cigebr2d(ConTxt, scope, top, m, n, reinterpret_cast<int*>(A), lda, rsrc, csrc);
    }
	if (std::is_same<T, float>::value) {
        Csgebr2d(ConTxt, scope, top, m, n, reinterpret_cast<float*>(A), lda, rsrc, csrc);
    }
	if (std::is_same<T, double>::value) {
        Cdgebr2d(ConTxt, scope, top, m, n, reinterpret_cast<double*>(A), lda, rsrc, csrc);
    }
	if (std::is_same<T, std::complex<float>>::value) {
        Ccgebr2d(ConTxt, scope, top, m, n, reinterpret_cast<std::complex<float>*>(A), lda, rsrc, csrc);
    }
	if (std::is_same<T, std::complex<double>>::value) {
        Czgebr2d(ConTxt, scope, top, m, n, reinterpret_cast<std::complex<double>*>(A), lda, rsrc, csrc);
    }
}


#ifdef __MPI
#include <mpi.h>
extern "C"
{
    int Csys2blacs_handle(MPI_Comm SysCtxt);
    MPI_Comm Cblacs2sys_handle(int BlacsCtxt);
}
#endif // __MPI

#endif
