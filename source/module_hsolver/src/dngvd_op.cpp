#include "module_hsolver/include/dngvd_op.h"

#include <algorithm>

namespace hsolver
{

template <>
void dngvx_op<double, psi::DEVICE_CPU>::operator()(const psi::DEVICE_CPU* d,
                                                   const int nstart,
                                                   const int ldh,
                                                   const std::complex<double>* A,
                                                   const std::complex<double>* B,
                                                   const int m,
                                                   double* W,
                                                   std::complex<double>* V)
{
    int info = 0;
    int lwork = 0;
    int nb = LapackConnector::ilaenv(1, "ZHETRD", "L", nstart, -1, -1, -1);
    if (nb < 1)
    {
        nb = std::max(1, nstart);
    }

    if (nb == 1 || nb >= nstart)
    {
        lwork = 2 * nstart; // qianrui fix a bug 2021-7-25 : lwork should be at least max(1,2*n)
    }
    else
    {
        lwork = (nb + 1) * nstart;
    }
    // std::complex<double>* work = new std::complex<double>[2 * lwork];
    std::complex<double>* work = new std::complex<double>[lwork];
    double* rwork = new double[7 * nstart];
    int* iwork = new int[5 * nstart];
    int* ifail = new int[nstart];
    ModuleBase::GlobalFunc::ZEROS(work, lwork); // qianrui change it, only first lwork numbers are used in zhegvx
    ModuleBase::GlobalFunc::ZEROS(rwork, 7 * nstart);
    ModuleBase::GlobalFunc::ZEROS(iwork, 5 * nstart);
    ModuleBase::GlobalFunc::ZEROS(ifail, nstart);

    // The A and B storage space is (nstart * ldh), and the data that really participates in the zhegvx
    // operation is (nstart * nstart). In this function, the data that A and B participate in the operation will 
    // be extracted into the new local variables aux and bux (the internal of the function).
    // V is the output of the function, the storage space is also (nstart * ldh), and the data size of valid V 
    // obtained by the zhegvx operation is (nstart * nstart) and stored in zux (internal to the function). When 
    // the function is output, the data of zux will be mapped to the corresponding position of V.
    LapackConnector::zhegvx(
        1, // ITYPE = 1:  A*x = (lambda)*B*x
        'V', // JOBZ = 'V':  Compute eigenvalues and eigenvectors.
        'I', // RANGE = 'I': the IL-th through IU-th eigenvalues will be found.
        'L', // UPLO = 'L':  Lower triangles of A and B are stored.
        nstart, // N = base
        A, // A is COMPLEX*16 array  dimension (LDA, N)
        nstart, // LDA = base
        B, // B is COMPLEX*16 array, dimension (LDB, N)
        nstart, // LDB = base
        0.0, // Not referenced if RANGE = 'A' or 'I'.
        0.0, // Not referenced if RANGE = 'A' or 'I'.
        1, // IL: If RANGE='I', the index of the smallest eigenvalue to be returned. 1 <= IL <= IU <= N,
        m, // IU: If RANGE='I', the index of the largest eigenvalue to be returned. 1 <= IL <= IU <= N,
        0.0, // ABSTOL
        m, // M: The total number of eigenvalues found.  0 <= M <= N. if RANGE = 'I', M = IU-IL+1.
        W, // W store eigenvalues
        V, // store eigenvector
        nstart, // LDZ: The leading dimension of the array Z.
        work,
        lwork,
        rwork,
        iwork,
        ifail,
        info,
        ldh);

    assert(0 == info);

    delete[] work;
    delete[] rwork;
    delete[] iwork;
    delete[] ifail;
};

template <>
void dngv_op<double, psi::DEVICE_CPU>::operator()(const psi::DEVICE_CPU* d,
                                                  const int nstart,
                                                  const int ldh,
                                                  const std::complex<double>* A,
                                                  const std::complex<double>* B,
                                                  double* W,
                                                  std::complex<double>* V)
{
    int info = 0;
    int lwork = 0;
    int nb = LapackConnector::ilaenv(1, "ZHETRD", "U", nstart, -1, -1, -1);
    if (nb < 1)
    {
        nb = std::max(1, nstart);
    }
    if (nb == 1 || nb >= nstart)
    {
        lwork = 2 * nstart; // mohan modify 2009-08-02
    }
    else
    {
        lwork = (nb + 1) * nstart;
    }
    std::complex<double>* work = new std::complex<double>[lwork];
    int rwork_dim = 3 * nstart - 2;
    double* rwork = new double[rwork_dim];
    ModuleBase::GlobalFunc::ZEROS(work, lwork);
    ModuleBase::GlobalFunc::ZEROS(rwork, rwork_dim);

    for (int i = 0; i < nstart * ldh; i++)
    {
        V[i] = A[i];
    }
    
    // The A and B storage space is (nstart * ldh), and the data that really participates in the zhegvx
    // operation is (nstart * nstart). In this function, the data that A and B participate in the operation will 
    // be extracted into the new local variables aux and bux (the internal of the function).
    // V is the output of the function, the storage space is also (nstart * ldh), and the data size of valid V 
    // obtained by the zhegvx operation is (nstart * nstart) and stored in zux (internal to the function). When 
    // the function is output, the data of zux will be mapped to the corresponding position of V.
    LapackConnector::zhegv(1, 'V', 'U', nstart, V, nstart, B, nstart, W, work, lwork, rwork, info, ldh);

    assert(0 == info);

    delete[] work;
    delete[] rwork;
}

} // namespace hsolver