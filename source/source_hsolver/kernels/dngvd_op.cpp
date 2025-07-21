#include "source_hsolver/kernels/dngvd_op.h"

#include <algorithm>
#include <fstream>
#include <iostream>

namespace hsolver
{

template <typename T>
struct dngvd_op<T, base_device::DEVICE_CPU>
{
    using Real = typename GetTypeReal<T>::type;
    void operator()(const base_device::DEVICE_CPU* d,
                    const int nstart,
                    const int ldh,
                    T* hcc,
                    T* scc,
                    Real* eigenvalue,
                    T* vcc)
    {
        for (int i = 0; i < nstart * ldh; i++)
        {
            vcc[i] = hcc[i];
        }

        //===========================
        // calculate all eigenvalues
        //===========================
        LapackConnector::hegvd(
            LapackConnector::ColMajor,
            1,
            'V',
            'U',
            nstart,
            vcc,
            ldh,
            scc,
            ldh,
            eigenvalue);
    }
};

template <typename T>
struct dngv_op<T, base_device::DEVICE_CPU>
{
    using Real = typename GetTypeReal<T>::type;
    void operator()(const base_device::DEVICE_CPU* d,
                    const int nbase,
                    const int ldh,
                    T* hcc,
                    T* scc,
                    Real* eigenvalue,
                    T* vcc)
    {
        for (int i = 0; i < nbase * ldh; i++)
        {
            vcc[i] = hcc[i];
        }

        //===========================
        // calculate all eigenvalues
        //===========================
        LapackConnector::hegv(LapackConnector::ColMajor, 1, 'V', 'U', nbase, vcc, ldh, scc, ldh, eigenvalue);
    }
};

template <typename T>
struct dnevx_op<T, base_device::DEVICE_CPU>
{
    using Real = typename GetTypeReal<T>::type;
    void operator()(const base_device::DEVICE_CPU* /*ctx*/,
                    const int nstart,
                    const int ldh,
                    T* hcc,     // hcc
                    const int nbands, // nbands
                    Real* eigenvalue, // eigenvalue
                    T* vcc)           // vcc
    {
        T* aux = new T[nstart * ldh];
        for (int ii = 0; ii < nstart * ldh; ii++)
        {
            aux[ii] = hcc[ii];
        }
        std::vector<int> ifail(nstart);
        int m = 0;
        // The A and B storage space is (nstart * ldh), and the data that really participates in the zhegvx
        // operation is (nstart * nstart). In this function, the data that A and B participate in the operation will
        // be extracted into the new local variables aux and bux (the internal of the function).
        // V is the output of the function, the storage space is also (nstart * ldh), and the data size of valid V
        // obtained by the zhegvx operation is (nstart * nstart) and stored in zux (internal to the function). When
        // the function is output, the data of zux will be mapped to the corresponding position of V.
        LapackConnector::heevx(
            LapackConnector::ColMajor, // matrix_layout
            'V',        // JOBZ = 'V':  Compute eigenvalues and eigenvectors.
            'I',        // RANGE = 'I': the IL-th through IU-th eigenvalues will be found.
            'L',        // UPLO = 'L':  Lower triangles of A and B are stored.
            nstart,     // N = base
            aux,        // A is COMPLEX*16 array  dimension (LDA, N)
            ldh,        // LDA = base
            0.0,        // Not referenced if RANGE = 'A' or 'I'.
            0.0,        // Not referenced if RANGE = 'A' or 'I'.
            1,          // IL: If RANGE='I', the index of the smallest eigenvalue to be returned. 1 <= IL <= IU <= N,
            nbands,     // IU: If RANGE='I', the index of the largest eigenvalue to be returned. 1 <= IL <= IU <= N,
            0.0,        // ABSTOL
            &m,     // M: The total number of eigenvalues found.  0 <= M <= N. if RANGE = 'I', M = IU-IL+1.
            eigenvalue, // W store eigenvalues
            vcc,        // store eigenvector
            ldh,        // LDZ: The leading dimension of the array Z.
            ifail.data());
    }
};

template <typename T>
struct dngvx_op<T, base_device::DEVICE_CPU>
{
    using Real = typename GetTypeReal<T>::type;
    void operator()(const base_device::DEVICE_CPU* d,
                    const int nbase,
                    const int ldh,
                    T* hcc,
                    T* scc,
                    const int m,
                    Real* eigenvalue,
                    T* vcc)
    {
        int mm = m;

        std::vector<int> ifail(nbase);

        LapackConnector::hegvx(
            LapackConnector::ColMajor, // matrix_layout
            1,     // ITYPE = 1:  A*x = (lambda)*B*x
            'V',   // JOBZ = 'V':  Compute eigenvalues and eigenvectors.
            'I',   // RANGE = 'I': the IL-th through IU-th eigenvalues will be found.
            'U',   // UPLO = 'L':  Lower triangles of A and B are stored.
            nbase, // N = base
            hcc,   // A is COMPLEX*16 array  dimension (LDA, N)
            ldh,   // LDA = base
            scc,
            ldh,
            0.0,        // Not referenced if RANGE = 'A' or 'I'.
            0.0,        // Not referenced if RANGE = 'A' or 'I'.
            1,          // IL: If RANGE='I', the index of the smallest eigenvalue to be returned. 1 <= IL <= IU <= N,
            m,          // IU: If RANGE='I', the index of the largest eigenvalue to be returned. 1 <= IL <= IU <= N,
            0.0,        // ABSTOL
            &mm,         // M: The total number of eigenvalues found.  0 <= M <= N. if RANGE = 'I', M = IU-IL+1.
            eigenvalue, // W store eigenvalues
            vcc,        // store eigenvector
            ldh,        // LDZ: The leading dimension of the array Z.
            ifail.data());
    }
};

template struct dngvd_op<std::complex<float>, base_device::DEVICE_CPU>;
template struct dngvd_op<std::complex<double>, base_device::DEVICE_CPU>;

template struct dnevx_op<std::complex<float>, base_device::DEVICE_CPU>;
template struct dnevx_op<std::complex<double>, base_device::DEVICE_CPU>;

template struct dngvx_op<std::complex<float>, base_device::DEVICE_CPU>;
template struct dngvx_op<std::complex<double>, base_device::DEVICE_CPU>;

template struct dngv_op<std::complex<float>, base_device::DEVICE_CPU>;
template struct dngv_op<std::complex<double>, base_device::DEVICE_CPU>;
#ifdef __LCAO
template struct dngvd_op<double, base_device::DEVICE_CPU>;
template struct dnevx_op<double, base_device::DEVICE_CPU>;
template struct dngvx_op<double, base_device::DEVICE_CPU>;
template struct dngv_op<double, base_device::DEVICE_CPU>;
#endif
} // namespace hsolver