#include "lapack_op.h"

#include <cassert>
#include <algorithm>

namespace container {
namespace op {

template <typename T>
struct dngvd_op<T, DEVICE_CPU> {
    void operator()(
            const int nstart,
            const int ldh,
            const std::complex<T> *hcc,
            const std::complex<T> *scc,
            T *eigenvalue,
            std::complex<T> *vcc)
    {
        for (int i = 0; i < nstart * ldh; i++) {
            vcc[i] = hcc[i];
        }
        int info = 0;
        int lwork = 2 * nstart + nstart * nstart;
        // std::complex<T> *work = new std::complex<T>[lwork];
        // ModuleBase::GlobalFunc::ZEROS(work, lwork);
        Tensor work(DataTypeToEnum<std::complex<T>>::value, DeviceType::CpuDevice, {lwork});
        work.zero();

        int lrwork = 1 + 5 * nstart + 2 * nstart * nstart;
        // T *rwork = new T[lrwork];
        // ModuleBase::GlobalFunc::ZEROS(rwork, lrwork);
        Tensor rwork(DataTypeToEnum<T>::value, DeviceType::CpuDevice, {lrwork});
        rwork.zero();

        int liwork = 3 + 5 * nstart;
        // int *iwork = new int[liwork];
        // ModuleBase::GlobalFunc::ZEROS(iwork, liwork);
        Tensor iwork(DataType::DT_INT, DeviceType::CpuDevice, {liwork});
        iwork.zero();
        //===========================
        // calculate all eigenvalues
        //===========================
        LapackConnector::xhegvd(1, 'V', 'U', nstart, vcc, ldh, scc, ldh, eigenvalue, work.data<std::complex<T>>(), lwork, rwork.data<T>(), lrwork, iwork.data<int>(), liwork, info);

        assert(0 == info);
    }
};


template <typename T>
struct dnevx_op<T, DEVICE_CPU> {
    void operator()(
            const int nstart,
            const int ldh,
            const std::complex<T>* hcc, // hcc
            const int nbands, // nbands
            T* eigenvalue, // eigenvalue
            std::complex<T>* vcc) // vcc
    {
        std::complex<T>* aux = new std::complex<T>[nstart * ldh];
        for (int ii = 0; ii < nstart * ldh; ii++) {
            aux[ii] = hcc[ii];
        }

        int info = 0;
        int lwork = 0;
        int nb = LapackConnector::ilaenv(1, "ZHETRD", "L", nstart, -1, -1, -1);
        if (nb < 1) {
            nb = std::max(1, nstart);
        }

        if (nb == 1 || nb >= nstart) {
            lwork = 2 * nstart; // qianrui fix a bug 2021-7-25 : lwork should be at least max(1,2*n)
        } else {
            lwork = (nb + 1) * nstart;
        }
        // T *rwork = new T[7 * nstart];
        // int *iwork = new int[5 * nstart];
        // int *ifail = new int[nstart];
        // ModuleBase::GlobalFunc::ZEROS(rwork, 7 * nstart);
        // ModuleBase::GlobalFunc::ZEROS(iwork, 5 * nstart);
        // ModuleBase::GlobalFunc::ZEROS(ifail, nstart);
        // important part:
        // In davidson, the size of work is different from dnevx_op in diagH_subspace.
        // std::complex<T> *work = new std::complex<T>[2 * lwork];
        Tensor rwork(DataTypeToEnum<T>::value, DeviceType::CpuDevice, {7 * nstart});
        Tensor iwork(DataType::DT_INT, DeviceType::CpuDevice, {5 * nstart});
        Tensor ifail(DataType::DT_INT, DeviceType::CpuDevice, {nstart});
        Tensor work(DataTypeToEnum<std::complex<T>>::value, DeviceType::CpuDevice, {2 * lwork});
        // ModuleBase::GlobalFunc::ZEROS(work, lwork); // qianrui change it, only first lwork numbers are used in zhegvx
        rwork.zero();
        iwork.zero();
        ifail.zero();
        work.zero();
        // The A and B storage space is (nstart * ldh), and the data that really participates in the zhegvx
        // operation is (nstart * nstart). In this function, the data that A and B participate in the operation will
        // be extracted into the new local variables aux and bux (the internal of the function).
        // V is the output of the function, the storage space is also (nstart * ldh), and the data size of valid V
        // obtained by the zhegvx operation is (nstart * nstart) and stored in zux (internal to the function). When
        // the function is output, the data of zux will be mapped to the corresponding position of V.
        LapackConnector::xheevx(
                1, // ITYPE = 1:  A*x = (lambda)*B*x
                'V', // JOBZ = 'V':  Compute eigenvalues and eigenvectors.
                'I', // RANGE = 'I': the IL-th through IU-th eigenvalues will be found.
                'L', // UPLO = 'L':  Lower triangles of A and B are stored.
                nstart, // N = base
                aux, // A is COMPLEX*16 array  dimension (LDA, N)
                ldh, // LDA = base
                0.0, // Not referenced if RANGE = 'A' or 'I'.
                0.0, // Not referenced if RANGE = 'A' or 'I'.
                1, // IL: If RANGE='I', the index of the smallest eigenvalue to be returned. 1 <= IL <= IU <= N,
                nbands, // IU: If RANGE='I', the index of the largest eigenvalue to be returned. 1 <= IL <= IU <= N,
                0.0, // ABSTOL
                nbands, // M: The total number of eigenvalues found.  0 <= M <= N. if RANGE = 'I', M = IU-IL+1.
                eigenvalue, // W store eigenvalues
                vcc, // store eigenvector
                ldh, // LDZ: The leading dimension of the array Z.
                work.data<std::complex<T>>(),
                lwork,
                rwork.data<T>(),
                iwork.data<int>(),
                ifail.data<int>(),
                info);

        delete[] aux;

        assert(0 == info);
    }
};

template struct dngvd_op<float, DEVICE_CPU>;
template struct dngvd_op<double, DEVICE_CPU>;

template struct dnevx_op<float, DEVICE_CPU>;
template struct dnevx_op<double, DEVICE_CPU>;
} // namespace container
} // namespace op