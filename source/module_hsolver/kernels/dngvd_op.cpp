#include "module_hsolver/kernels/dngvd_op.h"

#include <algorithm>

namespace hsolver {

template <typename FPTYPE>
struct dngvd_op<FPTYPE, psi::DEVICE_CPU> {
    void operator()(
        const psi::DEVICE_CPU *d,
        const int nstart,
        const int ldh,
        const std::complex<FPTYPE> *hcc,
        const std::complex<FPTYPE> *scc,
        FPTYPE *eigenvalue,
        std::complex<FPTYPE> *vcc)
    {
        for (int i = 0; i < nstart * ldh; i++) {
            vcc[i] = hcc[i];
        }
        int info = 0;
        int lwork = 2 * nstart + nstart * nstart;
        std::complex<FPTYPE> *work = new std::complex<FPTYPE>[lwork];
        ModuleBase::GlobalFunc::ZEROS(work, lwork);

        int lrwork = 1 + 5 * nstart + 2 * nstart * nstart;
        FPTYPE *rwork = new FPTYPE[lrwork];
        ModuleBase::GlobalFunc::ZEROS(rwork, lrwork);

        int liwork = 3 + 5 * nstart;
        int *iwork = new int[liwork];
        ModuleBase::GlobalFunc::ZEROS(iwork, liwork);

        //===========================
        // calculate all eigenvalues
        //===========================
        LapackConnector::xhegvd(1, 'V', 'U', nstart, vcc, ldh, scc, ldh, eigenvalue, work, lwork, rwork, lrwork, iwork, liwork, info);

        if (info != 0) {
            std::cout << "Error: xhegvd failed!" << std::endl;
        }
        assert(0 == info);

        delete[] work;
        delete[] rwork;
        delete[] iwork;
    }
};


template <typename FPTYPE>
struct dnevx_op<FPTYPE, psi::DEVICE_CPU> {
    void operator()(
        const psi::DEVICE_CPU* /*ctx*/,
        const int nstart,
        const int ldh,
        const std::complex<FPTYPE>* hcc, // hcc
        const int nbands, // nbands
        FPTYPE* eigenvalue, // eigenvalue
        std::complex<FPTYPE>* vcc) // vcc
    {
        std::complex<FPTYPE>* aux = new std::complex<FPTYPE>[nstart * ldh];
        for (int ii = 0; ii < nstart * ldh; ii++) {
            aux[ii] = hcc[ii];
        }

        int info = 0;
        int lwork = -1;
        std::complex<FPTYPE> *work = new std::complex<FPTYPE>[1];
        FPTYPE *rwork = new FPTYPE[7 * nstart];
        int *iwork = new int[5 * nstart];
        int *ifail = new int[nstart];
        
        // When lwork = -1, the demension of work will be assumed
        // Assume the denmension of work by output work[0]
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
            work,
            lwork,
            rwork,
            iwork,
            ifail,
            info);
       
    
        lwork = int(work[0].real());
        delete[] work; work = new std::complex<FPTYPE>[lwork];

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
            work,
            lwork,
            rwork,
            iwork,
            ifail,
            info);

        delete[] aux;
        delete[] work;
        delete[] rwork;
        delete[] iwork;
        delete[] ifail;

        assert(0 == info);
    }
};


template struct dngvd_op<float, psi::DEVICE_CPU>;
template struct dngvd_op<double, psi::DEVICE_CPU>;

template struct dnevx_op<float, psi::DEVICE_CPU>;
template struct dnevx_op<double, psi::DEVICE_CPU>;

} // namespace hsolver