// Fallback for DSP builds without C_PACK (MT_CPACK_DIR unset).
//
// diago_lapack.cpp calls cpack_zhegvx_driver under __DSP. The symbol normally
// comes from libc_pack.a; when C_PACK is not linked, this translation unit
// provides it by forwarding to the configured LAPACK provider's zhegvx_. CMake adds this
// file only when USE_DSP is on and MT_CPACK_DIR is unset, so it never coexists
// with libc_pack.a. See issue #8247.

#include "source_base/module_external/lapack_connector.h"

#ifdef __DSP

void cpack_zhegvx_driver(const int* itype,
                         const char* jobz,
                         const char* range,
                         const char* uplo,
                         const int* n,
                         std::complex<double>* a,
                         const int* lda,
                         std::complex<double>* b,
                         const int* ldb,
                         const double* vl,
                         const double* vu,
                         const int* il,
                         const int* iu,
                         const double* abstol,
                         int* m,
                         double* w,
                         std::complex<double>* z,
                         const int* ldz,
                         std::complex<double>* work,
                         const int* lwork,
                         double* rwork,
                         int* iwork,
                         int* ifail,
                         int* info)
{
    zhegvx_(itype,
            jobz,
            range,
            uplo,
            n,
            a,
            lda,
            b,
            ldb,
            vl,
            vu,
            il,
            iu,
            abstol,
            m,
            w,
            z,
            ldz,
            work,
            lwork,
            rwork,
            iwork,
            ifail,
            info);
}

#endif
