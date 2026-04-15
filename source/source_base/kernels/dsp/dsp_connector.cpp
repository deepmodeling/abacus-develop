#include "dsp_connector.h"

#include <complex>
#include <iostream>

extern "C"
{
#define complex_double ignore_complex_double
#include <mt_hthread_blas.h> // include faster mtblas kernels
#undef complex_double
#include <mtblas_interface.h> // include normal mtblas kernels that automatically operate memory, but slower.
}
namespace mtfunc
{
void dspInitHandle(int id)
{
    mt_blas_init(id);
    std::cout << " ** DSP inited on cluster " << id << " **" << std::endl;
} // Use this at the beginning of the program to start a dsp cluster

void dspDestoryHandle(int id)
{
    hthread_dev_close(id);
    std::cout << " ** DSP closed on cluster " << id << " **" << std::endl;
} // Close dsp cluster at the end of the program

// MTBlas secretly removed its MTBLAS_TRANSPOSE data type and used the original CBLAS_TRANSPOSE. So this function is modified.

CBLAS_TRANSPOSE convertBLASTranspose(const char* blasTrans)
{
    switch (blasTrans[0])
    {
    case 'N':
    case 'n':
        return CblasNoTrans;
    case 'T':
    case 't':
        return CblasTrans;
    case 'C':
    case 'c':
        return CblasConjTrans;
    default:
        std::cout << "Invalid BLAS transpose parameter!! Use default instead." << std::endl;
        return CblasNoTrans;
    }
} // Used to convert normal transpost char to cblas transpose flag

void* malloc_ht(size_t bytes, int cluster_id)
{
    void* ptr = hthread_malloc((int)cluster_id, bytes, HT_MEM_RW);
    return ptr;
} // Malloc on dsp. Used to replace original malloc



void free_ht(void* ptr)
{
    hthread_free(ptr);
} // Free on dsp. Used to replace original free


void sgemm_mt_(const char* transa,
               const char* transb,
               const int* m,
               const int* n,
               const int* k,
               const float* alpha,
               const float* a,
               const int* lda,
               const float* b,
               const int* ldb,
               const float* beta,
               float* c,
               const int* ldc,
               int cluster_id)
{
    mtblas_sgemm(CBLAS_ORDER::CblasColMajor,
                 convertBLASTranspose(transa),
                 convertBLASTranspose(transb),
                 *m,
                 *n,
                 *k,
                 *alpha,
                 a,
                 *lda,
                 b,
                 *ldb,
                 *beta,
                 c,
                 *ldc,
                 cluster_id);
} // zgemm that needn't malloc_ht or free_ht

void dgemm_mt_(const char* transa,
               const char* transb,
               const int* m,
               const int* n,
               const int* k,
               const double* alpha,
               const double* a,
               const int* lda,
               const double* b,
               const int* ldb,
               const double* beta,
               double* c,
               const int* ldc,
               int cluster_id)
{
    mtblas_dgemm(CBLAS_ORDER::CblasColMajor,
                 convertBLASTranspose(transa),
                 convertBLASTranspose(transb),
                 *m,
                 *n,
                 *k,
                 *alpha,
                 a,
                 *lda,
                 b,
                 *ldb,
                 *beta,
                 c,
                 *ldc,
                 cluster_id);
} // cgemm that needn't malloc_ht or free_ht

void zgemm_mt_(const char* transa,
               const char* transb,
               const int* m,
               const int* n,
               const int* k,
               const std::complex<double>* alpha,
               const std::complex<double>* a,
               const int* lda,
               const std::complex<double>* b,
               const int* ldb,
               const std::complex<double>* beta,
               std::complex<double>* c,
               const int* ldc,
               int cluster_id)
{
    mtblas_zgemm(CBLAS_ORDER::CblasColMajor,
                 convertBLASTranspose(transa),
                 convertBLASTranspose(transb),
                 *m,
                 *n,
                 *k,
                 (const void*)alpha,
                 (const void*)a,
                 *lda,
                 (const void*)b,
                 *ldb,
                 (const void*)beta,
                 (void*)c,
                 *ldc,
                 cluster_id);
} // zgemm that needn't malloc_ht or free_ht

void cgemm_mt_(const char* transa,
               const char* transb,
               const int* m,
               const int* n,
               const int* k,
               const std::complex<float>* alpha,
               const std::complex<float>* a,
               const int* lda,
               const std::complex<float>* b,
               const int* ldb,
               const std::complex<float>* beta,
               std::complex<float>* c,
               const int* ldc,
               int cluster_id)
{
    mtblas_cgemm(CBLAS_ORDER::CblasColMajor,
                 convertBLASTranspose(transa),
                 convertBLASTranspose(transb),
                 *m,
                 *n,
                 *k,
                 (const void*)alpha,
                 (const void*)a,
                 *lda,
                 (const void*)b,
                 *ldb,
                 (const void*)beta,
                 (void*)c,
                 *ldc,
                 cluster_id);
} // cgemm that needn't malloc_ht or free_ht

void sgemv_mt_(const char* transa,
               const int* m,
               const int* n,
               const float* alpha,
               const float* a,
               const int* lda,
               const float* x,
               const int* incx,
               const float* beta,
               float* y,
               const int* incy,
               int cluster_id)
{
    mtblas_sgemv(CBLAS_ORDER::CblasColMajor,
                 convertBLASTranspose(transa),
                 *m,
                 *n,
                 *alpha,
                 a,
                 *lda,
                 x,
                 *incx,
                 *beta,
                 y,
                 *incy,
                 cluster_id);
}

void dgemv_mt_(const char* transa,
               const int* m,
               const int* n,
               const double* alpha,
               const double* a,
               const int* lda,
               const double* x,
               const int* incx,
               const double* beta,
               double* y,
               const int* incy,
               int cluster_id)
{
    mtblas_dgemv(CBLAS_ORDER::CblasColMajor,
                 convertBLASTranspose(transa),
                 *m,
                 *n,
                 *alpha,
                 a,
                 *lda,
                 x,
                 *incx,
                 *beta,
                 y,
                 *incy,
                 cluster_id);
}

void zgemv_mt_(const char* transa,
               const int* m,
               const int* n,
               const std::complex<double>* alpha,
               const std::complex<double>* a,
               const int* lda,
               const std::complex<double>* x,
               const int* incx,
               const std::complex<double>* beta,
               std::complex<double>* y,
               const int* incy,
               int cluster_id)
{
    mtblas_zgemv(CBLAS_ORDER::CblasColMajor,
                 convertBLASTranspose(transa),
                 *m,
                 *n,
                 (const void*)alpha,
                 (const void*)a,
                 *lda,
                 (const void*)x,
                 *incx,
                 (const void*)beta,
                 (void*)y,
                 *incy,
                 cluster_id);
}

void cgemv_mt_(const char* transa,
               const int* m,
               const int* n,
               const std::complex<float>* alpha,
               const std::complex<float>* a,
               const int* lda,
               const std::complex<float>* x,
               const int* incx,
               const std::complex<float>* beta,
               std::complex<float>* y,
               const int* incy,
               int cluster_id)
{
    mtblas_cgemv(CBLAS_ORDER::CblasColMajor,
                 convertBLASTranspose(transa),
                 *m,
                 *n,
                 (const void*)alpha,
                 (const void*)a,
                 *lda,
                 (const void*)x,
                 *incx,
                 (const void*)beta,
                 (void*)y,
                 *incy,
                 cluster_id);
}

// Used to replace original free

void sgemm_mth_(const char* transa,
                const char* transb,
                const int* m,
                const int* n,
                const int* k,
                const float* alpha,
                const float* a,
                const int* lda,
                const float* b,
                const int* ldb,
                const float* beta,
                float* c,
                const int* ldc,
                int cluster_id)
{
    mt_hthread_sgemm(CBLAS_ORDER::CblasColMajor,
                     convertBLASTranspose(transa),
                     convertBLASTranspose(transb),
                     *m,
                     *n,
                     *k,
                     *alpha,
                     a,
                     *lda,
                     b,
                     *ldb,
                     *beta,
                     c,
                     *ldc,
                     cluster_id);
} // zgemm that needn't malloc_ht or free_ht

void dgemm_mth_(const char* transa,
                const char* transb,
                const int* m,
                const int* n,
                const int* k,
                const double* alpha,
                const double* a,
                const int* lda,
                const double* b,
                const int* ldb,
                const double* beta,
                double* c,
                const int* ldc,
                int cluster_id)
{
    mt_hthread_dgemm(CBLAS_ORDER::CblasColMajor,
                     convertBLASTranspose(transa),
                     convertBLASTranspose(transb),
                     *m,
                     *n,
                     *k,
                     *alpha,
                     a,
                     *lda,
                     b,
                     *ldb,
                     *beta,
                     c,
                     *ldc,
                     cluster_id);
} // cgemm that needn't malloc_ht or free_ht

void zgemm_mth_(const char* transa,
                const char* transb,
                const int* m,
                const int* n,
                const int* k,
                const std::complex<double>* alpha,
                const std::complex<double>* a,
                const int* lda,
                const std::complex<double>* b,
                const int* ldb,
                const std::complex<double>* beta,
                std::complex<double>* c,
                const int* ldc,
                int cluster_id)
{
    std::complex<double>* alp = (std::complex<double>*)malloc_ht(sizeof(std::complex<double>), cluster_id);
    *alp = *alpha;
    std::complex<double>* bet = (std::complex<double>*)malloc_ht(sizeof(std::complex<double>), cluster_id);
    *bet = *beta;
    mt_hthread_zgemm(CBLAS_ORDER::CblasColMajor,
                     convertBLASTranspose(transa),
                     convertBLASTranspose(transb),
                     *m,
                     *n,
                     *k,
                     alp,
                     a,
                     *lda,
                     b,
                     *ldb,
                     bet,
                     c,
                     *ldc,
                     cluster_id);
    free_ht(alp);
    free_ht(bet);
} // zgemm that needn't malloc_ht or free_ht

void zgemm_pack_mth_(const char* transa,
                     const char* transb,
                     const int* m,
                     const int* n,
                     const int* k,
                     const std::complex<double>* alpha,
                     const std::complex<double>* a,
                     const int* lda,
                     const std::complex<double>* b,
                     const int* ldb,
                     const std::complex<double>* beta,
                     std::complex<double>* c,
                     const int* ldc,
                     int cluster_id)
{
    const bool transa_not = (transa[0] == 'N' || transa[0] == 'n');
    const bool transb_not = (transb[0] == 'N' || transb[0] == 'n');
    // const size_t a_elems = static_cast<size_t>(*lda) * (transa_not ? static_cast<size_t>(*k) : static_cast<size_t>(*m));
    // const size_t b_elems = static_cast<size_t>(*ldb) * (transb_not ? static_cast<size_t>(*n) : static_cast<size_t>(*k));
    const size_t c_elems = static_cast<size_t>(*ldc) * static_cast<size_t>(*n);

    // std::complex<double>* A_dsp = static_cast<std::complex<double>*>(malloc_ht(a_elems * sizeof(std::complex<double>), cluster_id));
    // std::complex<double>* B_dsp = static_cast<std::complex<double>*>(malloc_ht(b_elems * sizeof(std::complex<double>), cluster_id));
    std::complex<double>* C_dsp = static_cast<std::complex<double>*>(malloc_ht(c_elems * sizeof(std::complex<double>), cluster_id));
    std::complex<double>* alp = static_cast<std::complex<double>*>(malloc_ht(sizeof(std::complex<double>), cluster_id));
    std::complex<double>* bet = static_cast<std::complex<double>*>(malloc_ht(sizeof(std::complex<double>), cluster_id));

    // memcpy(A_dsp, a, a_elems * sizeof(std::complex<double>));
    // memcpy(B_dsp, b, b_elems * sizeof(std::complex<double>));
    memcpy(C_dsp, c, c_elems * sizeof(std::complex<double>));
    *alp = *alpha;
    *bet = *beta;

    mt_hthread_zgemm(CBLAS_ORDER::CblasColMajor,
                     convertBLASTranspose(transa),
                     convertBLASTranspose(transb),
                     *m,
                     *n,
                     *k,
                     alp,
                     a,
                    //  A_dsp,
                     *lda,
                     b,
                    //  B_dsp,
                     *ldb,
                     bet,
                    //  c,
                     C_dsp,
                     *ldc,
                     cluster_id);
    memcpy(c, C_dsp, c_elems * sizeof(std::complex<double>));

    // free_ht(A_dsp);
    // free_ht(B_dsp);
    free_ht(C_dsp);
    free_ht(alp);
    free_ht(bet);
}

void cgemm_mth_(const char* transa,
                const char* transb,
                const int* m,
                const int* n,
                const int* k,
                const std::complex<float>* alpha,
                const std::complex<float>* a,
                const int* lda,
                const std::complex<float>* b,
                const int* ldb,
                const std::complex<float>* beta,
                std::complex<float>* c,
                const int* ldc,
                int cluster_id)
{
    std::complex<float>* alp = (std::complex<float>*)malloc_ht(sizeof(std::complex<float>), cluster_id);
    *alp = *alpha;
    std::complex<float>* bet = (std::complex<float>*)malloc_ht(sizeof(std::complex<float>), cluster_id);
    *bet = *beta;

    mt_hthread_cgemm(CBLAS_ORDER::CblasColMajor,
                     convertBLASTranspose(transa),
                     convertBLASTranspose(transb),
                     *m,
                     *n,
                     *k,
                     (const void*)alp,
                     (const void*)a,
                     *lda,
                     (const void*)b,
                     *ldb,
                     (const void*)bet,
                     (void*)c,
                     *ldc,
                     cluster_id);

    free_ht(alp);
    free_ht(bet);
} // cgemm that needn't malloc_ht or free_ht

void cgemm_pack_mth_(const char* transa,
                     const char* transb,
                     const int* m,
                     const int* n,
                     const int* k,
                     const std::complex<float>* alpha,
                     const std::complex<float>* a,
                     const int* lda,
                     const std::complex<float>* b,
                     const int* ldb,
                     const std::complex<float>* beta,
                     std::complex<float>* c,
                     const int* ldc,
                     int cluster_id)
{
    const bool transa_not = (transa[0] == 'N' || transa[0] == 'n');
    const bool transb_not = (transb[0] == 'N' || transb[0] == 'n');
    const size_t a_elems = static_cast<size_t>(*lda) * (transa_not ? static_cast<size_t>(*k) : static_cast<size_t>(*m));
    const size_t b_elems = static_cast<size_t>(*ldb) * (transb_not ? static_cast<size_t>(*n) : static_cast<size_t>(*k));
    const size_t c_elems = static_cast<size_t>(*ldc) * static_cast<size_t>(*n);

    std::complex<float>* A_dsp = static_cast<std::complex<float>*>(malloc_ht(a_elems * sizeof(std::complex<float>), cluster_id));
    std::complex<float>* B_dsp = static_cast<std::complex<float>*>(malloc_ht(b_elems * sizeof(std::complex<float>), cluster_id));
    std::complex<float>* C_dsp = static_cast<std::complex<float>*>(malloc_ht(c_elems * sizeof(std::complex<float>), cluster_id));
    std::complex<float>* alp = static_cast<std::complex<float>*>(malloc_ht(sizeof(std::complex<float>), cluster_id));
    std::complex<float>* bet = static_cast<std::complex<float>*>(malloc_ht(sizeof(std::complex<float>), cluster_id));

    memcpy(A_dsp, a, a_elems * sizeof(std::complex<float>));
    memcpy(B_dsp, b, b_elems * sizeof(std::complex<float>));
    memcpy(C_dsp, c, c_elems * sizeof(std::complex<float>));
    *alp = *alpha;
    *bet = *beta;

    mt_hthread_cgemm(CBLAS_ORDER::CblasColMajor,
                     convertBLASTranspose(transa),
                     convertBLASTranspose(transb),
                     *m,
                     *n,
                     *k,
                     (const void*)alp,
                     (const void*)A_dsp,
                     *lda,
                     (const void*)B_dsp,
                     *ldb,
                     (const void*)bet,
                     (void*)C_dsp,
                     *ldc,
                     cluster_id);

    memcpy(c, C_dsp, c_elems * sizeof(std::complex<float>));

    free_ht(A_dsp);
    free_ht(B_dsp);
    free_ht(C_dsp);
    free_ht(alp);
    free_ht(bet);
}

void sgemv_mth_(const char* transa,
                const int* m,
                const int* n,
                const float* alpha,
                const float* a,
                const int* lda,
                const float* x,
                const int* incx,
                const float* beta,
                float* y,
                const int* incy,
                int cluster_id)
{
    mt_hthread_sgemv(CBLAS_ORDER::CblasColMajor,
                     convertBLASTranspose(transa),
                     *m,
                     *n,
                     *alpha,
                     a,
                     *lda,
                     x,
                     *incx,
                     *beta,
                     y,
                     *incy,
                     cluster_id);
}

void dgemv_mth_(const char* transa,
                const int* m,
                const int* n,
                const double* alpha,
                const double* a,
                const int* lda,
                const double* x,
                const int* incx,
                const double* beta,
                double* y,
                const int* incy,
                int cluster_id)
{
    mt_hthread_dgemv(CBLAS_ORDER::CblasColMajor,
                     convertBLASTranspose(transa),
                     *m,
                     *n,
                     *alpha,
                     a,
                     *lda,
                     x,
                     *incx,
                     *beta,
                     y,
                     *incy,
                     cluster_id);
}

void zgemv_mth_(const char* transa,
                const int* m,
                const int* n,
                const std::complex<double>* alpha,
                const std::complex<double>* a,
                const int* lda,
                const std::complex<double>* x,
                const int* incx,
                const std::complex<double>* beta,
                std::complex<double>* y,
                const int* incy,
                int cluster_id)
{
    std::complex<double>* alp = (std::complex<double>*)malloc_ht(sizeof(std::complex<double>), cluster_id);
    *alp = *alpha;
    std::complex<double>* bet = (std::complex<double>*)malloc_ht(sizeof(std::complex<double>), cluster_id);
    *bet = *beta;

    mt_hthread_zgemv(CBLAS_ORDER::CblasColMajor,
                     convertBLASTranspose(transa),
                     *m,
                     *n,
                     (const void*)alp,
                     (const void*)a,
                     *lda,
                     (const void*)x,
                     *incx,
                     (const void*)bet,
                     (void*)y,
                     *incy,
                     cluster_id);

    free_ht(alp);
    free_ht(bet);
}

void cgemv_mth_(const char* transa,
                const int* m,
                const int* n,
                const std::complex<float>* alpha,
                const std::complex<float>* a,
                const int* lda,
                const std::complex<float>* x,
                const int* incx,
                const std::complex<float>* beta,
                std::complex<float>* y,
                const int* incy,
                int cluster_id)
{
    std::complex<float>* alp = (std::complex<float>*)malloc_ht(sizeof(std::complex<float>), cluster_id);
    *alp = *alpha;
    std::complex<float>* bet = (std::complex<float>*)malloc_ht(sizeof(std::complex<float>), cluster_id);
    *bet = *beta;

    mt_hthread_cgemv(CBLAS_ORDER::CblasColMajor,
                     convertBLASTranspose(transa),
                     *m,
                     *n,
                     (const void*)alp,
                     (const void*)a,
                     *lda,
                     (const void*)x,
                     *incx,
                     (const void*)bet,
                     (void*)y,
                     *incy,
                     cluster_id);

    free_ht(alp);
    free_ht(bet);
}
} // namespace mtfunc
