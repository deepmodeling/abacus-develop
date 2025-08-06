#include "dsp_connector.h"

#include <complex>
#include <iostream>

extern "C"
{
#define complex_double ignore_complex_double
#include <mt_hthread_blas.h> // MTBLAS_TRANSPOSE etc
#undef complex_double
#include <mtblas_interface.h> // gemm
}
namespace mtfunc
{
std::complex<double>* alp=nullptr;
std::complex<double>* bet=nullptr;
std::complex<float>*  alp_f=nullptr;
std::complex<float>*  bet_f=nullptr;

void dspInitHandle(int id)
{
    mt_blas_init(id);
    std::cout << " ** DSP inited on cluster " << id << " **" << std::endl;
    mtfunc::alp=(std::complex<double>*)mtfunc::malloc_ht(sizeof(std::complex<double>), id);
    mtfunc::bet=(std::complex<double>*)mtfunc::malloc_ht(sizeof(std::complex<double>), id);
    mtfunc::alp_f=(std::complex<float>*)mtfunc::malloc_ht(sizeof(std::complex<float>), id);
    mtfunc::bet_f=(std::complex<float>*)mtfunc::malloc_ht(sizeof(std::complex<float>), id);
} // Use this at the beginning of the program to start a dsp cluster

void dspDestoryHandle(int id)
{
    hthread_dev_close(id);
    std::cout << " ** DSP closed on cluster " << id << " **" << std::endl;
    mtfunc::free_ht(mtfunc::alp);
    mtfunc::free_ht(mtfunc::bet);
    mtfunc::free_ht(mtfunc::alp_f);
    mtfunc::free_ht(mtfunc::bet_f);
    mtfunc::alp = nullptr;
    mtfunc::bet = nullptr;
    mtfunc::alp_f = nullptr;
    mtfunc::bet_f = nullptr;
} // Close dsp cluster at the end

MTBLAS_TRANSPOSE convertBLASTranspose(const char* blasTrans)
{
    switch (blasTrans[0])
    {
    case 'N':
    case 'n':
        return MtblasNoTrans;
    case 'T':
    case 't':
        return MtblasTrans;
    case 'C':
    case 'c':
        return MtblasConjTrans;
    default:
        std::cout << "Invalid BLAS transpose parameter!! Use default instead." << std::endl;
        return MtblasNoTrans;
    }
} // Used to convert normal transpost char to mtblas transpose flag

void* malloc_ht(size_t bytes, int cluster_id)
{
    void* ptr = hthread_malloc((int)cluster_id, bytes, HT_MEM_RW);
    return ptr;
}

// Used to replace original malloc

void free_ht(void* ptr)
{
    hthread_free(ptr);
}

// Used to replace original free

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
    mtblas_sgemm(MTBLAS_ORDER::MtblasColMajor,
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
    mtblas_dgemm(MTBLAS_ORDER::MtblasColMajor,
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
    mtblas_zgemm(MTBLAS_ORDER::MtblasColMajor,
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
    mtblas_cgemm(MTBLAS_ORDER::MtblasColMajor,
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
    mt_hthread_sgemm(MTBLAS_ORDER::MtblasColMajor,
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
    mt_hthread_dgemm(MTBLAS_ORDER::MtblasColMajor,
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
    *alp = *alpha;
    *bet = *beta;
    mt_hthread_zgemm(MTBLAS_ORDER::MtblasColMajor,
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

} // zgemm that needn't malloc_ht or free_ht

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
    alp_f = alpha;
    bet_f = beta;

    mt_hthread_cgemm(MTBLAS_ORDER::MtblasColMajor,
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
} // cgemm that needn't malloc_ht or free_ht
} // namespace mtfunc