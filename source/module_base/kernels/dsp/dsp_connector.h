#ifndef DSP_CONNECTOR_H
#define DSP_CONNECTOR_H
#ifdef __DSP

#include "module_base/module_device/device.h"
#include "module_base/module_device/memory_op.h"
#include "module_hsolver/diag_comm_info.h"

// Base dsp functions
void dspInitHandle(int id);
void dspDestoryHandle(int id);
void *malloc_ht(size_t bytes, int cluster_id);
void free_ht(void* ptr);


// mtblas functions

void sgemm_mt_(const char *transa, const char *transb,
	const int *m, const int *n, const int *k,
	const float *alpha, const float *a, const int *lda,
	const float *b, const int *ldb, const float *beta,
	float *c, const int *ldc, int cluster_id);

void dgemm_mt_(const char *transa, const char *transb,
	const int *m, const int *n, const int *k,
	const double *alpha,const double *a, const int *lda,
	const double *b, const int *ldb, const double *beta,
	double *c, const int *ldc, int cluster_id);

void zgemm_mt_(const char *transa, const char *transb,
	const int *m, const int *n, const int *k,
	const std::complex<double> *alpha, const std::complex<double> *a, const int *lda,
	const std::complex<double> *b, const int *ldb, const std::complex<double> *beta,
	std::complex<double> *c, const int *ldc, int cluster_id);

void cgemm_mt_(const char *transa, const char *transb,
	const int *m, const int *n, const int *k,
	const std::complex<float> *alpha, const std::complex<float> *a, const int *lda,
	const std::complex<float> *b, const int *ldb, const std::complex<float> *beta,
	std::complex<float> *c, const int *ldc, int cluster_id);


void sgemm_mth_(const char *transa, const char *transb,
	const int *m, const int *n, const int *k,
	const float *alpha, const float *a, const int *lda,
	const float *b, const int *ldb, const float *beta,
	float *c, const int *ldc, int cluster_id);

void dgemm_mth_(const char *transa, const char *transb,
	const int *m, const int *n, const int *k,
	const double *alpha,const double *a, const int *lda,
	const double *b, const int *ldb, const double *beta,
	double *c, const int *ldc, int cluster_id);

void zgemm_mth_(const char *transa, const char *transb,
	const int *m, const int *n, const int *k,
	const std::complex<double> *alpha, const std::complex<double> *a, const int *lda,
	const std::complex<double> *b, const int *ldb, const std::complex<double> *beta,
	std::complex<double> *c, const int *ldc, int cluster_id);

void cgemm_mth_(const char *transa, const char *transb,
	const int *m, const int *n, const int *k,
	const std::complex<float> *alpha, const std::complex<float> *a, const int *lda,
	const std::complex<float> *b, const int *ldb, const std::complex<float> *beta,
	std::complex<float> *c, const int *ldc, int cluster_id);

//#define zgemm_ zgemm_mt

// The next is dsp utils. It may be moved to other files if this file get too huge

template <typename T>
void dsp_dav_subspace_reduce(T* hcc, T* scc, int nbase, int nbase_x, int notconv, MPI_Comm diag_comm){

	using syncmem_complex_op = base_device::memory::synchronize_memory_op<T, base_device::DEVICE_CPU, base_device::DEVICE_CPU>;

	auto* swap = new T[notconv * nbase_x];
    auto* target = new T[notconv * nbase_x];
    syncmem_complex_op()(swap, hcc + nbase * nbase_x, notconv * nbase_x);
    if (base_device::get_current_precision(swap) == "single")
    {
        MPI_Reduce(swap,
                    target,
                    notconv * nbase_x,
                    MPI_COMPLEX,
                    MPI_SUM,
                    0,
                    diag_comm);
    }
    else
    {
        MPI_Reduce(swap,
                    target,
                    notconv * nbase_x,
                    MPI_DOUBLE_COMPLEX,
                    MPI_SUM,
                    0,
                    diag_comm);
    }

    syncmem_complex_op()(hcc + nbase * nbase_x, target, notconv * nbase_x);
    syncmem_complex_op()(swap, scc + nbase * nbase_x, notconv * nbase_x);

    if (base_device::get_current_precision(swap) == "single")
    {
        MPI_Reduce(swap,
                    target,
                    notconv * nbase_x,
                    MPI_COMPLEX,
                    MPI_SUM,
                    0,
                    diag_comm);
    }
    else
    {
        MPI_Reduce(swap,
                    target,
                    notconv * nbase_x,
                    MPI_DOUBLE_COMPLEX,
                    MPI_SUM,
                    0,
                    diag_comm);
    }

    syncmem_complex_op()(scc + nbase * nbase_x, target, notconv * nbase_x);
    delete[] swap;
    delete[] target;
}


#endif
#endif