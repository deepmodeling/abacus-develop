#ifndef DSP_CONNECTOR_H
#define DSP_CONNECTOR_H
#ifdef __DSP

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

#endif
#endif