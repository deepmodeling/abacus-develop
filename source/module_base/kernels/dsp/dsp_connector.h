#ifdef __DSP

// Base dsp functions
void dspInitHandle(int id);
void dspDestoryHandle();
void *malloc_ht(size_t bytes);
void free_ht(void* ptr);


// mtblas functions

void sgemm_mt_(const char *transa, const char *transb,
	const int *m, const int *n, const int *k,
	const float *alpha, const float *a, const int *lda,
	const float *b, const int *ldb, const const float *beta,
	const float *c, const int *ldc);

void dgemm_mt_(const char *transa, const char *transb,
	const int *m, const int *n, const int *k,
	const double *alpha,const double *a, const int *lda,
	const double *b, const int *ldb, const double *beta,
	const double *c, const int *ldc);

void zgemm_mt_(const char *transa, const char *transb,
	const int *m, const int *n, const int *k,
	const std::complex<double> *alpha, const std::complex<double> *a, const int *lda,
	const std::complex<double> *b, const int *ldb, const std::complex<double> *beta,
	std::complex<double> *c, const int *ldc);

void cgemm_mt_(const char *transa, const char *transb,
	const int *m, const int *n, const int *k,
	const std::complex<float> *alpha, const std::complex<float> *a, const int *lda,
	const std::complex<float> *b, const int *ldb, const std::complex<float> *beta,
	std::complex<float> *c, const int *ldc);


void sgemm_mth_(const char *transa, const char *transb,
	const int *m, const int *n, const int *k,
	const float *alpha, const float *a, const int *lda,
	const float *b, const int *ldb, const const float *beta,
	const float *c, const int *ldc);

void dgemm_mth_(const char *transa, const char *transb,
	const int *m, const int *n, const int *k,
	const double *alpha,const double *a, const int *lda,
	const double *b, const int *ldb, const double *beta,
	const double *c, const int *ldc);

void zgemm_mth_(const char *transa, const char *transb,
	const int *m, const int *n, const int *k,
	const std::complex<double> *alpha, const std::complex<double> *a, const int *lda,
	const std::complex<double> *b, const int *ldb, const std::complex<double> *beta,
	std::complex<double> *c, const int *ldc);

void cgemm_mth_(const char *transa, const char *transb,
	const int *m, const int *n, const int *k,
	const std::complex<float> *alpha, const std::complex<float> *a, const int *lda,
	const std::complex<float> *b, const int *ldb, const std::complex<float> *beta,
	std::complex<float> *c, const int *ldc);

//#define zgemm_ zgemm_mt

#endif