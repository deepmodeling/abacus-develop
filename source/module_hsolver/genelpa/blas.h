#pragma once
//blas
void dcopy_(const int *n, const double *x, const int *incx, double *y, const int *incy);
void zcopy_(const int *n, const double _Complex *x, const int *incx, double _Complex *y, const int *incy);
void dgemm_(const char *transa, const char *transb, const int *m, const int *n, const int *k,
            const double *alpha, double *a, const int *lda, 
                                 double *b, const int *ldb,
            const double *beta,  double *c, const int *ldc);
void dsymm_(char *side, char *uplo, int *m, int *n, 
            const double *alpha, double *a,  int *lda,  
                                 double *b, int *ldb, 
            const double *beta,  double *c, int *ldc);
void dtrsm_(char *side, char *uplo, char *transa, char *diag, int *m, int *n,
            const double *alpha, double *a, int *lda, 
                                 double *b, int *ldb);
//void zcopy_(int *n, double _Complex *x, int *incx, double _Complex *y, int *incy);
void zgemm_(const char *transa, const char *transb, const int *m, const int *n, const int *k,
            const double _Complex *alpha, double _Complex *a, const int *lda, 
                                          double _Complex *b, const int *ldb,
            const double _Complex *beta,  double _Complex *c, const int *ldc);
void zsymm_(char *side, char *uplo, int *m, int *n, 
            const double _Complex *alpha, double _Complex *a,  int *lda,  
                                          double _Complex *b, int *ldb, 
            const double _Complex *beta,  double _Complex *c, int *ldc);
void ztrsm_(char *side, char *uplo, char *transa, char *diag, int *m, int *n,
            double _Complex *alpha, double _Complex *a, int *lda, 
                                    double _Complex *b, int *ldb);