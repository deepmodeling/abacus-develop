#ifndef LAPACK_HPP
#define LAPACK_HPP
#include <iostream>
extern "C"
{
    // =================================================================================
    // gvd:
    void dsygvd_(const int* itype, const char* jobz, const char* uplo, const int* n,
        double* a, const int* lda,
        const double* b, const int* ldb, double* w,
        double* work, int* lwork,
        int* iwork, int* liwork, int* info);

    void chegvd_(const int* itype, const char* jobz, const char* uplo, const int* n,
             std::complex<float>* a, const int* lda,
             const std::complex<float>* b, const int* ldb, float* w,
             std::complex<float>* work, int* lwork, float* rwork, int* lrwork,
             int* iwork, int* liwork, int* info);

    void zhegvd_(const int* itype, const char* jobz, const char* uplo, const int* n,
                 std::complex<double>* a, const int* lda, 
                 const std::complex<double>* b, const int* ldb, double* w,
                 std::complex<double>* work, int* lwork, double* rwork, int* lrwork,
                 int* iwork, int* liwork, int* info);
    // =================================================================================

    // =================================================================================
    // evx
    void dsyevx_(const char* jobz, const char* range, const char* uplo, const int* n,
        double* a, const int* lda,
        const double* vl, const double* vu, const int* il, const int* iu, const double* abstol,
        const int* m, double* w, double* z, const int* ldz,
        double* work, const int* lwork, double* rwork, int* iwork, int* ifail, int* info);

    void cheevx_(const char* jobz, const char* range, const char* uplo, const int* n,
             std::complex<float> *a, const int* lda,
             const float* vl, const float* vu, const int* il, const int* iu, const float* abstol,
             const int* m, float* w, std::complex<float> *z, const int *ldz,
             std::complex<float> *work, const int* lwork, float* rwork, int* iwork, int* ifail, int* info);

    void zheevx_(const char* jobz, const char* range, const char* uplo, const int* n, 
                 std::complex<double> *a, const int* lda,
                 const double* vl, const double* vu, const int* il, const int* iu, const double* abstol, 
                 const int* m, double* w, std::complex<double> *z, const int *ldz, 
                 std::complex<double> *work, const int* lwork, double* rwork, int* iwork, int* ifail, int* info);
    // =================================================================================


    // =================================================================================
    // gvx
    void dsygvx_(const int* itype, const char* jobz, const char* range, const char* uplo,
        const int* n, double* A, const int* lda, double* B, const int* ldb,
        const double* vl, const double* vu, const int* il, const int* iu,
        const double* abstol, const int* m, double* w, double* Z, const int* ldz,
        double* work, const int* lwork, int* iwork, int* ifail, int* info);

    void chegvx_(const int* itype,const char* jobz,const char* range,const char* uplo,
             const int* n,std::complex<float> *a,const int* lda,std::complex<float> *b,
             const int* ldb,const float* vl,const float* vu,const int* il,
             const int* iu,const float* abstol,const int* m,float* w,
             std::complex<float> *z,const int *ldz,std::complex<float> *work,const int* lwork,
             float* rwork,int* iwork,int* ifail,int* info);

    void zhegvx_(const int* itype,const char* jobz,const char* range,const char* uplo,
                 const int* n,std::complex<double> *a,const int* lda,std::complex<double> *b,
                 const int* ldb,const double* vl,const double* vu,const int* il,
                 const int* iu,const double* abstol,const int* m,double* w,
                 std::complex<double> *z,const int *ldz,std::complex<double> *work,const int* lwork,
                 double* rwork,int* iwork,int* ifail,int* info);
    // =================================================================================

    // =================================================================================
    // gv
    void zhegv_(const int* itype,const char* jobz,const char* uplo,const int* n,
                std::complex<double>* a,const int* lda,std::complex<double>* b,const int* ldb,
                double* w,std::complex<double>* work,int* lwork,double* rwork,int* info);
    void chegv_(const int* itype,const char* jobz,const char* uplo,const int* n,
                std::complex<float>* a,const int* lda,std::complex<float>* b,const int* ldb,
                float* w,std::complex<float>* work,int* lwork,float* rwork,int* info);
	void dsygv_(const int* itype, const char* jobz,const char* uplo, const int* n,
				double* a,const int* lda,double* b,const int* ldb,
	 			double* w,double* work,int* lwork,int* info);
    // =================================================================================

}

class LapackWrapper
{
  private:
  public:
    // wrap function of fortran lapack routine zhegvd. (pointer version)
    static inline void xhegvd(const int itype,
                              const char jobz,
                              const char uplo,
                              const int n,
                              double* a,
                              const int lda,
                              const double* b,
                              const int ldb,
                              double* w,
                              double* work,
                              int lwork,
                              double* rwork,
                              int lrwork,
                              int* iwork,
                              int liwork,
                              int& info)
    {
        // call the fortran routine
        dsygvd_(&itype, &jobz, &uplo, &n, a, &lda, b, &ldb, w, work, &lwork, iwork, &liwork, &info);
    }

    // wrap function of fortran lapack routine zhegvd. (pointer version)
    static inline void xhegvd(const int itype,
                              const char jobz,
                              const char uplo,
                              const int n,
                              std::complex<float>* a,
                              const int lda,
                              const std::complex<float>* b,
                              const int ldb,
                              float* w,
                              std::complex<float>* work,
                              int lwork,
                              float* rwork,
                              int lrwork,
                              int* iwork,
                              int liwork,
                              int& info)
    {
        // call the fortran routine
        chegvd_(&itype, &jobz, &uplo, &n, a, &lda, b, &ldb, w, work, &lwork, rwork, &lrwork, iwork, &liwork, &info);
    }

    // wrap function of fortran lapack routine zhegvd.
    static inline void xhegvd(const int itype,
                              const char jobz,
                              const char uplo,
                              const int n,
                              std::complex<double>* a,
                              const int lda,
                              const std::complex<double>* b,
                              const int ldb,
                              double* w,
                              std::complex<double>* work,
                              int lwork,
                              double* rwork,
                              int lrwork,
                              int* iwork,
                              int liwork,
                              int& info)
    {
        // call the fortran routine
        zhegvd_(&itype, &jobz, &uplo, &n, a, &lda, b, &ldb, w, work, &lwork, rwork, &lrwork, iwork, &liwork, &info);
    }

    // wrap function of fortran lapack routine dsyevx.
    static inline void xheevx(const int itype,
                              const char jobz,
                              const char range,
                              const char uplo,
                              const int n,
                              double* a,
                              const int lda,
                              const double vl,
                              const double vu,
                              const int il,
                              const int iu,
                              const double abstol,
                              const int m,
                              double* w,
                              double* z,
                              const int ldz,
                              double* work,
                              const int lwork,
                              double* rwork,
                              int* iwork,
                              int* ifail,
                              int& info)
    {
        dsyevx_(&jobz,
                &range,
                &uplo,
                &n,
                a,
                &lda,
                &vl,
                &vu,
                &il,
                &iu,
                &abstol,
                &m,
                w,
                z,
                &ldz,
                work,
                &lwork,
                rwork,
                iwork,
                ifail,
                &info);
    }

    // wrap function of fortran lapack routine cheevx.
    static inline void xheevx(const int itype,
                              const char jobz,
                              const char range,
                              const char uplo,
                              const int n,
                              std::complex<float>* a,
                              const int lda,
                              const float vl,
                              const float vu,
                              const int il,
                              const int iu,
                              const float abstol,
                              const int m,
                              float* w,
                              std::complex<float>* z,
                              const int ldz,
                              std::complex<float>* work,
                              const int lwork,
                              float* rwork,
                              int* iwork,
                              int* ifail,
                              int& info)
    {
        cheevx_(&jobz,
                &range,
                &uplo,
                &n,
                a,
                &lda,
                &vl,
                &vu,
                &il,
                &iu,
                &abstol,
                &m,
                w,
                z,
                &ldz,
                work,
                &lwork,
                rwork,
                iwork,
                ifail,
                &info);
    }

    // wrap function of fortran lapack routine zheevx.
    static inline void xheevx(const int itype,
                              const char jobz,
                              const char range,
                              const char uplo,
                              const int n,
                              std::complex<double>* a,
                              const int lda,
                              const double vl,
                              const double vu,
                              const int il,
                              const int iu,
                              const double abstol,
                              const int m,
                              double* w,
                              std::complex<double>* z,
                              const int ldz,
                              std::complex<double>* work,
                              const int lwork,
                              double* rwork,
                              int* iwork,
                              int* ifail,
                              int& info)
    {
        zheevx_(&jobz,
                &range,
                &uplo,
                &n,
                a,
                &lda,
                &vl,
                &vu,
                &il,
                &iu,
                &abstol,
                &m,
                w,
                z,
                &ldz,
                work,
                &lwork,
                rwork,
                iwork,
                ifail,
                &info);
    }

    // wrap function of fortran lapack routine xhegvx ( pointer version ).
    static inline void xhegvx(const int itype,
                              const char jobz,
                              const char range,
                              const char uplo,
                              const int n,
                              std::complex<float>* a,
                              const int lda,
                              std::complex<float>* b,
                              const int ldb,
                              const float vl,
                              const float vu,
                              const int il,
                              const int iu,
                              const float abstol,
                              const int m,
                              float* w,
                              std::complex<float>* z,
                              const int ldz,
                              std::complex<float>* work,
                              const int lwork,
                              float* rwork,
                              int* iwork,
                              int* ifail,
                              int& info)
    {
        chegvx_(&itype,
                &jobz,
                &range,
                &uplo,
                &n,
                a,
                &lda,
                b,
                &ldb,
                &vl,
                &vu,
                &il,
                &iu,
                &abstol,
                &m,
                w,
                z,
                &ldz,
                work,
                &lwork,
                rwork,
                iwork,
                ifail,
                &info);
    }

    // wrap function of fortran lapack routine xhegvx ( pointer version ).
    static inline void xhegvx(const int itype,
                              const char jobz,
                              const char range,
                              const char uplo,
                              const int n,
                              std::complex<double>* a,
                              const int lda,
                              std::complex<double>* b,
                              const int ldb,
                              const double vl,
                              const double vu,
                              const int il,
                              const int iu,
                              const double abstol,
                              const int m,
                              double* w,
                              std::complex<double>* z,
                              const int ldz,
                              std::complex<double>* work,
                              const int lwork,
                              double* rwork,
                              int* iwork,
                              int* ifail,
                              int& info)
    {
        zhegvx_(&itype,
                &jobz,
                &range,
                &uplo,
                &n,
                a,
                &lda,
                b,
                &ldb,
                &vl,
                &vu,
                &il,
                &iu,
                &abstol,
                &m,
                w,
                z,
                &ldz,
                work,
                &lwork,
                rwork,
                iwork,
                ifail,
                &info);
    }
    // wrap function of fortran lapack routine xhegvx ( pointer version ).
    static inline void xhegvx(const int itype,
                              const char jobz,
                              const char range,
                              const char uplo,
                              const int n,
                              double* a,
                              const int lda,
                              double* b,
                              const int ldb,
                              const double vl,
                              const double vu,
                              const int il,
                              const int iu,
                              const double abstol,
                              const int m,
                              double* w,
                              double* z,
                              const int ldz,
                              double* work,
                              const int lwork,
                              double* rwork,
                              int* iwork,
                              int* ifail,
                              int& info)
    {
        dsygvx_(&itype, &jobz, &range, &uplo, &n, a, &lda, b, &ldb, &vl,
            &vu, &il, &iu, &abstol, &m, w, z, &ldz, work, &lwork, iwork, ifail, &info);
    }

    // wrap function of fortran lapack routine xhegvx ( pointer version ).
    static inline void xhegv(const int itype,
                             const char jobz,
                             const char uplo,
                             const int n,
                             double* a,
                             const int lda,
                             double* b,
                             const int ldb,
                             double* w,
                             double* work,
                             int lwork,
                             double* rwork,
                             int& info)
    {
        // TODO
    }

    // wrap function of fortran lapack routine xhegvx ( pointer version ).
    static inline void xhegv(const int itype,
                             const char jobz,
                             const char uplo,
                             const int n,
                             std::complex<float>* a,
                             const int lda,
                             std::complex<float>* b,
                             const int ldb,
                             float* w,
                             std::complex<float>* work,
                             int lwork,
                             float* rwork,
                             int& info)
    {
        // TODO
    }
    // wrap function of fortran lapack routine xhegvx ( pointer version ).
    static inline void xhegv(const int itype,
                             const char jobz,
                             const char uplo,
                             const int n,
                             std::complex<double>* a,
                             const int lda,
                             std::complex<double>* b,
                             const int ldb,
                             double* w,
                             std::complex<double>* work,
                             int lwork,
                             double* rwork,
                             int& info)
    {
        zhegv_(&itype, &jobz, &uplo, &n, a, &lda, b, &ldb, w, work, &lwork, rwork, &info);
    }
};
#endif // LAPACK_HPP