#include "cal_edm.h"
#include "source_base/module_external/scalapack_connector.h"
#include "source_base/module_external/blas_connector.h"
namespace LR
{
    // $X_{\mu i}=\sum_a c_{\mu a} X_{ai}$
    template<>
    void cal_X_ao_occ(const double* const X, const Parallel_2D& px,
        const double* const c, const Parallel_2D& pc,
        double* const X_ao_occ, const Parallel_2D& px_ao_occ)
    {
        const int nocc = px.get_global_col_size();
        const int nvirt = px.get_global_row_size();
        const int naos = pc.get_global_row_size();
        const double alpha = 1.0;
        const double beta = 0.0;
        const char transa = 'N', transb = 'N';
#ifdef __MPI
        const int i1 = 1;
        const int ivirt = nocc + 1;
        pdgemm_(&transa, &transb, &naos, &nocc, &nvirt,
            &alpha, c, &i1, &ivirt, pc.desc,
            X, &i1, &i1, px.desc,
            &beta, X_ao_occ, &i1, &i1, px_ao_occ.desc);
#else
        dgemm_(&transa, &transb, &naos, &nocc, &nvirt,
            &alpha, c + nocc * naos, &naos,
            X, &nvirt,
            &beta, X_ao_occ, &naos);
#endif
    }
    template<>
    void cal_X_ao_occ(const std::complex<double>* const X, const Parallel_2D& px,
        const std::complex<double>* const c, const Parallel_2D& pc,
        std::complex<double>* const X_ao_occ, const Parallel_2D& px_ao_occ)
    {
        const int nocc = px.get_global_col_size();
        const int nvirt = px.get_global_row_size();
        const int naos = pc.get_global_row_size();
        const std::complex<double> alpha(1.0, 0.0), beta(0.0, 0.0);
        const char transa = 'N', transb = 'N';
#ifdef __MPI
        const int i1 = 1;
        const int ivirt = nocc + 1;
        pzgemm_(&transa, &transb, &naos, &nocc, &nvirt,
            &alpha, c, &i1, &ivirt, pc.desc,
            X, &i1, &i1, px.desc,
            &beta, X_ao_occ, &i1, &i1, px_ao_occ.desc);
#else
        zgemm_(&transa, &transb, &naos, &nocc, &nvirt,
            &alpha, c + nocc * naos, &naos,
            X, &nvirt,
            &beta, X_ao_occ, &naos);
#endif
    }
    // D=X1*X2^T
    // $D_{\mu\nu} = \sum_i X1_{\mu i}X2_{\nu i}$
    template<>
    void matdot(const double* const vec1, const double* const vec2, const Parallel_2D& pvec,
        double* const dm, const Parallel_2D& pmat)
    {
        const int nocc = pvec.get_global_col_size();
        const int naos = pvec.get_global_row_size();
        const double alpha = 1.0;
        const double beta = 0.0;
        const char transa = 'N', transb = 'T';
#ifdef __MPI
        const int i1 = 1;
        pdgemm_(&transa, &transb, &naos, &naos, &nocc,
            &alpha, vec1, &i1, &i1, pvec.desc,
            vec2, &i1, &i1, pvec.desc,
            &beta, dm, &i1, &i1, pmat.desc);
#else
        dgemm_(&transa, &transb, &naos, &naos, &nocc,
            &alpha, vec1, &naos,
            vec2, &naos,
            &beta, dm, &naos);
#endif
    }

    template<>
    void matdot(const std::complex<double>* const vec1, const std::complex<double>* const vec2, const Parallel_2D& pvec,
        std::complex<double>* const dm, const Parallel_2D& pmat)
    {
        const int nocc = pvec.get_global_col_size();
        const int naos = pvec.get_global_row_size();
        const std::complex<double> alpha(1.0, 0.0), beta(0.0, 0.0);
        const char transa = 'N', transb = 'C';
#ifdef __MPI
        const int i1 = 1;
        pzgemm_(&transa, &transb, &naos, &naos, &nocc,
            &alpha, vec1, &i1, &i1, pvec.desc,
            vec2, &i1, &i1, pvec.desc,
            &beta, dm, &i1, &i1, pmat.desc);
#else
        zgemm_(&transa, &transb, &naos, &naos, &nocc,
            &alpha, vec1, &naos,
            vec2, &naos,
            &beta, dm, &naos);
#endif
    }
} // namespace LR
