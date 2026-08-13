#include "cal_edm_from_multipliers.h"
#include "source_base/module_external/scalapack_connector.h"
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
        const double alpha = 1.0, beta = 0.0;
        const int i1 = 1, ivirt = nocc + 1;
        const char transa = 'N', transb = 'N';
        pdgemm_(&transa, &transb, &naos, &nocc, &nvirt,
            &alpha, c, &i1, &ivirt, pc.desc,
            X, &i1, &i1, px.desc,
            &beta, X_ao_occ, &i1, &i1, px_ao_occ.desc);
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
        const int i1 = 1, ivirt = nocc + 1;
        const char transa = 'N', transb = 'N';
        pzgemm_(&transa, &transb, &naos, &nocc, &nvirt,
            &alpha, c, &i1, &ivirt, pc.desc,
            X, &i1, &i1, px.desc,
            &beta, X_ao_occ, &i1, &i1, px_ao_occ.desc);
    }
    // D=X1*X2^T
    // $D_{\mu\nu} = \sum_i X1_{\mu i}X2_{\nu i}$
    template<>
    void matdot(const double* const vec1, const double* const vec2, const Parallel_2D& pvec,
        double* const dm, const Parallel_2D& pmat)
    {
        const int nocc = pvec.get_global_col_size();
        const int naos = pvec.get_global_row_size();
        const double alpha = 1.0, beta = 0.0;
        const int i1 = 1;
        const char transa = 'N', transb = 'T';
        pdgemm_(&transa, &transb, &naos, &naos, &nocc,
            &alpha, vec1, &i1, &i1, pvec.desc,
            vec2, &i1, &i1, pvec.desc,
            &beta, dm, &i1, &i1, pmat.desc);
    }

    template<>
    void matdot(const std::complex<double>* const vec1, const std::complex<double>* const vec2, const Parallel_2D& pvec,
        std::complex<double>* const dm, const Parallel_2D& pmat)
    {
        const int nocc = pvec.get_global_col_size();
        const int naos = pvec.get_global_row_size();
        const std::complex<double> alpha(1.0, 0.0), beta(0.0, 0.0);
        const int i1 = 1;
        const char transa = 'N', transb = 'C';
        pzgemm_(&transa, &transb, &naos, &naos, &nocc,
            &alpha, vec1, &i1, &i1, pvec.desc,
            vec2, &i1, &i1, pvec.desc,
            &beta, dm, &i1, &i1, pmat.desc);
    }
} // namespace LR
