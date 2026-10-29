#include "cal_edm.h"
#include "source_base/module_external/scalapack_connector.h"
#include "source_base/module_external/blas_connector.h"
#include <type_traits>

namespace LR
{
    // X_ao_occ = C_virtual X; all matrices are column-major.
    template<typename T>
    void cal_X_ao_occ(const T* const X, const BlockCyclicDistributionInformation& px,
        const T* const c, const BlockCyclicDistributionInformation& pc,
        T* const X_ao_occ, const BlockCyclicDistributionInformation& px_ao_occ)
    {
        const int nocc = px.get_global_col_size();
        const int nvirt = px.get_global_row_size();
        const int naos = pc.get_global_row_size();
        const T alpha = T(1.0);
        const T beta = T(0.0);
        const char normal = 'N';
#ifdef __MPI
        const int first = 1;
        const int virtual_start = nocc + 1;
        ScalapackConnector::gemm(normal, normal, naos, nocc, nvirt,
            alpha, c, first, virtual_start, pc.get_desc(),
            X, first, first, px.get_desc(), beta, X_ao_occ, first, first, px_ao_occ.get_desc());
#else
        const T* const c_virtual = c + nocc * naos;
        BlasConnector::gemm_cm(normal, normal, naos, nocc, nvirt,
            alpha, c_virtual, naos, X, nvirt, beta, X_ao_occ, naos);
#endif
    }

    // D = vec1 vec2^dagger (transpose for real coefficients).
    template<typename T>
    void matdot(const T* const vec1, const T* const vec2, const BlockCyclicDistributionInformation& pvec,
        T* const dm, const BlockCyclicDistributionInformation& pmat)
    {
        const int nocc = pvec.get_global_col_size();
        const int naos = pvec.get_global_row_size();
        const T alpha = T(1.0);
        const T beta = T(0.0);
        const char normal = 'N';
        const char adjoint = std::is_same<T, double>::value ? 'T' : 'C';
#ifdef __MPI
        const int first = 1;
        ScalapackConnector::gemm(normal, adjoint, naos, naos, nocc,
            alpha, vec1, first, first, pvec.get_desc(),
            vec2, first, first, pvec.get_desc(), beta, dm, first, first, pmat.get_desc());
#else
        BlasConnector::gemm_cm(normal, adjoint, naos, naos, nocc,
            alpha, vec1, naos, vec2, naos, beta, dm, naos);
#endif
    }

    template void cal_X_ao_occ<double>(const double*, const BlockCyclicDistributionInformation&, const double*,
        const BlockCyclicDistributionInformation&, double*, const BlockCyclicDistributionInformation&);
    template void cal_X_ao_occ<std::complex<double>>(const std::complex<double>*, const BlockCyclicDistributionInformation&,
        const std::complex<double>*, const BlockCyclicDistributionInformation&, std::complex<double>*, const BlockCyclicDistributionInformation&);
    template void matdot<double>(const double*, const double*, const BlockCyclicDistributionInformation&, double*, const BlockCyclicDistributionInformation&);
    template void matdot<std::complex<double>>(const std::complex<double>*, const std::complex<double>*,
        const BlockCyclicDistributionInformation&, std::complex<double>*, const BlockCyclicDistributionInformation&);
}
