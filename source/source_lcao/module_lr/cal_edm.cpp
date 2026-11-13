#include "cal_edm.h"
#include "source_base/module_external/scalapack_connector.h"
#include "source_base/module_external/blas_connector.h"
#include <type_traits>

namespace LR
{
    // X_ao_occ = C_virtual X; all matrices are column-major.
    template<typename T>
    void cal_X_ao_occ(const T* const X, const BlockCyclicDistribution& px,
        const T* const c, const BlockCyclicDistribution& pc,
        T* const X_ao_occ, const BlockCyclicDistribution& px_ao_occ)
    {
        const int nocc = px.ncglb();
        const int nvirt = px.nrglb();
        const int naos = pc.nrglb();
        const T alpha = T(1.0);
        const T beta = T(0.0);
        const char normal = 'N';
#ifdef __MPI
        const int first = 1;
        const int virtual_start = nocc + 1;
        ScalapackConnector::gemm(normal, normal, naos, nocc, nvirt,
            alpha, c, first, virtual_start, pc.desc(),
            X, first, first, px.desc(), beta, X_ao_occ, first, first, px_ao_occ.desc());
#else
        const T* const c_virtual = c + nocc * naos;
        BlasConnector::gemm_cm(normal, normal, naos, nocc, nvirt,
            alpha, c_virtual, naos, X, nvirt, beta, X_ao_occ, naos);
#endif
    }

    // D = vec1 vec2^dagger (transpose for real coefficients).
    template<typename T>
    void matdot(const T* const vec1, const T* const vec2, const BlockCyclicDistribution& pvec,
        T* const dm, const BlockCyclicDistribution& pmat)
    {
        const int nocc = pvec.ncglb();
        const int naos = pvec.nrglb();
        const T alpha = T(1.0);
        const T beta = T(0.0);
        const char normal = 'N';
        const char adjoint = std::is_same<T, double>::value ? 'T' : 'C';
#ifdef __MPI
        const int first = 1;
        ScalapackConnector::gemm(normal, adjoint, naos, naos, nocc,
            alpha, vec1, first, first, pvec.desc(),
            vec2, first, first, pvec.desc(), beta, dm, first, first, pmat.desc());
#else
        BlasConnector::gemm_cm(normal, adjoint, naos, naos, nocc,
            alpha, vec1, naos, vec2, naos, beta, dm, naos);
#endif
    }

    template void cal_X_ao_occ<double>(const double*, const BlockCyclicDistribution&, const double*,
        const BlockCyclicDistribution&, double*, const BlockCyclicDistribution&);
    template void cal_X_ao_occ<std::complex<double>>(const std::complex<double>*, const BlockCyclicDistribution&,
        const std::complex<double>*, const BlockCyclicDistribution&, std::complex<double>*, const BlockCyclicDistribution&);
    template void matdot<double>(const double*, const double*, const BlockCyclicDistribution&, double*, const BlockCyclicDistribution&);
    template void matdot<std::complex<double>>(const std::complex<double>*, const std::complex<double>*,
        const BlockCyclicDistribution&, std::complex<double>*, const BlockCyclicDistribution&);
}
