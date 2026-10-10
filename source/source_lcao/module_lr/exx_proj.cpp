#include "exx_proj.h"

#include "source_base/module_external/blas_connector.h"

namespace LR
{
void project_exx(const double* h,
                 const double* left,
                 const double* right,
                 const int naos,
                 const int nleft,
                 const int nright,
                 const double factor,
                 double* scratch,
                 double* result)
{
    if (nleft == 0 || nright == 0) { return; }
    const double one = 1.0;
    const double zero = 0.0;
    // Transpose, rather than symmetrize: K[D^X] need not be symmetric.
    BlasConnector::gemm_cm('T', 'N', naos, nleft, naos, one,
                           h, naos, left, naos, zero, scratch, naos);
    BlasConnector::gemm_cm('T', 'N', nright, nleft, naos, factor,
                           right, naos, scratch, naos, one, result, nright);
}

void project_exx(const std::complex<double>* h,
                 const std::complex<double>* left,
                 const std::complex<double>* right,
                 const int naos,
                 const int nleft,
                 const int nright,
                 const double factor,
                 std::complex<double>* scratch,
                 std::complex<double>* result)
{
    if (nleft == 0 || nright == 0) { return; }
    const std::complex<double> one(1.0, 0.0);
    const std::complex<double> zero(0.0, 0.0);
    const std::complex<double> scale(factor, 0.0);
    BlasConnector::gemm_cm('N', 'N', naos, nleft, naos, one,
                           h, naos, left, naos, zero, scratch, naos);
    BlasConnector::gemm_cm('C', 'N', nright, nleft, naos, scale,
                           right, naos, scratch, naos, one, result, nright);
}
}
