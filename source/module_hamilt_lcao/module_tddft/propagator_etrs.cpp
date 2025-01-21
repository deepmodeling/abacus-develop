#include "module_base/lapack_connector.h"
#include "module_base/module_container/ATen/kernels/blas.h"
#include "module_base/module_container/ATen/kernels/lapack.h"
#include "module_base/module_container/ATen/kernels/memory.h" // memory operations (Tensor)
#include "module_base/module_device/memory_op.h"              // memory operations
#include "module_base/scalapack_connector.h"
#include "module_parameter/parameter.h"
#include "propagator.h"

#include <complex>
#include <iostream>

namespace module_tddft
{
#ifdef __MPI
void Propagator::compute_propagator_etrs(const int nlocal,
                                         const std::complex<double>* Stmp,
                                         const std::complex<double>* Htmp,
                                         const std::complex<double>* H_laststep,
                                         std::complex<double>* U_operator,
                                         std::ofstream& ofs_running,
                                         const int print_matrix) const
{
    std::vector<std::complex<double>> U1(this->ParaV->nloc);
    std::vector<std::complex<double>> U2(this->ParaV->nloc);
    int tag = 2;
    compute_propagator_taylor(nlocal, Stmp, Htmp, U1.data(), ofs_running, print_matrix, tag);
    compute_propagator_taylor(nlocal, Stmp, H_laststep, U2.data(), ofs_running, print_matrix, tag);
    ScalapackConnector::gemm('N',
                             'N',
                             nlocal,
                             nlocal,
                             nlocal,
                             1.0,
                             U1.data(),
                             1,
                             1,
                             this->ParaV->desc,
                             U2.data(),
                             1,
                             1,
                             this->ParaV->desc,
                             0.0,
                             U_operator,
                             1,
                             1,
                             this->ParaV->desc);
}
#endif // __MPI
} // namespace module_tddft
