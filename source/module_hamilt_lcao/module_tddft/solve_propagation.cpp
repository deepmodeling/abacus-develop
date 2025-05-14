#include "solve_propagation.h"

#include <iostream>

#include "module_base/lapack_connector.h"
#include "module_base/scalapack_connector.h"

namespace module_tddft
{
#ifdef __MPI
void solve_propagation(const Parallel_Orbitals* pv,
                        const int nband,
                        const int nlocal,
                        const double dt,
                        const std::complex<double>* Stmp,
                        const std::complex<double>* Htmp,
                        const std::complex<double>* psi_k_laststep,
                        std::complex<double>* psi_k)
{
    // (1) init A,B and copy Htmp to A & B
    std::complex<double>* operator_A = new std::complex<double>[pv->nloc];
    ModuleBase::GlobalFunc::ZEROS(operator_A, pv->nloc);
    BlasConnector::copy(pv->nloc, Htmp, 1, operator_A, 1);

    std::complex<double>* operator_B = new std::complex<double>[pv->nloc];
    ModuleBase::GlobalFunc::ZEROS(operator_B, pv->nloc);
    BlasConnector::copy(pv->nloc, Htmp, 1, operator_B, 1);
    
    // ->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // (2) compute operator_A & operator_B by GEADD
    // operator_A = Stmp + i*para * Htmp;   beta2 = para = 0.25 * dt
    // operator_B = Stmp - i*para * Htmp;     beta1 = - para = -0.25 * dt
    std::complex<double> alpha = {1.0, 0.0};
    std::complex<double> beta1 = {0.0, -0.25 * dt};
    std::complex<double> beta2 = {0.0, 0.25 * dt};

    ScalapackConnector::geadd('N',
                              nlocal,
                              nlocal,
                              alpha,
                              Stmp,
                              1,
                              1,
                              pv->desc,
                              beta2,
                              operator_A,
                              1,
                              1,
                              pv->desc);
    ScalapackConnector::geadd('N',
                              nlocal,
                              nlocal,
                              alpha,
                              Stmp,
                              1,
                              1,
                              pv->desc,
                              beta1,
                              operator_B,
                              1,
                              1,
                              pv->desc);
    // ->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // (3) b = operator_B @ psi_k_laststep
    std::complex<double>* tmp_b = new std::complex<double>[pv->nloc_wfc];
    ScalapackConnector::gemm('N',
                        'N',
                        nlocal,
                        nband,
                        nlocal,
                        1.0,
                        operator_B,
                        1,
                        1,
                        pv->desc,
                        psi_k_laststep,
                        1,
                        1,
                        pv->desc_wfc,
                        0.0,
                        tmp_b,
                        1,
                        1,
                        pv->desc_wfc);
    //get ipiv
    int* ipiv = new int[pv->nloc];
    int info = 0;
    // (4) solve Ac=b
    ScalapackConnector::gesv(nlocal,
                            nband,
                            operator_A,
                            1,
                            1,
                            pv->desc,
                            ipiv,
                            tmp_b,
                            1,
                            1,
                            pv->desc_wfc,
                            &info);

    //copy solution to psi_k
    BlasConnector::copy(pv->nloc_wfc, tmp_b, 1, psi_k, 1);

    delete []tmp_b;
    delete []ipiv;
    delete []operator_A;
    delete []operator_B;
}
#endif // __MPI
} // namespace module_tddft
