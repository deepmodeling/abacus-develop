#include "solve_propagation.h"

#include "source_base/block_cyclic_distribution_information.h"
#include "source_base/constants.h"
#include "source_base/global_function.h"
#include "source_base/module_external/blas_connector.h"
#include "source_base/module_external/scalapack_connector.h"

#include <iostream>

namespace module_rt
{
#ifdef __MPI
void solve_propagation(const BlockCyclicDistributionInformation& matrix_layout,
                        const BlockCyclicDistributionInformation& wfc_layout,
                        const int nband,
                        const int nlocal,
                        const double dt,
                        const std::complex<double>* Stmp,
                        const std::complex<double>* Htmp,
                        const std::complex<double>* psi_k_laststep,
                        std::complex<double>* psi_k)
{
    // (1) init A,B and copy Htmp to A & B
    std::complex<double>* operator_A = new std::complex<double>[matrix_layout.size()];
    ModuleBase::GlobalFunc::ZEROS(operator_A, matrix_layout.size());
    BlasConnector::copy(matrix_layout.size(), Htmp, 1, operator_A, 1);

    std::complex<double>* operator_B = new std::complex<double>[matrix_layout.size()];
    ModuleBase::GlobalFunc::ZEROS(operator_B, matrix_layout.size());
    BlasConnector::copy(matrix_layout.size(), Htmp, 1, operator_B, 1);

    const double dt_au = dt / ModuleBase::AU_to_FS;
    
    // ->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // (2) compute operator_A & operator_B by GEADD
    // operator_A = Stmp + i*para * Htmp;   beta2 = para = 0.25 * dt
    // operator_B = Stmp - i*para * Htmp;     beta1 = - para = -0.25 * dt
    std::complex<double> alpha = {1.0, 0.0};
    std::complex<double> beta1 = {0.0, -0.25 * dt_au};
    std::complex<double> beta2 = {0.0, 0.25 * dt_au};

    ScalapackConnector::geadd('N', nlocal, nlocal, alpha, Stmp, 1, 1, matrix_layout.desc(), beta2, operator_A, 1, 1, matrix_layout.desc());
    ScalapackConnector::geadd('N', nlocal, nlocal, alpha, Stmp, 1, 1, matrix_layout.desc(), beta1, operator_B, 1, 1, matrix_layout.desc());
    // ->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // (3) b = operator_B @ psi_k_laststep
    std::complex<double>* tmp_b = new std::complex<double>[wfc_layout.size()];
    ScalapackConnector::gemm('N',
                        'N',
                        nlocal,
                        nband,
                        nlocal,
                        1.0,
                        operator_B,
                        1,
                        1,
                        matrix_layout.desc(),
                        psi_k_laststep,
                        1,
                        1,
                        wfc_layout.desc(),
                        0.0,
                        tmp_b,
                        1,
                        1,
                        wfc_layout.desc());
    //get ipiv
    int* ipiv = new int[matrix_layout.size()];
    int info = 0;
    // (4) solve Ac=b
    ScalapackConnector::gesv(nlocal, nband, operator_A, 1, 1, matrix_layout.desc(), ipiv, tmp_b, 1, 1, wfc_layout.desc(), &info);

    // copy solution to psi_k
    BlasConnector::copy(wfc_layout.size(), tmp_b, 1, psi_k, 1);

    delete[] tmp_b;
    delete[] ipiv;
    delete[] operator_A;
    delete[] operator_B;
}

void solve_propagation(const BlockCyclicDistributionInformation& matrix_layout,
                       const BlockCyclicDistributionInformation& wfc_layout,
                       const int nband,
                       const int nlocal,
                       const double dt,
                       const std::complex<double>* Stmp,
                       const std::complex<double>* Htmp,
                       const std::complex<double>* P_k, ///< receives P_k
                       const std::complex<double>* psi_k_laststep,
                       std::complex<double>* psi_k)
{
    // Print message for debugging, should be removed later
    std::cout << "Entering solve_propagation with moving gauge P_k..." << std::endl;
    // (1) init A, B and compute HPtmp = Htmp + P_k
    std::complex<double>* operator_A = new std::complex<double>[matrix_layout.size()];
    std::complex<double>* operator_B = new std::complex<double>[matrix_layout.size()];

    // Add up Htmp and P_k to get the effective Hamiltonian matrix for moving spatial gauge
    for (int i = 0; i < matrix_layout.size(); ++i)
    {
        operator_A[i] = Htmp[i] + P_k[i];
        operator_B[i] = Htmp[i] + P_k[i];
    }

    const double dt_au = dt / ModuleBase::AU_to_FS;

    // ->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // (2) compute operator_A & operator_B by GEADD
    // operator_A = Stmp + i*para * (Htmp + P_k);
    // operator_B = Stmp - i*para * (Htmp + P_k);
    std::complex<double> alpha = {1.0, 0.0};
    std::complex<double> beta1 = {0.0, -0.25 * dt_au};
    std::complex<double> beta2 = {0.0, 0.25 * dt_au};

    ScalapackConnector::geadd('N', nlocal, nlocal, alpha, Stmp, 1, 1, matrix_layout.desc(), beta2, operator_A, 1, 1, matrix_layout.desc());
    ScalapackConnector::geadd('N', nlocal, nlocal, alpha, Stmp, 1, 1, matrix_layout.desc(), beta1, operator_B, 1, 1, matrix_layout.desc());

    // ->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // (3) b = operator_B @ psi_k_laststep
    std::complex<double>* tmp_b = new std::complex<double>[wfc_layout.size()];
    ScalapackConnector::gemm('N',
                             'N',
                             nlocal,
                            nband,
                             nlocal,
                             1.0,
                             operator_B,
                            1,
                            1,
                            matrix_layout.desc(),
                             psi_k_laststep,
                             1,
                             1,
                             wfc_layout.desc(),
                             0.0,
                            tmp_b,
                            1,
                            1,
                             wfc_layout.desc());

    // get ipiv
    int* ipiv = new int[matrix_layout.size()];
    int info = 0;

    // (4) solve Ac=b
    ScalapackConnector::gesv(nlocal, nband, operator_A, 1, 1, matrix_layout.desc(), ipiv, tmp_b, 1, 1, wfc_layout.desc(), &info);

    //copy solution to psi_k
    BlasConnector::copy(wfc_layout.size(), tmp_b, 1, psi_k, 1);

    delete []tmp_b;
    delete []ipiv;
    delete []operator_A;
    delete []operator_B;
}
#endif // __MPI
} // namespace module_rt
