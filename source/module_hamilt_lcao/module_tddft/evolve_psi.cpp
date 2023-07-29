#include "evolve_psi.h"

#include <complex>

#include "bandenergy.h"
#include "middle_hamilt.h"
#include "module_base/lapack_connector.h"
#include "module_base/scalapack_connector.h"
#include "module_hamilt_lcao/hamilt_lcaodft/hamilt_lcao.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_io/input.h"
#include "norm_psi.h"
#include "propagator.h"
#include "upsi.h"

namespace module_tddft
{
void evolve_psi(const int nband,
                const int nlocal,
                const Parallel_Orbitals* pv,
                hamilt::Hamilt<double>* p_hamilt,
                std::complex<double>* psi_k,
                std::complex<double>* psi_k_laststep,
                std::complex<double>* H_laststep,
                std::complex<double>* S_laststep,
                double* ekb,
                int htype,
                int propagator)
{
    ModuleBase::TITLE("Evolve_psi", "evolve_psi");
    time_t time_start = time(NULL);
    GlobalV::ofs_running << " Start Time : " << ctime(&time_start);

#ifdef __MPI

    int print_matrix = 0;
    hamilt::MatrixBlock<std::complex<double>> h_mat, s_mat;
    p_hamilt->matrix(h_mat, s_mat);

    std::complex<double>* Stmp = new std::complex<double>[pv->nloc];
    ModuleBase::GlobalFunc::ZEROS(Stmp, pv->nloc);
    BlasConnector::copy(pv->nloc, s_mat.p, 1, Stmp, 1);

    std::complex<double>* Htmp = new std::complex<double>[pv->nloc];
    ModuleBase::GlobalFunc::ZEROS(Htmp, pv->nloc);
    BlasConnector::copy(pv->nloc, h_mat.p, 1, Htmp, 1);

    std::complex<double>* Hold = new std::complex<double>[pv->nloc];
    ModuleBase::GlobalFunc::ZEROS(Hold, pv->nloc);
    BlasConnector::copy(pv->nloc, h_mat.p, 1, Hold, 1);

    std::complex<double>* U_operator = new std::complex<double>[pv->nloc];
    ModuleBase::GlobalFunc::ZEROS(U_operator, pv->nloc);

    // (1)->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    /// @brief compute H(t+dt/2)
    /// @input H_laststep, Htmp, print_matrix
    /// @output Htmp
    if (htype == 1 && propagator != 2)
    {
        half_Hmatrix(pv, nband, nlocal, Htmp, Stmp, H_laststep, S_laststep, print_matrix);
    }

    // (2)->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    /// @brief compute U_operator
    /// @input Stmp, Htmp, print_matrix
    /// @output U_operator
    Propagator prop(propagator, pv);
    prop.compute_propagator(nlocal, Stmp, Htmp, H_laststep, U_operator, print_matrix);

    // (3)->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    /// @brief apply U_operator to the wave function of the previous step for new wave function
    /// @input U_operator, psi_k_laststep, print_matrix
    /// @output psi_k
    upsi(pv, nband, nlocal, U_operator, psi_k_laststep, psi_k, print_matrix);

    // (4)->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    /// @brief normalize psi_k
    /// @input Stmp, psi_not_norm, psi_k, print_matrix
    /// @output psi_k
    norm_psi(pv, nband, nlocal, Stmp, psi_k, print_matrix);

    // (5)->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    /// @brief compute ekb
    /// @input Htmp, psi_k
    /// @output ekb
    compute_ekb(pv, nband, nlocal, Hold, psi_k, ekb);

    delete[] Stmp;
    delete[] Htmp;
    delete[] Hold;
    delete[] U_operator;

#endif

    time_t time_end = time(NULL);
    ModuleBase::GlobalFunc::OUT_TIME("evolve(std::complex)", time_start, time_end);

    return;
}
} // namespace module_tddft