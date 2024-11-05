#pragma once
#include "module_basis/module_pw/pw_basis.h"

namespace LR
{
    /// @brief Calculate the exchange-correlation (XC) kernel ($f_{xc}=\delta^2E_xc/\delta\rho^2$) and store its components.
    class KernelXC
    {
    public:
        KernelXC(const ModulePW::PW_Basis& rho_basis) : rho_basis_(rho_basis) {};
        ~KernelXC() {};

#ifdef USE_LIBXC
        /// @brief Calculate the XC kernel using libxc.
        void f_xc_libxc(const int& nspin, const double& omega, const double& tpiba, const double* const* const rho_gs, const double* const rho_core = nullptr);
#endif

        // See https://libxc.gitlab.io/manual/libxc-5.1.x/ for the naming convention of the following members.
        // std::map<std::string, std::vector<double>> kernel_set_; // [kernel_type][nrxx][nspin]
        std::vector<double> vrho;
        std::vector<double> vsigma;
        std::vector<double> v2rho2;
        std::vector<double> v2rhosigma;
        std::vector<double> v2sigma2;
        // std::map<std::string, std::vector<ModuleBase::Vector3<double>>> grad_kernel_set_;// [kernel_type][nrxx][nspin],  intermediate terms for GGA
        std::vector<ModuleBase::Vector3<double>> drho_gs;
        std::vector<ModuleBase::Vector3<double>> v2rhosigma_2drho;  ///< $f^{\rho\sigma}*\nabla\rho *2$
        std::vector<ModuleBase::Vector3<double>> v2sigma2_4drho; ///< $f^{\sigma\sigma}*\nabla\rho *4$

    protected:
        const ModulePW::PW_Basis& rho_basis_;
    };
}

