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

        // Setters and getters for certain components, such as v2rho2, v2rhosigma, v2sigma2,  or other custom combinations.
        void set_kernel(const std::string& name, const std::vector<double>& vec) { this->kernel_set_[name] = vec; }
        void set_kernel(const std::string& name, const std::vector<double>&& vec) { this->kernel_set_[name] = std::move(vec); }
        const std::vector<double>& get_kernel(const std::string& name) { return kernel_set_[name]; }
        const std::vector<ModuleBase::Vector3<double>>& get_grad_kernel(const std::string& name) { return grad_kernel_set_[name]; }

    protected:

        const ModulePW::PW_Basis& rho_basis_;
        std::map<std::string, std::vector<double>> kernel_set_; // [kernel_type][nrxx][nspin]
        std::map<std::string, std::vector<ModuleBase::Vector3<double>>> grad_kernel_set_;// [kernel_type][nrxx][nspin],  intermediate terms for GGA
    };
}

