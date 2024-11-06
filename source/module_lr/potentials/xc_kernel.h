#pragma once
#include "module_basis/module_pw/pw_basis.h"
#include "module_cell/unitcell.h"
#include "module_hamilt_pw/hamilt_pwdft/parallel_grid.h"
#include "module_elecstate/module_charge/charge.h"
#define CREF(x) const std::vector<double>& x = x##_;
#define CREF3(x) const std::vector<ModuleBase::Vector3<double>>& x = x##_;
namespace LR
{
    /// @brief Calculate the exchange-correlation (XC) kernel ($f_{xc}=\delta^2E_xc/\delta\rho^2$) and store its components.
    class KernelXC
    {
    public:
        KernelXC(const ModulePW::PW_Basis& rho_basis,
            const UnitCell& ucell,
            const Charge& chg_gs,
            const Parallel_Grid& pgrid,
            const int& nspin,
            const std::string& kernel_name,
            const std::vector<std::string>& lr_init_xc_kernel);
        ~KernelXC() {};

        // const references
        CREF(vrho);CREF(vsigma); CREF(v2rho2); CREF(v2rhosigma); CREF(v2sigma2);
        CREF3(drho_gs); CREF3(v2rhosigma_2drho); CREF3(v2sigma2_4drho);

    private:
#ifdef USE_LIBXC
        /// @brief Calculate the XC kernel using libxc.
        void f_xc_libxc(const int& nspin, const double& omega, const double& tpiba, const double* const* const rho_gs, const double* const rho_core = nullptr);
#endif
        // See https://libxc.gitlab.io/manual/libxc-5.1.x/ for the naming convention of the following members.
        // std::map<std::string, std::vector<double>> kernel_set_; // [kernel_type][nrxx][nspin]
        std::vector<double> vrho_;
        std::vector<double> vsigma_;
        std::vector<double> v2rho2_;
        std::vector<double> v2rhosigma_;
        std::vector<double> v2sigma2_;
        // std::map<std::string, std::vector<ModuleBase::Vector3<double>>> grad_kernel_set_;// [kernel_type][nrxx][nspin],  intermediate terms for GGA
        std::vector<ModuleBase::Vector3<double>> drho_gs_;
        std::vector<ModuleBase::Vector3<double>> v2rhosigma_2drho_;  ///< $f^{\rho\sigma}*\nabla\rho *2$
        std::vector<ModuleBase::Vector3<double>> v2sigma2_4drho_; ///< $f^{\sigma\sigma}*\nabla\rho *4$

        const ModulePW::PW_Basis& rho_basis_;
    };
}

