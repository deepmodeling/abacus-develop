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
        using Tvec = std::vector<double>;
        using Tvec3 = std::vector<ModuleBase::Vector3<double>>;
    public:
        KernelXC(const ModulePW::PW_Basis& rho_basis,
            const UnitCell& ucell,
            const Charge& chg_gs,
            const Parallel_Grid& pgrid,
            const int& nspin,
            const std::string& kernel_name,
            const std::vector<std::string>& lr_init_xc_kernel,
            const bool openshell = false);
        ~KernelXC() {};

        // const references
        CREF(vrho);CREF(vsigma); CREF(v2rho2); CREF(v2rhosigma); CREF(v2sigma2);
        CREF3(v2rhosigma_2drho); CREF3(v2sigma2_4drho);
        CREF3(v2rhosigma_drho_singlet); CREF3(v2rhosigma_drho_triplet); CREF3(v2sigma2_drho_singlet);CREF3(v2sigma2_drho_triplet);
        CREF3(v2rhosigma_drho_uu); CREF3(v2rhosigma_drho_ud); CREF3(v2rhosigma_drho_du); CREF3(v2rhosigma_drho_dd);
        CREF3(v2sigma2_drho_uu_u); CREF3(v2sigma2_drho_uu_d); CREF3(v2sigma2_drho_ud_u); CREF3(v2sigma2_drho_ud_d);
        CREF3(v2sigma2_drho_du_u); CREF3(v2sigma2_drho_du_d); CREF3(v2sigma2_drho_dd_u); CREF3(v2sigma2_drho_dd_d);
        const std::vector<std::vector<ModuleBase::Vector3<double>>>& drho_gs = drho_gs_;
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
        // for nspin=1, gga
        std::vector<std::vector<ModuleBase::Vector3<double>>> drho_gs_; ///< $\nabla\rho$ of the ground state, size: nspin*nrxx
        std::vector<ModuleBase::Vector3<double>> v2rhosigma_2drho_;  ///< $f^{\rho\sigma}*\nabla\rho *2$
        std::vector<ModuleBase::Vector3<double>> v2sigma2_4drho_; ///< $f^{\sigma\sigma}*\nabla\rho *4$

        // for nspin=2 (close-shell), gga kernels
        Tvec3 v2rhosigma_drho_singlet_; ///< $2(f^{\rho_u\sigma_{uu}}+f^{\rho_u\sigma_{ud}}+f^{\rho_u\sigma_{dd}})*\nabla\rho)$
        Tvec3 v2rhosigma_drho_triplet_; ///< $2(f^{\rho_u\sigma_{uu}}-f^{\rho_u\sigma_{ud}})*\nabla\rho)$
        Tvec3 v2sigma2_drho_singlet_;   /// < $(4f^{\sigma_{uu}\sigma_{uu}}+6f^{\sigma_{uu}\sigma_{ud}}+4f^{\sigma_{uu}\sigma_{dd}}+2f^{\sigma_{ud}\sigma_{ud}}+2f^{\sigma_{ud}\sigma_{dd}})\nabla\rho$
        Tvec3 v2sigma2_drho_triplet_;   /// < $(4f^{\sigma_{uu}\sigma_{uu}}+2f^{\sigma_{uu}\sigma_{ud}}-4f^{\sigma_{uu}\sigma_{dd}}-2f^{\sigma_{ud}\sigma_{dd}})\nabla\rho$

        // for nspin=2 (open shell), gga kernels
        Tvec3 v2rhosigma_drho_uu_;  ///< $2f^{\rho_u\sigma_{uu}}\nabla\rho_u+f^{\rho_u\sigma_{ud}}\nabla\rho_d$
        Tvec3 v2rhosigma_drho_ud_;  ///< $f^{\rho_u\sigma_{ud}}\nabla\rho_u+2f^{\rho_u\sigma_{dd}}\nabla\rho_d$
        Tvec3 v2rhosigma_drho_du_;  ///< $2f^{\rho_d\sigma_{uu}}\nabla\rho_u+f^{\rho_d\sigma_{ud}}\nabla\rho_d$
        Tvec3 v2rhosigma_drho_dd_;  ///< $f^{\rho_d\sigma_{ud}}\nabla\rho_u+2f^{\rho_d\sigma_{dd}}\nabla\rho_d$
        Tvec3 v2sigma2_drho_uu_u_;   /// $4f^{\sigma_{uu}\sigma_{uu}}\nabla\rho_u+2f^{\sigma_{uu}\sigma_{ud}\nabla\rho_d$
        Tvec3 v2sigma2_drho_uu_d_;   /// $2f^{\sigma_{uu}\sigma_{ud}}\nabla\rho_u+f^{\sigma_{ud}\sigma_{ud}}\nabla\rho_d$
        Tvec3 v2sigma2_drho_ud_u_;   /// $2f^{\sigma_{uu}\sigma_{ud}}\nabla\rho_u+f^{\sigma_{ud}\sigma_{ud}}\nabla\rho_d$
        Tvec3 v2sigma2_drho_ud_d_;   /// $4f^{\sigma_{uu}\sigma_{dd}}\nabla\rho_u+2f^{\sigma_{ud}\sigma_{dd}\nabla\rho_d$
        Tvec3 v2sigma2_drho_du_u_;   /// $4f^{\sigma_{uu}\sigma_{dd}}\nabla\rho_d+2f^{\sigma_{uu}\sigma_{ud}\nabla\rho_u$
        Tvec3 v2sigma2_drho_du_d_;   /// $2f^{\sigma_{ud}\sigma_{dd}}\nabla\rho_d+f^{\sigma_{ud}\sigma_{ud}}\nabla\rho_u$
        Tvec3 v2sigma2_drho_dd_u_;   /// $2f^{\sigma_{ud}\sigma_{dd}}\nabla\rho_d+f^{\sigma_{ud}\sigma_{ud}}\nabla\rho_u$
        Tvec3 v2sigma2_drho_dd_d_;   /// $4f^{\sigma_{dd}\sigma_{dd}}\nabla\rho_d+2f^{\sigma_{ud}\sigma_{dd}\nabla\rho_u$
        const ModulePW::PW_Basis& rho_basis_;
        const bool openshell_ = false;
    };
}

