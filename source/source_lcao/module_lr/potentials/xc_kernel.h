#pragma once
#include "source_basis/module_pw/pw_basis.h"
#include "source_cell/unitcell.h"
#include "source_base/parallel_grid.h"
#include "source_estate/module_charge/charge.h"
#define CREF(x) const std::vector<double>& x = x##_
#define CREF3(x) const std::vector<ModuleBase::Vector3<double>>& x = x##_
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
        ~KernelXC() {}

        // const references
        CREF(vrho);CREF(vsigma); CREF(v2rho2); CREF(v2rhosigma); CREF(v2sigma2);
        CREF3(v2rhosigma_2drho); CREF3(v2sigma2_4drho);
        CREF3(v2rhosigma_drho_singlet); CREF3(v2rhosigma_drho_triplet); CREF3(v2sigma2_drho_singlet);CREF3(v2sigma2_drho_triplet);
        CREF3(v2rhosigma_drho_uu); CREF3(v2rhosigma_drho_ud); CREF3(v2rhosigma_drho_du); CREF3(v2rhosigma_drho_dd);
        CREF3(v2sigma2_drho_uu_u); CREF3(v2sigma2_drho_uu_d); CREF3(v2sigma2_drho_ud_u); CREF3(v2sigma2_drho_ud_d);
        CREF3(v2sigma2_drho_du_u); CREF3(v2sigma2_drho_du_d); CREF3(v2sigma2_drho_dd_u); CREF3(v2sigma2_drho_dd_d);
        CREF(v3rho3); CREF(v3rho2sigma); CREF(v3rhosigma2); CREF(v3sigma3);
        /// @brief The coefficients of
        ///     $v^{(2)}(r)=\iint dr'dr''\,g^{xc}(r,r',r'')\rho^1(r')\rho^1(r'')$
        /// for ONE spin combination, stored exactly as they appear in the final formula:
        ///     $v^{(2)} = a_{s^2}s^2 + a_{st}\,s\,t + a_{t^2}t^2 + a_q\,q
        ///        - \nabla\cdot[\,\boldsymbol{e}_{s^2}s^2 + \boldsymbol{e}_{st}\,s\,t
        ///          + \boldsymbol{e}_{t^2}t^2 + \boldsymbol{e}_q\,q
        ///          + (c_s\,s + c_t\,t)\,\nabla\rho^1\,]$
        /// with $s=\rho^1$, $t=\nabla\rho\cdot\nabla\rho^1$, $q=\nabla\rho^1\cdot\nabla\rho^1$.
        ///
        /// Every numeric factor and every spin sum is folded in here, so `PotGradXCLR::cal_v_eff`
        /// is a literal transcription of the formula with no arithmetic of its own, and nspin=1,
        /// singlet and triplet all run through the same code. For LDA only `a_s2` is filled.
        /// $v^{(2)} = A - \nabla\cdot E$ 
        struct GxcCoef
        {
            std::vector<double> a_s2, a_st, a_t2, a_q;  ///< the local part $A$
            std::vector<double> c_s, c_t;               ///< the two $\nabla\rho^1$-weighted scalars
            std::vector<ModuleBase::Vector3<double>> e_s2, e_st, e_t2, e_q;  ///< under the divergence
        };
        /// nspin=1 has no singlet/triplet distinction, so it always returns the one set that is built.
        const GxcCoef& gxc(const bool triplet) const { return (nspin_ == 1 || !triplet) ? gxc_s_ : gxc_t_; }

        CREF3(v3rho2sigma_2drho); CREF3(v3rhosigma2_8drho); CREF3(v3sigma3_8drho);
        const bool& openshell = openshell_;
        const std::vector<std::vector<ModuleBase::Vector3<double>>>& drho_gs = drho_gs_;
    private:
#ifdef __LIBXC
        /// @brief Calculate the XC kernel using libxc.
        void f_xc_libxc(const int& nspin, const double& omega, const double& tpiba, const double* const* const rho_gs, const double* const rho_core = nullptr);
        /// calculate the input rho, grad rho, and sigma for libxc 
        void get_rho_drho_sigma(const int& nspin,
            const double& tpiba,
            const double* const* const rho_gs,
            const double* const rho_core,
            const bool& is_gga,
            std::vector<double>& rho,
            std::vector<std::vector<ModuleBase::Vector3<double>>>& gradrho,
            std::vector<double>& sigma);
#endif
        // See https://libxc.gitlab.io/manual/libxc-5.1.x/ for the naming convention of the following members.
        // std::map<std::string, std::vector<double>> kernel_set_; // [kernel_type][nrxx][nspin]

        // ================================== XC kernels ============================================
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
        // ================================== XC kernels ============================================
        // ================================== XC kernel Gradiants ====================================
        Tvec v3rho3_;
        Tvec v3rho2sigma_;
        Tvec v3rhosigma2_;
        Tvec v3sigma3_;
        // for nspin=1, gga: the third-order kernels already contracted with $\nabla\rho$, as they
        // appear in the divergence term of $v^{(2)}$ (LR-Grad-formulas/GGA-kxc-to-v积分公式.md, sec. 3).
        // The numeric coefficient is baked into the name so the potential code reads off the formula.
        // The fourth vector the formula needs, $4f^{\sigma\sigma}\nabla\rho$, is `v2sigma2_4drho_`
        // above -- it is shared with the $f\to v$ path.
        Tvec3 v3rho2sigma_2drho_;   ///< $2g^{\rho\rho\sigma}\nabla\rho$
        Tvec3 v3rhosigma2_8drho_;   ///< $8g^{\rho\sigma\sigma}\nabla\rho$
        Tvec3 v3sigma3_8drho_;      ///< $8g^{\sigma\sigma\sigma}\nabla\rho$

        // The two spin combinations of $v^{(2)}$'s coefficients (see `GxcCoef` above).
        // `gxc_t_` stays empty for nspin=1, where there is no triplet.
        GxcCoef gxc_s_;
        GxcCoef gxc_t_;
        /// @brief Fill `dst` for one spin combination. All the spin algebralives here, driven by the weight
        /// vectors that distinguish singlet from triplet -- the two differ only in those weights.
        void build_gxc_coef(GxcCoef& dst, const bool triplet, const int& nspin, const bool& is_gga,
            const std::vector<ModuleBase::Vector3<double>>& drho);
        int nspin_ = 1;
        // ================================== XC kernel Gradiants ====================================
        const ModulePW::PW_Basis& rho_basis_;
        const bool openshell_ = false;
    };
}

