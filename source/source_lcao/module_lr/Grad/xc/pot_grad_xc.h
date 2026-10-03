#pragma once
#include "source_lcao/module_lr/potentials/xc_kernel.h"
#include "source_lcao/module_lr/potentials/pot_lr_base.h"

namespace LR
{
    /// the "potential" contributing to RHS of Z-vector equation
    /// from the derivative of xc kernel
    class PotGradXCLR : public PotLRBase
    {
    public:
        /// constructor for exchange-correlation kernel
        /// `triplet` selects which spin combination of $g^{xc}$ to use. It matters: 
        /// the triplet kernel $K^T_{xc}=f_{uu}-f_{ud}$ is NOT zero for a local functional.
        PotGradXCLR(const KernelXC& xc_kernel_in, const ModulePW::PW_Basis& rho_basis, const UnitCell& ucell,
            const int& nrxx, const bool triplet = false);
        ~PotGradXCLR() {}
        virtual void cal_v_eff(double** rho, const UnitCell& ucell, ModuleBase::matrix& v_eff, const std::vector<int>& ispin_op = { 0,0 }) const override;
        /// @brief Open-shell (spin-unrestricted) $v^{(2)}_\tau$.
        /// Unlike the closed-shell case this is NOT a per-(sl,sr) block object: the free spin
        /// $\tau$ is the output index, but BOTH channels of the transition density enter the
        /// same expression, so `rho1` must carry `rho1[0]` and `rho1[1]` together.
        ///     $v^{(2)}_\tau = A_\tau - \nabla\cdot\boldsymbol{E}_\tau$
        void cal_v_eff_openshell(const double* const* const rho1, const UnitCell& ucell,
            ModuleBase::matrix& v_eff, const int tau) const;
        /// kernel components from PotHxcLR
        const KernelXC& xc_kernel_components_;
        const bool triplet_ = false;

    private:
        /// Scratch, shared by every `PotGradXCLR` and grown on demand. These used to be
        /// allocated (and value-initialized) on every call. Safe to share because
        /// `cal_v_eff` is only ever entered from a single thread (all the OpenMP is inside).
        struct Scratch
        {
            std::vector<ModuleBase::Vector3<double>> drho1[2];  ///< $\nabla\rho^1_\sigma$
            std::vector<ModuleBase::Vector3<double>> gdot;      ///< integrand of the divergence
            std::vector<double> vtmp;                           ///< local part $A$
            std::vector<double> div;                            ///< $-\nabla\cdot\boldsymbol{E}$
            void alloc(const int nrxx, const bool two_channel, const bool gga);
        };
        static Scratch& scratch();
    };

}