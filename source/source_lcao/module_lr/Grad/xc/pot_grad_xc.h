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
        /// kernel components from PotHxcLR
        const KernelXC& xc_kernel_components_;
        const bool triplet_ = false;

    private:
        /// Scratch, shared by every `PotGradXCLR` and grown on demand. These used to be
        /// allocated (and value-initialized) on every call. Safe to share because
        /// `cal_v_eff` is only ever entered from a single thread (all the OpenMP is inside).
        struct Scratch
        {
            std::vector<ModuleBase::Vector3<double>> drho1;  ///< $\nabla\rho^1$
            std::vector<ModuleBase::Vector3<double>> gdot;   ///< integrand of the divergence
            std::vector<double> vtmp;                        ///< local part $A$
            void alloc(const int nrxx, const bool gga);
        };
        static Scratch& scratch();
    };

}