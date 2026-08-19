#pragma once
#include "source_estate/module_pot/h_hartree_pw.h"
#include "xc_kernel.h"
#include "pot_lr_base.h"
#include <unordered_map>
#include <memory>

namespace LR
{
    class PotHxcLR : public PotLRBase
    {
    public:
        /// S1: K^Hartree + K^xc
        /// S2_singlet: 2*K^Hartree + K^xc_{upup} + K^xc_{updown}
        /// S2_triplet: K^xc_{upup} - K^xc_{updown}
        /// S2_updown: K^Hartree + (K^xc_{upup}, K^xc_{updown},  K^xc_{downup} or K^xc_{downdown}), according to `ispin_op` (for spin-polarized systems)
        /// S2_gs: the nspin=2 counterpart of S1, i.e. the *ground-state* Hxc kernel
        ///     K^Hartree + (K^xc_{upup} + K^xc_{updown})/2 = S2_singlet / 2.
        ///     Used for `pot_hxc_gs` in LR gradients, where the convention is
        ///     `K_Hxc(singlet) = 2 * pot_hxc_gs` (see `cal_multiplier_w_from_z.h`).
        ///     The 1/2 on the xc part is not a convention but the chain rule: the derivative is
        ///     taken w.r.t. the *total* density matrix, and $\partial v_u/\partial\rho =
        ///     (f_{uu}+f_{ud})/2$ because $\rho_u=\rho_d=\rho/2$. The Hartree part needs no
        ///     halving, which is exactly why S1 and S2_gs share the same Hartree weight.
        ///     Do NOT use S1 here when nspin=2: `KernelXC` is built with `PARAM.inp.nspin`, so the
        ///     kernel arrays carry 3 spin components per grid point while the S1 integrand indexes
        ///     them as if there were 1.
        enum SpinType { S1 = 0, S2_singlet = 1, S2_triplet = 2, S2_updown = 3, S2_gs = 4 };
        /// XCType here is to determin the method of integration from kernel to potential, not the way calculating the kernel
        enum XCType { None = 0, LDA = 1, GGA = 2, HYB_GGA = 4 };
        /// constructor for exchange-correlation kernel
        PotHxcLR(const std::string& xc_kernel, const ModulePW::PW_Basis& rho_basis,
            const UnitCell& ucell, const Charge& chg_gs/*ground state*/, const Parallel_Grid& pgrid,
            const SpinType& st = SpinType::S1, const std::vector<std::string>& lr_init_xc_kernel = { "default" });
        ~PotHxcLR() {}
        virtual void cal_v_eff(double** rho, const UnitCell& ucell, ModuleBase::matrix& v_eff, const std::vector<int>& ispin_op = { 0,0 })  const override;

        // const references
        const KernelXC& xc_kernel_components = xc_kernel_components_;
    private:
        std::unique_ptr<elecstate::PotHartree> pot_hartree_;
        /// different components of local and semi-local xc kernels:
        /// LDA: v2rho2
        /// GGA: v2rho2, v2rhosigma, v2sigma2
        /// meta-GGA: v2rho2, v2rhosigma, v2sigma2, v2rholap, v2rhotau, v2sigmalap, v2sigmatau, v2laptau, v2lap2, v2tau2
        const KernelXC xc_kernel_components_;
        const std::string xc_kernel_;
        const SpinType spin_type_ = SpinType::S1;
        XCType xc_type_ = XCType::None;

        // enum class as key for unordered_map is not supported in C++11 sometimes
        // https://github.com/llvm/llvm-project/issues/49601
        // struct SpinHash { std::size_t operator()(const SpinType& s) const { return std::hash<int>()(static_cast<int>(s)); } };
        using Tfunc = std::function<void(const double* const /**<[in] rho*/,
            ModuleBase::matrix& /**<[out] v_eff */,
            const std::vector<int>& ispin_op) >;
        // std::unordered_map<SpinType, Tfunc, SpinHash> kernel_to_potential_;
        std::map<SpinType, Tfunc> kernel_to_potential_;

        void set_integral_func(const SpinType& s, const XCType& xc);
    };

} // namespace LR
