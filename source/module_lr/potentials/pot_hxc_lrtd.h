#pragma once
#include "module_elecstate/potentials/H_Hartree_pw.h"
#include "xc_kernel.h"
#include <unordered_map>
#include <memory>

namespace LR
{
    class PotHxcLR
    {
    public:
        /// S1: K^Hartree + K^xc
        /// S2_singlet: 2*K^Hartree + K^xc_{upup} + K^xc_{updown}
        /// S2_triplet: K^xc_{upup} - K^xc_{updown}
        /// S2_updown: K^Hartree + (K^xc_{upup}, K^xc_{updown},  K^xc_{downup} or K^xc_{downdown}), according to `ispin_op` (for spin-polarized systems)
        enum SpinType { S1 = 0, S2_singlet = 1, S2_triplet = 2, S2_updown = 3 };
        /// XCType here is to determin the method of integration from kernel to potential, not the way calculating the kernel
        enum XCType { None = 0, LDA = 1, GGA = 2, HYB_GGA = 4 };
        /// constructor for exchange-correlation kernel
        PotHxcLR(const std::string& xc_kernel, const ModulePW::PW_Basis& rho_basis,
            const UnitCell& ucell, const Charge& chg_gs/*ground state*/, const Parallel_Grid& pgrid,
            const SpinType& st = SpinType::S1, const std::vector<std::string>& lr_init_xc_kernel = { "default" });
        ~PotHxcLR() {}
        void cal_v_eff(double** rho, const UnitCell& ucell, ModuleBase::matrix& v_eff, const std::vector<int>& ispin_op = { 0,0 });
        const int& nrxx = nrxx_;
    private:
        const ModulePW::PW_Basis& rho_basis_;
        const int nspin_ = 1;
        const int nrxx_ = 1;
        std::unique_ptr<elecstate::PotHartree> pot_hartree_;
        /// different components of local and semi-local xc kernels:
        /// LDA: v2rho2
        /// GGA: v2rho2, v2rhosigma, v2sigma2
        /// meta-GGA: v2rho2, v2rhosigma, v2sigma2, v2rholap, v2rhotau, v2sigmalap, v2sigmatau, v2laptau, v2lap2, v2tau2
        const KernelXC xc_kernel_components_;
        const std::string xc_kernel_;
        const double& tpiba_;
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
