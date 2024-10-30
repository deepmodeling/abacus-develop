#pragma once
#include "module_elecstate/potentials/pot_base.h"
#include "module_elecstate/potentials/H_Hartree_pw.h"
#include "kernel.h"
#include <unordered_map>
#include <memory>
namespace LR
{
    class PotHxcLR : public elecstate::PotBase
    {
    public:
        /// S1: K^Hartree + K^xc
        /// S2_singlet: 2*K^Hartree + K^xc_{upup} + K^xc_{updown}
        /// S2_triplet: K^xc_{upup} - K^xc_{updown}
        enum SpinType { S1 = 0, S2_singlet = 1, S2_triplet = 2, S2_updown = 3 };
        enum XCType { None = 0, LDA = 1, GGA = 2, HYB_GGA = 4 };
        /// constructor for exchange-correlation kernel
        PotHxcLR(const std::string& xc_kernel_in, const ModulePW::PW_Basis* rho_basis_in, const UnitCell* ucell_in, const Charge* chg_gs/*ground state*/,
            const SpinType& st_in = SpinType::S1);
        ~PotHxcLR() {}
        void cal_v_eff(const Charge* chg/*excited state*/, const UnitCell* ucell, ModuleBase::matrix& v_eff) override {};
        void cal_v_eff(double** rho, const UnitCell* ucell, ModuleBase::matrix& v_eff, const std::vector<int>& ispin_op = { 0,0 });
        int nrxx;
    private:
        int nspin;
        std::unique_ptr<elecstate::PotHartree> pot_hartree;
        /// different components of local and semi-local xc kernels:
        /// LDA: v2rho2
        /// GGA: v2rho2, v2rhosigma, v2sigma2
        /// meta-GGA: v2rho2, v2rhosigma, v2sigma2, v2rholap, v2rhotau, v2sigmalap, v2sigmatau, v2laptau, v2lap2, v2tau2
        KernelXC xc_kernel_components_;
        const std::string xc_kernel;
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
