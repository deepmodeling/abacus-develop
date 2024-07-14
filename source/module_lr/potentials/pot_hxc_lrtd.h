#pragma once
#include "module_elecstate/potentials/pot_base.h"
#include "module_elecstate/potentials/H_Hartree_pw.h"
#include "kernel.h"
#include <unordered_map>
namespace LR
{
    class PotHxcLR : public elecstate::PotBase
    {
        /// S1: K^Hartree + K^xc
        /// S2_singlet: 2*K^Hartree + K^xc_{upup} + K^xc_{updown}
        /// S2_triplet: K^xc_{upup} - K^xc_{updown}
        enum SpinType { S1 = 0, S2_singlet = 1, S2_triplet = 2 } spin_type_ = SpinType::S1;
        enum XCType { None = 0, LDA = 1, GGA = 2, HYB_GGA = 4 } xc_type_ = XCType::None;
    public:
        /// constructor for exchange-correlation kernel
        PotHxcLR(const std::string& xc_kernel_in, const ModulePW::PW_Basis* rho_basis_in, const UnitCell* ucell_in, const Charge* chg_gs/*ground state*/,
            SpinType st_in = S1);
        ~PotHxcLR() {}
        void cal_v_eff(const Charge* chg/*excited state*/, const UnitCell* ucell, ModuleBase::matrix& v_eff) override {};
        void cal_v_eff(double** rho, const UnitCell* ucell, ModuleBase::matrix& v_eff);
        int nrxx;
        int nspin;
        elecstate::PotHartree* pot_hartree;
        /// different components of local and semi-local xc kernels:
        /// LDA: v2rho2
        /// GGA: v2rho2, v2rhosigma, v2sigma2
        /// meta-GGA: v2rho2, v2rhosigma, v2sigma2, v2rholap, v2rhotau, v2sigmalap, v2sigmatau, v2laptau, v2lap2, v2tau2
        KernelXC xc_kernel_components_;
        const std::string xc_kernel;
        const double& tpiba_;

        using Tfunc = std::function<void(const double* const /**<[in] rho*/,
            ModuleBase::matrix& /**<[out] v_eff */)>;
        std::unordered_map<SpinType, Tfunc> kernel_to_potential_;

        void set_integral_func(SpinType s, XCType xc);
    };

} // namespace LR
