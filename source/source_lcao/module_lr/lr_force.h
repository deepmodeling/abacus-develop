#ifndef ABACUS_LR_FORCE_H
#define ABACUS_LR_FORCE_H
#include <memory>
#include "force_funcs.h"
#include "source_lcao/module_lr/potentials/pot_hxc_lrtd.h"
#include "source_lcao/module_lr/potentials/pot_grad_xc.h"
// free functions, usefull for both ground and excited state
#ifdef __EXX
#include "source_lcao/module_ri/exx_lri.h"
#include "exx_force_dm.h"
using TAC = std::pair<int, std::array<int, 3>>;
#endif
namespace LR
{
    /// `dm_gs` carries the ground-state occupations, at nspin=1 they are 2 for fully occupied
    /// bands, so any term that contracts against ONE channel of `dm_gs` needs this factor.
    /// Closed-shell bookkeeping only: the open-shell path sums the two channels explicitly and
    /// must not apply it.
inline double gs_dm_channel_factor(const int& nspin)
{
    return (nspin == 1) ? 0.5 : 1.0;
}

    /// Force-term switches and functional name, owned by the caller's resolved INPUT.
    struct ForceSettings
    {
        int nspin;
        bool test_force;
        bool vl_in_h;
        bool vh_in_h;
        std::string dft_functional;
    };

    template<typename TK>
    class LR_Force
    {
    public:
        LR_Force(const UnitCell& ucell,
            const std::vector<ModuleBase::Vector3<double>>& kvec_d,
            const Parallel_Orbitals& pv,
            const ModulePW::PW_Basis& rhodpw,
            const ModulePW::PW_Basis& rhopw,
            const pseudopot_cell_vl& locpp,
            const Structure_Factor& sf,
            const Grid_Driver& gd,
            const TwoCenterBundle& two_center_bundle,
            const ForceSettings& settings,
            std::ofstream& running_log
#ifdef __EXX
            , std::weak_ptr<Exx_LRI<TK>> exx_lri_in,
            const double& alpha
#endif
        )///<  for 2-center integrals
            : ucell_(ucell), kvec_d_(kvec_d), pv_(pv),
            rhodpw_(rhodpw), rhopw_(rhopw), sf_(sf), locpp_(locpp), gd_(gd),
            two_center_bundle_(two_center_bundle)
#ifdef __EXX
            , exx_lri_(exx_lri_in), alpha_(alpha)
#endif
            , ofs_running_(running_log), nspin_(settings.nspin), test_force_(settings.test_force),
            vl_in_h_(settings.vl_in_h), vh_in_h_(settings.vh_in_h), dft_functional_(settings.dft_functional)
        {
        }

        /// 1. $Tr[H_{GS}^x * (T+D^Z)]$, where GS=groud state and $(T+D^Z)$ is the relaxed difference density matrix
        ModuleBase::matrix cal_force_hamilt_gs_dm_relaxed_diff(const module_dm::DensityMatrix<TK, double>& relaxed_diff_dm,
            const module_dm::DensityMatrix<TK, double>& dm_gs, const bool reproduce_gs = false,
            const PotHxcLR* pot_hxc_gs = nullptr);

        /// 2. $Tr[S^x * (EDM)]
        ModuleBase::matrix cal_force_overlap_edm(const module_dm::DensityMatrix<TK, double>& edm);

        /// 3. $\sum_{mnkl}(mn|f_{Hxc}|kl)^x *D^X *D^X$
        ModuleBase::matrix cal_force_hxc_dmtrans(const module_dm::DensityMatrix<TK, double>& dm_trans, const PotHxcLR& pot_hxc);

        /// 3b. the $g^{xc}$ half of $\partial_x K^{S/T}[D^X]D^X$:
        ///     $\int v^{(2)}[\rho^X,\rho^X](r)\,\partial_x\rho^\text{gs}(r)|_\text{basis}$
        ModuleBase::matrix cal_force_gxc_dmtrans(const module_dm::DensityMatrix<TK, double>& dm_trans,
            const module_dm::DensityMatrix<TK, double>& dm_gs, const PotGradXCLR& pot_grad);

        /// 3b'. open-shell version: $\sum_\tau\int v^{(2)}_\tau[\rho^X,\rho^X]\,
        ///      \partial_x\rho^\text{gs}_\tau|_\text{basis}$. Both transition-density channels
        ///      enter each $v^{(2)}_\tau$, so this cannot be a per-channel loop over the above.
        ModuleBase::matrix cal_force_gxc_dmtrans_openshell(const module_dm::DensityMatrix<TK, double>& dm_trans,
            const module_dm::DensityMatrix<TK, double>& dm_gs, const PotGradXCLR& pot_grad);

#ifdef __EXX
        // auto* lrexx_ptr = dynamic_cast<RI::LR<int, std::array<int, 3>, 3, TK>*>(&exx_lri_in.get());
        /// 4. $\alpha \sum_{mnkl}(mk|nl)^x *D^X *D^X$
        ModuleBase::matrix cal_force_exx_dm_trans(
            const std::map<int, std::map<TAC, RI::Tensor<TK>>>& dm_trans,
            const double& alpha,
            const std::string& spin_suffix = "");
        ModuleBase::matrix cal_force_exx_gs_dm_relaxed_diff(
            const std::map<int, std::map<TAC, RI::Tensor<TK>>>& dm_gs,
            const std::map<int, std::map<TAC, RI::Tensor<TK>>>& relaxed_diff_dm,
            const double& alpha,
            const std::string& spin_suffix = "");
#endif

        // test functions
        /// reproduce the force of the ground state
        ModuleBase::matrix reproduce_force_gs(const K_Vectors& kv,
            const module_dm::DensityMatrix<TK, double>& dm_gs,
            const module_dm::DensityMatrix<TK, double>& edm_gs);

        /// repreduce the ground state local term
        ModuleBase::matrix reproduce_force_gs_loc(const module_dm::DensityMatrix<TK, double>& dm_gs,
            const elecstate::Potential& pot_gs);

    protected:
        const UnitCell& ucell_;
        const std::vector<ModuleBase::Vector3<double>>& kvec_d_;
        const Parallel_Orbitals& pv_;
        const ModulePW::PW_Basis& rhodpw_;
        const ModulePW::PW_Basis& rhopw_;
        const pseudopot_cell_vl& locpp_;
        const Structure_Factor& sf_;
        const Grid_Driver& gd_;
        const TwoCenterBundle& two_center_bundle_;
#ifdef __EXX
        std::weak_ptr<Exx_LRI<TK>> exx_lri_;
        const double alpha_;
#endif
        /// Explicit configuration and log borrowed/passed at the ESolver wiring point.
        std::ofstream& ofs_running_;
        const int nspin_;
        const bool test_force_;
        const bool vl_in_h_;
        const bool vh_in_h_;
        const std::string dft_functional_;

        // PotXC borrows these buffers; they outlive every local potential created here.
        double etxc_ = 0.0;
        double vtxc_ = 0.0;

        void dm_to_charge(const module_dm::DensityMatrix<TK, double>& dm, Charge& chr_out);
        std::unique_ptr<elecstate::Potential> dm_to_hxc_potential(const module_dm::DensityMatrix<TK, double>& dm);
        std::unique_ptr<elecstate::Potential> local_potential();
    };
}

#endif
