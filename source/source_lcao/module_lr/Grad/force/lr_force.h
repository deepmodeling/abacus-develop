#include "force_funcs_lcao.h"
#include "source_lcao/module_lr/potentials/pot_hxc_lrtd.h"
#include "source_lcao/module_lr/Grad/xc/pot_grad_xc.h"
// free functions, usefull for both ground and excited state
#ifdef __EXX
#include "source_lcao/module_ri/exx_lri.h"
#include "exx_force_two_dm.h"
using TAC = std::pair<int, std::array<int, 3>>;
#endif
namespace LR
{
    /// `dm_gs` carries the ground-state occupations, at nspin=1 they are 2 for fully occupied
    /// bands, so any term that contracts against ONE channel of `dm_gs` needs this factor.
    /// Closed-shell bookkeeping only: the open-shell path sums the two channels explicitly and
    /// must not apply it.
    inline double gs_dm_channel_factor() { return (PARAM.inp.nspin == 1) ? 0.5 : 1.0; }

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
            const TwoCenterBundle& two_center_bundle
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
        {
        }

        /// 1. $Tr[H_{GS}^x * (T+D^Z)]$, where GS=groud state and $(T+D^Z)$ is the relaxed difference density matrix
        ModuleBase::matrix cal_force_hamilt_gs_dm_relaxed_diff(const elecstate::DensityMatrix<TK, double>& relaxed_diff_dm,
            const elecstate::DensityMatrix<TK, double>& dm_gs, const bool reproduce_gs = false,
            const PotHxcLR* pot_hxc_gs = nullptr);

        /// 2. $Tr[S^x * (EDM)]
        ModuleBase::matrix cal_force_overlap_edm(const elecstate::DensityMatrix<TK, double>& edm);

        /// 3. $\sum_{mnkl}(mn|f_{Hxc}|kl)^x *D^X *D^X$
        ModuleBase::matrix cal_force_hxc_dmtrans(const elecstate::DensityMatrix<TK, double>& dm_trans, const PotHxcLR& pot_hxc);

        /// 3b. the $g^{xc}$ half of $\partial_x K^{S/T}[D^X]D^X$:
        ///     $\int v^{(2)}[\rho^X,\rho^X](r)\,\partial_x\rho^\text{gs}(r)|_\text{basis}$
        ModuleBase::matrix cal_force_gxc_dmtrans(const elecstate::DensityMatrix<TK, double>& dm_trans,
            const elecstate::DensityMatrix<TK, double>& dm_gs, const PotGradXCLR& pot_grad);

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
            const elecstate::DensityMatrix<TK, double>& dm_gs,
            const elecstate::DensityMatrix<TK, double>& edm_gs);

        /// repreduce the ground state local term
        ModuleBase::matrix reproduce_force_gs_loc(const elecstate::DensityMatrix<TK, double>& dm_gs,
            const elecstate::Potential& pot_gs);

        /// derivatives of 2-center integrates: dtau(S_ij) and dtau(h_{ij}) (set vh_in_h=0)
        void cal_H2_sz_center2_deriv(const std::vector<double>& orb_cutoffs, const K_Vectors& kv);
        /// 4-center integrates or their derivatives: (ij | kl) or dtau(ij | kl) (set vl_in_h=0)
        void cal_H2_sz_center4(const std::vector<double>& orb_cutoffs,
            const K_Vectors& kv, const bool is_grad = false);

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

        Charge dm_to_charge(const elecstate::DensityMatrix<TK, double>& dm);
        elecstate::Potential dm_to_hxc_potential(const elecstate::DensityMatrix<TK, double>& dm);
        elecstate::Potential local_potential();
    };
}
