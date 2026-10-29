#include "source_pw/module_rt/td_prop_pw.h"

#include "source_pw/module_rt/prop_cn_pw.h"
#include "source_pw/module_rt/prop_etrs_pw.h"

#include <utility>

namespace pw
{
template <typename T, typename Device>
std::unique_ptr<TDPropPW<T, Device>> make_td_prop(const ModulePW::PW_Basis_K& basis,
                                                  const Input_para& input,
                                                  int npol,
                                                  const std::vector<int>& global_k,
                                                  const hsolver::diag_comm_info& comm,
                                                  std::ostream& log)
{
    const hsolver::OrthMethod orth_method = hsolver::parse_orth_method(input.td_orthonormal);
    std::unique_ptr<hsolver::TDOrbOps<T, Device>> orb_ops(
        new hsolver::TDOrbOps<T, Device>(npol, orth_method, input.td_out_stat, global_k, comm, log));
    std::unique_ptr<PWProp<T, Device>> prop;
    if (input.propagator == 0)
    {
        hsolver::CNControls cn_ctrl;
        cn_ctrl.linear.solver = hsolver::parse_linear_method(input.lin_solver);
        cn_ctrl.linear.tolerance = input.lin_thr;
        cn_ctrl.linear.max_iterations = input.lin_maxiter;
        cn_ctrl.linear.restart = input.lin_gmres_restart;
        cn_ctrl.linear.reconstruct = input.lin_reconstruct;
        cn_ctrl.preconditioner = hsolver::parse_pw_precond(input.lin_precond);
        cn_ctrl.cn_init = input.td_cn_init;
        cn_ctrl.kinetic_enabled = input.t_in_h;
        prop.reset(new CNPropPW<T, Device>(basis, cn_ctrl, npol, global_k, comm, log));
    }
    else
    {
        hsolver::ExpControls exp_ctrl;
        if (input.td_exp_method == "taylor")
        {
            exp_ctrl.expansion = hsolver::ExpMethod::taylor;
        }
        exp_ctrl.taylor_order = input.td_exp_order;
        exp_ctrl.rel_tol = input.td_exp_thr;
        exp_ctrl.krylov_max = input.td_exp_krylov_dim;
        prop.reset(new ETRSPropPW<T, Device>(exp_ctrl, npol, global_k, comm, log));
    }

    return std::unique_ptr<TDPropPW<T, Device>>(new TDPropPW<T, Device>(basis, std::move(orb_ops), std::move(prop)));
}

template <typename T, typename Device>
TDPropPW<T, Device>::TDPropPW(const ModulePW::PW_Basis_K& basis,
                              std::unique_ptr<hsolver::TDOrbOps<T, Device>> orb_ops,
                              std::unique_ptr<PWProp<T, Device>> prop)
    : basis_(basis), orb_ops_(std::move(orb_ops)), prop_(std::move(prop))
{
}

template <typename T, typename Device>
PropStage TDPropPW<T, Device>::stage(int scf_iter) const
{
    return scf_iter == 1 ? PropStage::predictor : PropStage::corrector;
}
template <typename T, typename Device>
void TDPropPW<T, Device>::begin_trial()
{
    prop_->begin_trial();
}
template <typename T, typename Device>
void TDPropPW<T, Device>::finish_trial()
{
    prop_->discard_trial();
}
template <typename T, typename Device>
void TDPropPW<T, Device>::discard_trial()
{
    prop_->discard_trial();
}
template <typename T, typename Device>
void TDPropPW<T, Device>::invalidate_basis()
{
    prop_->invalidate_basis();
}
template <typename T, typename Device>
void TDPropPW<T, Device>::evolve(hamilt::HamiltPW<T, Device>& hamilt,
                                 elecstate::Potential& potential,
                                 const TDHistoryPW<T, Device>& history,
                                 const elecstate::TDFieldManager& field,
                                 psi::Psi<T, Device>* current,
                                 int elec_step,
                                 int scf_iter,
                                 bool log_detail)
{
    using Real = typename T::value_type;
    TDHamiltPW<T, Device> td_hamilt(hamilt, basis_);
    prop_->evolve(td_hamilt, history, potential, field, current, elec_step, scf_iter, log_detail);
    orb_ops_->correct(current, elec_step, scf_iter, log_detail);
    const Real* endpoint = potential.template get_veff_smooth_data<Real>();
    const ModuleBase::Vector3<double>& vecpot_ha = field.A_right_ha();
    td_hamilt.bind(endpoint, vecpot_ha);
}
template class TDPropPW<std::complex<float>, base_device::DEVICE_CPU>;
template std::unique_ptr<TDPropPW<std::complex<float>, base_device::DEVICE_CPU>> make_td_prop(const ModulePW::PW_Basis_K&,
                                                                                              const Input_para&,
                                                                                              int,
                                                                                              const std::vector<int>&,
                                                                                              const hsolver::diag_comm_info&,
                                                                                              std::ostream&);
template class TDPropPW<std::complex<double>, base_device::DEVICE_CPU>;
template std::unique_ptr<TDPropPW<std::complex<double>, base_device::DEVICE_CPU>> make_td_prop(const ModulePW::PW_Basis_K&,
                                                                                               const Input_para&,
                                                                                               int,
                                                                                               const std::vector<int>&,
                                                                                               const hsolver::diag_comm_info&,
                                                                                               std::ostream&);
#if defined(__CUDA) || defined(__ROCM)
template class TDPropPW<std::complex<float>, base_device::DEVICE_GPU>;
template std::unique_ptr<TDPropPW<std::complex<float>, base_device::DEVICE_GPU>> make_td_prop(const ModulePW::PW_Basis_K&,
                                                                                              const Input_para&,
                                                                                              int,
                                                                                              const std::vector<int>&,
                                                                                              const hsolver::diag_comm_info&,
                                                                                              std::ostream&);
template class TDPropPW<std::complex<double>, base_device::DEVICE_GPU>;
template std::unique_ptr<TDPropPW<std::complex<double>, base_device::DEVICE_GPU>> make_td_prop(const ModulePW::PW_Basis_K&,
                                                                                               const Input_para&,
                                                                                               int,
                                                                                               const std::vector<int>&,
                                                                                               const hsolver::diag_comm_info&,
                                                                                               std::ostream&);
#endif
} // namespace pw
