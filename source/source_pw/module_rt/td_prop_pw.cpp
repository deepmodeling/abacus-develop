#include "source_pw/module_rt/td_prop_pw.h"

#include "source_pw/module_rt/prop_cn_pw.h"
#include "source_pw/module_rt/prop_etrs_pw.h"

namespace pw
{
TDPropOpts td_prop_options(const Input_para& input, int npol, const std::vector<int>& global_k)
{
    TDPropOpts prop_opts;
    prop_opts.orb.npol = npol;
    prop_opts.orb.global_k_indices = global_k;
    prop_opts.orb.orthonormal = hsolver::parse_orth_method(input.td_orthonormal);
    prop_opts.orb.out_stat = input.td_out_stat;
    if (input.propagator == 0)
    {
        prop_opts.prop_kind = PropKind::cn;
        prop_opts.cn.npol = npol;
        prop_opts.cn.global_k_indices = global_k;
        prop_opts.cn.linear.method = hsolver::parse_linear_method(input.lin_solver);
        prop_opts.cn.linear.tolerance = input.lin_thr;
        prop_opts.cn.linear.max_iterations = input.lin_maxiter;
        prop_opts.cn.linear.restart = input.lin_gmres_restart;
        prop_opts.cn.linear.reconstruct = input.lin_reconstruct;
        prop_opts.cn.preconditioner = hsolver::parse_pw_precond(input.lin_precond);
        prop_opts.cn.cn_init = input.td_cn_init;
        prop_opts.cn.kinetic_enabled = input.t_in_h;
    }
    else
    {
        prop_opts.prop_kind = PropKind::etrs;
        if (input.td_exp_method == "taylor")
        {
            prop_opts.exp.method = hsolver::ExpMethod::taylor;
        }
        prop_opts.exp.taylor_order = input.td_exp_order;
        prop_opts.exp.rel_tol = input.td_exp_thr;
        prop_opts.exp.krylov_max = input.td_exp_krylov_dim;
    }
    return prop_opts;
}

template <typename T, typename Device>
TDPropPW<T, Device>::TDPropPW(const ModulePW::PW_Basis_K& basis,
                              const TDPropOpts& prop_opts,
                              const hsolver::diag_comm_info& comm,
                              std::ostream& log)
    : basis_(basis), orb_ops_(prop_opts.orb, comm, log)
{
    if (prop_opts.prop_kind == PropKind::cn)
    {
        prop_.reset(new CNPropPW<T, Device>(basis, prop_opts.cn, comm, log));
    }
    else
    {
        prop_.reset(new ETRSPropPW<T, Device>(prop_opts, comm, log));
    }
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
    PropArgs<T, Device> args{td_hamilt, history, potential, field, elec_step, scf_iter, log_detail};
    prop_->evolve(args, current);
    orb_ops_.correct(current, elec_step, scf_iter, log_detail);
    const Real* endpoint = potential.template get_veff_smooth_data<Real>();
    const auto& vecpot_ha = field.A_right_ha();
    td_hamilt.bind(endpoint, vecpot_ha);
}
template class TDPropPW<std::complex<float>, base_device::DEVICE_CPU>;
template class TDPropPW<std::complex<double>, base_device::DEVICE_CPU>;
#if defined(__CUDA) || defined(__ROCM)
template class TDPropPW<std::complex<float>, base_device::DEVICE_GPU>;
template class TDPropPW<std::complex<double>, base_device::DEVICE_GPU>;
#endif
} // namespace pw
