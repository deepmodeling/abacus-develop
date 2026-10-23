#include "source_pw/module_rt/prop_etrs_pw.h"

#include "source_base/parallel_device.h"
#include "source_base/tool_quit.h"

#include <iostream>
#include <sstream>

namespace pw
{
template <typename T, typename Device>
ETRSPropPW<T, Device>::ETRSPropPW(const TDPropOpts& prop_opts, const hsolver::diag_comm_info& comm, std::ostream& log)
    : pool_comm_(comm), orb_opts_(prop_opts.orb), log_(log)
{
    exponential_ = hsolver::make_exp_action<T, Device>(prop_opts.exp, comm);
    log_ << " PW RT-TDDFT propagator: ETRS; exponential: ";
    if (prop_opts.exp.method == hsolver::ExpMethod::taylor)
    {
        log_ << "taylor; order=" << prop_opts.exp.taylor_order << "; predictor=shared powers\n";
    }
    else
    {
        const double rel_tol = exponential_->rel_tol();
        log_ << "lanczos; tolerance=" << rel_tol << "; maximum dimension=" << prop_opts.exp.krylov_max << '\n';
    }
}
template <typename T, typename Device>
void ETRSPropPW<T, Device>::begin_trial()
{
    half_ready_ = false;
}
template <typename T, typename Device>
void ETRSPropPW<T, Device>::discard_trial()
{
    half_ready_ = false;
}
template <typename T, typename Device>
void ETRSPropPW<T, Device>::invalidate_basis()
{
    half_ready_ = false;
}
template <typename T, typename Device>
void ETRSPropPW<T, Device>::report(const hsolver::ExpResult& result, const PropArgs<T, Device>& args, int ik) const
{
    std::ostringstream text;
    text << " PW ETRS exponential: step=" << args.elec_step << " scf_iter=" << args.scf_iter
         << " global_k=" << orb_opts_.global_k_indices.at(ik) + 1 << " order_or_dimension=" << result.order_or_dim
         << " operator_calls=" << result.hpsi_calls << " operator_columns=" << result.hpsi_cols << " subdivisions=" << result.bisect_count;
    if (result.has_err_est)
    {
        text << " estimated_error=" << result.rel_err_est;
    }
    if (args.log_detail)
    {
        log_ << text.str() << '\n';
    }
    if (result.state != hsolver::ExpState::success)
    {
        text << " failed_band=" << result.failed_band + 1 << "; reduce td_dt or adjust exponential controls.";
        const std::string message = text.str();
        // stdout is disabled on non-world-root ranks; report from the failing pool.
        if (pool_comm_.rank == 0)
        {
            std::cerr << message << std::endl;
        }
#ifdef __MPI
        // Flush the pool root's diagnostic before any peer can terminate the MPI job.
        double diag_sent = 1.0;
        Parallel_Common::bcast_data(&diag_sent, 1, pool_comm_.comm, 0);
#endif
        ModuleBase::WARNING_QUIT("ETRSPropPW", message);
    }
}
template <typename T, typename Device>
void ETRSPropPW<T, Device>::evolve(PropArgs<T, Device>& args, psi::Psi<T, Device>* current)
{
    using Real = typename T::value_type;
    const bool predictor = args.scf_iter == 1;
    const auto& previous = args.history.previous();
    const int ld = current->get_nbasis();
    const int bands = current->get_nbands();
    const int64_t stride = static_cast<int64_t>(ld) * bands;
    if (predictor)
    {
        const int64_t size = stride * current->get_nk();
        if (half_psi_.NumElements() < size)
        {
            // Release obsolete storage before allocating a larger half-step buffer.
            half_psi_ = ct::Tensor();
        }
        hsolver::linear_buffer<T, Device>(&half_psi_, size);
        left_gen_ = args.history.generation();
        args.hamilt.bind(args.history.left_potential(), args.field.A_left_ha());
    }
    else
    {
        if (!half_ready_ || left_gen_ != args.history.generation())
        {
            ModuleBase::WARNING_QUIT("ETRSPropPW", "Missing left half-step for the current trial.");
        }
        const Real* endpoint = args.potential.template get_veff_smooth_data<Real>();
        args.hamilt.bind(endpoint, args.field.A_right_ha());
    }
    const double alpha = args.field.interval_dt_ha() / 4.0;
    for (int ik = 0; ik < current->get_nk(); ++ik)
    {
        args.hamilt.update_k(ik);
        current->fix_k(ik);
        previous.fix_k(ik);
        const hsolver::OrbLayout layout{ld, current->get_ngk(ik), orb_opts_.npol, bands};
        T* output = current->get_pointer();
        T* half = half_psi_.template data<T>() + ik * stride;
        hsolver::ExpResult result;
        if (predictor)
        {
            const T* input = previous.get_pointer();
            result = exponential_->apply_pair(args.hamilt, layout, alpha, input, half, output);
        }
        else
        {
            result = exponential_->apply(args.hamilt, layout, alpha, half, output);
        }
        report(result, args, ik);
    }
    half_ready_ = true;
}
template class ETRSPropPW<std::complex<float>, base_device::DEVICE_CPU>;
template class ETRSPropPW<std::complex<double>, base_device::DEVICE_CPU>;
#if defined(__CUDA) || defined(__ROCM)
template class ETRSPropPW<std::complex<float>, base_device::DEVICE_GPU>;
template class ETRSPropPW<std::complex<double>, base_device::DEVICE_GPU>;
#endif
} // namespace pw
