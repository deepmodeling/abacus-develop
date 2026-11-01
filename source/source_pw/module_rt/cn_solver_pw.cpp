#include "source_pw/module_rt/cn_solver_pw.h"

#include "source_base/module_device/memory_op.h"
#include "source_base/parallel_device.h"
#include "source_base/timer.h"
#include "source_base/tool_quit.h"
#include "source_basis/module_pw/pw_momentum.h"
#include "source_hsolver/module_linalg/kernels/linear_op.h"
#include "source_hsolver/module_linalg/orb_buffer.h"

#include <iomanip>
#include <iostream>
#include <sstream>

namespace hsolver
{
namespace
{

struct LinearMethodName
{
    const char* name;
    LinearMethod method;
};

constexpr LinearMethodName linear_method_names[]
    = {{"bicgstab", LinearMethod::bicgstab}, {"cgs", LinearMethod::cgs}, {"gmres", LinearMethod::gmres}};

struct PWPrecondName
{
    const char* name;
    PWPreconditioner preconditioner;
};

constexpr PWPrecondName pw_precond_names[] = {{"none", PWPreconditioner::none},
                                              {"kinetic", PWPreconditioner::kinetic},
                                              {"kinetic_recycle", PWPreconditioner::kinetic_recycle},
                                              {"kinetic_subspace", PWPreconditioner::kinetic_subspace}};

template <typename T, typename Device>
class ShiftedHOperator final : public LinearOperator<T, Device>
{
  public:
    ShiftedHOperator(const HSOperator<T, Device>& op, const T coefficient, const int dim, const int npw, const int npol)
        : op_(op), coefficient_(coefficient), dim_(dim), npw_(npw), npol_(npol)
    {
    }
    void apply(const T* x, T* y, const int ld, const int nvec) const override
    {
        op_.hpsi(x, y, ld, nvec);
        // Keep H zero on padding so the shifted operator retains only its identity term there.
        if (npol_ == 2)
        {
            const OrbLayout layout{ld, npw_, npol_, nvec};
            zero_pad<T, Device>(y, layout);
        }
        linear_op<T, Device>().batch(ld, dim_, nvec, y, x, y, T(1), coefficient_, nullptr, nullptr, nullptr);
    }

  private:
    const HSOperator<T, Device>& op_;
    const T coefficient_;
    const int dim_;
    const int npw_;
    const int npol_;
};

} // namespace

LinearMethod parse_linear_method(const std::string& name)
{
    for (const LinearMethodName& entry: linear_method_names)
    {
        if (name == entry.name)
        {
            return entry.method;
        }
    }
    const std::string message = "Unsupported linear solver: " + name;
    ModuleBase::WARNING_QUIT("CNSolverPW", message);
    return LinearMethod::bicgstab;
}

PWPreconditioner parse_pw_precond(const std::string& name)
{
    for (const PWPrecondName& entry: pw_precond_names)
    {
        if (name == entry.name)
        {
            return entry.preconditioner;
        }
    }
    const std::string message = "Unsupported preconditioner: " + name;
    ModuleBase::WARNING_QUIT("CNSolverPW", message);
    return PWPreconditioner::kinetic;
}

template <typename T, typename Device>
CNSolverPW<T, Device>::CNSolverPW(const ModulePW::PW_Basis_K& basis,
                                  const CNControls& cn_ctrl,
                                  int npol,
                                  const std::vector<int>& global_k,
                                  const diag_comm_info& comm,
                                  std::ostream& log)
    : basis_(basis), comm_(comm), cn_ctrl_(cn_ctrl), npol_(npol), global_k_(global_k), algebra_(comm), log_(log)
{
    initialize();
}

template <typename T, typename Device>
void CNSolverPW<T, Device>::initialize()
{
    const char* method = nullptr;
    for (const LinearMethodName& entry: linear_method_names)
    {
        if (cn_ctrl_.linear.solver == entry.method)
        {
            method = entry.name;
            break;
        }
    }
    if (method == nullptr)
    {
        ModuleBase::WARNING_QUIT("CNSolverPW", "Unsupported linear solver.");
        return;
    }
    const char* preconditioner = nullptr;
    for (const PWPrecondName& entry: pw_precond_names)
    {
        if (cn_ctrl_.preconditioner == entry.preconditioner)
        {
            preconditioner = entry.name;
            break;
        }
    }
    if (preconditioner == nullptr)
    {
        ModuleBase::WARNING_QUIT("CNSolverPW", "Unsupported preconditioner.");
        return;
    }
    if (cn_ctrl_.linear.solver != LinearMethod::gmres)
    {
        // Ignore the GMRES-only option before creating solver state or reporting effective settings.
        cn_ctrl_.linear.reconstruct = false;
    }
    linear_solver_.reset(new HSolverLinear<T, Device>(cn_ctrl_.linear, comm_));
    std::ostringstream info;
    info << " PW RT-TDDFT linear solver: " << method << "; tolerance: " << std::scientific << std::setprecision(6)
         << linear_solver_->tolerance() << "; maximum iterations: " << cn_ctrl_.linear.max_iterations
         << "; restart: " << cn_ctrl_.linear.restart << '\n'
         << "   Preconditioner: " << preconditioner << "; CN initial guess: " << cn_ctrl_.cn_init
         << "; residual reconstruction: " << cn_ctrl_.linear.reconstruct << '\n';
    log_ << info.str();
}

template <typename T, typename Device>
void CNSolverPW<T, Device>::prepare_buffers(const int nbands, const int nbasis)
{
    ModuleBase::timer::start("CNSolverPW", "prepare_buffers");
    using CtDevice = typename ct::PsiToContainer<Device>::type;
    const ct::DeviceType device = ct::DeviceTypeToEnum<CtDevice>::value;
    const int64_t elements = static_cast<int64_t>(nbands) * nbasis;
    const int64_t size = std::max<int64_t>(1, elements);
    if (rhs_.NumElements() < size || rhs_.data_type() != ct::DataTypeToEnum<T>::value || rhs_.device_type() != device)
    {
        rhs_ = ct::Tensor(ct::DataTypeToEnum<T>::value, device, {size});
        clean_ = ct::Tensor(ct::DataTypeToEnum<T>::value, device, {size});
    }
    ModuleBase::timer::end("CNSolverPW", "prepare_buffers");
}

template <typename T, typename Device>
void CNSolverPW<T, Device>::update_precond(const int ik,
                                           const int dim,
                                           const T coefficient,
                                           const ModuleBase::Vector3<double>& momentum_shift)
{
    ModuleBase::timer::start("CNSolverPW", "update_precond");
    using CtDevice = typename ct::PsiToContainer<Device>::type;
    const ct::DeviceType device = ct::DeviceTypeToEnum<CtDevice>::value;
    const int elements = std::max(1, dim);
    if (inverse_kinetic_.NumElements() < elements || inverse_kinetic_.data_type() != ct::DataTypeToEnum<T>::value
        || inverse_kinetic_.device_type() != device)
    {
        inverse_kinetic_ = ct::Tensor(ct::DataTypeToEnum<T>::value, device, {elements});
    }
    std::vector<T> inverse(dim, T(1));
    for (int ig = 0; ig < dim; ++ig)
    {
        // The two spin components share the same kinetic diagonal; padding keeps an identity preconditioner.
        const int physical_g = npol_ == 2 ? ig % basis_.npwk_max : ig;
        if (physical_g >= basis_.npwk[ik])
        {
            continue;
        }
        const double kinetic = cn_ctrl_.kinetic_enabled ? ModulePW::shifted_kinetic(basis_, ik, physical_g, momentum_shift) : 0.0;
        inverse[ig] = T(1) / (T(1) + coefficient * static_cast<Real>(kinetic));
    }
    if (dim > 0)
    {
        base_device::memory::synchronize_memory_op<T, Device, base_device::DEVICE_CPU>()(inverse_kinetic_.template data<T>(),
                                                                                         inverse.data(),
                                                                                         dim);
    }
    ModuleBase::timer::end("CNSolverPW", "update_precond");
}

template <typename T, typename Device>
bool CNSolverPW<T, Device>::tracks_state() const
{
    return cn_ctrl_.preconditioner == PWPreconditioner::kinetic_recycle || cn_ctrl_.linear.reconstruct;
}

template <typename T, typename Device>
void CNSolverPW<T, Device>::prepare_sequence(int nk, int ld, int bands, double dt, int step, int iteration)
{
    if (!tracks_state())
    {
        return;
    }
    const bool reset = states_.size() != static_cast<size_t>(nk) || dt != sequence_.dt || ld != sequence_.ld || bands != sequence_.bands
                       || step < sequence_.step || (step == sequence_.step && iteration <= sequence_.iteration);
    if (reset)
    {
        states_.clear();
        states_.resize(nk);
    }
    sequence_.dt = dt;
    sequence_.ld = ld;
    sequence_.bands = bands;
    sequence_.step = step;
    sequence_.iteration = iteration;
}

template <typename T, typename Device>
void CNSolverPW<T, Device>::invalidate_basis()
{
    for (KPointState& state: states_)
    {
        state.response.clear();
        state.solve_count = 0;
        state.independent_step = -1;
    }
}

template <typename T, typename Device>
bool CNSolverPW<T, Device>::require_audit(int step, KPointState* state) const
{
    if (!cn_ctrl_.linear.reconstruct)
    {
        return true;
    }
    constexpr unsigned int audit_period = 16;
    const bool required = state->solve_count % audit_period == 0 || state->independent_step == step;
    ++state->solve_count;
    return required;
}

template <typename T, typename Device>
void CNSolverPW<T, Device>::retry_kinetic(const LinearOperator<T, Device>& op,
                                          T* current,
                                          const BatchLayout& batch,
                                          LinearSolveResult* result)
{
    const int remaining_iter = cn_ctrl_.linear.max_iterations - result->iterations;
    const T* inverse = inverse_kinetic_.template data<T>();
    const LinearLowRank<T, Device> diagonal(algebra_, inverse, batch.dim);
    const T* rhs = rhs_.template data<T>();
    LinearSolveResult retry
        = linear_solver_->solve(op, diagonal, batch.ld, batch.bands, batch.dim, current, rhs, nullptr, true, remaining_iter, false);
    retry.iterations += result->iterations;
    retry.restarts += result->restarts;
    retry.operator_calls += result->operator_calls;
    retry.operator_columns += result->operator_columns;
    retry.true_checks += result->true_checks;
    retry.reconstruction_fallbacks += result->reconstruction_fallbacks;
    *result = retry;
}

template <typename T, typename Device>
void CNSolverPW<T, Device>::update_state(KPointState* state,
                                         const SolveDetails& details,
                                         int step,
                                         const BatchLayout& batch,
                                         const T* current)
{
    if (!state)
    {
        return;
    }
    if (details.retried || details.linear.reconstruction_fallbacks > 0)
    {
        state->response.clear();
    }
    if (details.linear.reconstruction_fallbacks > 0)
    {
        state->independent_step = step;
    }
    if (cn_ctrl_.preconditioner == PWPreconditioner::kinetic_recycle && details.projected
        && details.linear.status == LinearSolveStatus::converged && details.linear.reconstruction_fallbacks == 0 && !details.retried)
    {
        const T* seed = projection_.seed();
        const T* residual = projection_.residual();
        const double tolerance = linear_solver_->tolerance();
        state->response.update(algebra_, batch.ld, batch.dim, batch.bands, current, seed, residual, tolerance, &response_workspace_);
    }
}

template <typename T, typename Device>
typename CNSolverPW<T, Device>::SolveDetails CNSolverPW<T, Device>::solve_kpoint(const LinearOperator<T, Device>& op,
                                                                                 const T* previous,
                                                                                 T* current,
                                                                                 const BatchLayout& batch,
                                                                                 int step,
                                                                                 int iteration,
                                                                                 KPointState* state)
{
    SolveDetails details;
    const T* rhs = rhs_.template data<T>();
    const bool need_projection = (cn_ctrl_.cn_init && iteration == 1) || cn_ctrl_.preconditioner == PWPreconditioner::kinetic_recycle
                                 || cn_ctrl_.preconditioner == PWPreconditioner::kinetic_subspace;
    details.projected = need_projection && projection_.prepare(algebra_, batch.ld, batch.dim, batch.bands, previous, rhs);
    const T* initial_residual = nullptr;
    if (cn_ctrl_.cn_init && iteration == 1 && details.projected)
    {
        const T* seed = projection_.seed();
        if (batch.dim > 0)
        {
            base_device::memory::synchronize_memory_2d_op<T, Device, Device>()(current, batch.ld, seed, batch.ld, batch.dim, batch.bands);
        }
        initial_residual = projection_.residual();
        details.cn_initial = true;
    }
    const bool force_check = require_audit(step, state);
    const T* inverse = cn_ctrl_.preconditioner == PWPreconditioner::none ? nullptr : inverse_kinetic_.template data<T>();
    LinearLowRank<T, Device> preconditioner(algebra_, inverse, batch.dim);
    if (cn_ctrl_.preconditioner == PWPreconditioner::kinetic_subspace && details.projected)
    {
        const T* image = projection_.image();
        const LinearSmallLU& factor = projection_.factor();
        preconditioner.prepare_subspace(batch.ld, batch.bands, previous, image, factor, &correction_workspace_);
    }
    else if (cn_ctrl_.preconditioner == PWPreconditioner::kinetic_recycle && state->response.rank() > 0)
    {
        const int rank = state->response.rank();
        const T* directions = state->response.directions();
        const T* images = state->response.images();
        preconditioner.prepare_response(batch.ld, rank, directions, images, &correction_workspace_);
    }
    details.coarse_rank = preconditioner.rank();
    details.linear
        = linear_solver_->solve(op, preconditioner, batch.ld, batch.bands, batch.dim, current, rhs, initial_residual, force_check);
    if (details.linear.status != LinearSolveStatus::converged && details.coarse_rank > 0
        && details.linear.iterations < cn_ctrl_.linear.max_iterations)
    {
        retry_kinetic(op, current, batch, &details.linear);
        details.retried = true;
    }
    update_state(state, details, step, batch, current);
    return details;
}

template <typename T, typename Device>
void CNSolverPW<T, Device>::report_solve(const SolveDetails& details, int ik, int step, int iteration) const
{
    const LinearSolveResult& result = details.linear;
    const int evolution_step = step + 1;
    std::ostringstream record;
    record << std::scientific << std::setprecision(6) << " PW RT-TDDFT linear solve: evolution_step=" << evolution_step
           << " scf_iter=" << iteration << " global_k=" << global_k_.at(ik) + 1 << '\n'
           << "   iterations=" << result.iterations << " restarts=" << result.restarts << " operator_calls=" << result.operator_calls
           << " operator_columns=" << result.operator_columns << " pre_orth_residual=" << result.max_residual
           << " residual_kind=" << (result.reconstructed ? "reconstructed" : "independent") << " true_checks=" << result.true_checks << '\n'
           << "   reconstruction_fallbacks=" << result.reconstruction_fallbacks << " coarse_rank=" << details.coarse_rank
           << " kinetic_retry=" << details.retried << " cn_projected=" << details.projected << " cn_initial=" << details.cn_initial << '\n';
    log_ << record.str();
}

template <typename T, typename Device>
void CNSolverPW<T, Device>::solve(HSOperator<T, Device>& op,
                                  const psi::Psi<T, Device>& previous,
                                  psi::Psi<T, Device>* current,
                                  const double dt,
                                  const ModuleBase::Vector3<double>& momentum_shift,
                                  const int istep,
                                  const int iter,
                                  const bool detailed_output)
{
    ModuleBase::timer::start("CNSolverPW", "solve");
    const int bands = current->get_nbands();
    const int ld = current->get_nbasis();
    const int nk = current->get_nk();
    prepare_buffers(bands, ld);
    prepare_sequence(nk, ld, bands, dt, istep, iter);
    T* rhs = rhs_.template data<T>();
    // H is stored in Rydberg; CN includes the conversion to Hartree.
    const T coefficient(0.0, dt / 4.0);
    const T rhs_coefficient = -coefficient;
    for (int ik = 0; ik < nk; ++ik)
    {
        op.update_k(ik);
        current->fix_k(ik);
        previous.fix_k(ik);
        const int npw = current->get_ngk(ik);
        // A spinor spans both npwx-strided blocks, including their explicitly zeroed padding.
        const int dim = npol_ == 2 ? ld : npw;
        KPointState* state = tracks_state() ? &states_[ik] : nullptr;
        const ShiftedHOperator<T, Device> rhs_op(op, rhs_coefficient, dim, npw, npol_);
        const ShiftedHOperator<T, Device> lhs_op(op, coefficient, dim, npw, npol_);
        const OrbLayout layout{ld, npw, npol_, bands};
        const T* previous_data = clean_orbs<T, Device>(previous.get_pointer(), layout, &clean_);
        T* current_data = current->get_pointer();
        if (npol_ == 2)
        {
            zero_pad<T, Device>(current_data, layout);
        }
        rhs_op.apply(previous_data, rhs, ld, bands);
        if (cn_ctrl_.preconditioner != PWPreconditioner::none)
        {
            update_precond(ik, dim, coefficient, momentum_shift);
        }
        const BatchLayout batch{ld, dim, bands};
        const SolveDetails details = solve_kpoint(lhs_op, previous_data, current_data, batch, istep, iter, state);
        if (detailed_output && log_.good())
        {
            report_solve(details, ik, istep, iter);
        }
        if (details.linear.status != LinearSolveStatus::converged)
        {
            const LinearSolveResult& result = details.linear;
            std::ostringstream message;
            const int evolution_step = istep + 1;
            std::string band_label = "unknown";
            if (result.failed_band >= 0)
            {
                band_label = std::to_string(result.failed_band + 1);
            }
            message << "PW RT-TDDFT linear solve failed at electronic evolution step " << evolution_step << ", SCF iteration " << iter
                    << ", global_k=" << global_k_.at(ik) + 1 << ", band " << band_label << ", after " << result.iterations
                    << " iterations: " << linear_status_name(result.status) << "; residual = " << result.max_residual;
            ModuleBase::WARNING_QUIT("CNSolverPW", message.str());
        }
    }
    ModuleBase::timer::end("CNSolverPW", "solve");
}

template class CNSolverPW<std::complex<float>, base_device::DEVICE_CPU>;
template class CNSolverPW<std::complex<double>, base_device::DEVICE_CPU>;
#if defined(__CUDA) || defined(__ROCM)
template class CNSolverPW<std::complex<float>, base_device::DEVICE_GPU>;
template class CNSolverPW<std::complex<double>, base_device::DEVICE_GPU>;
#endif
} // namespace hsolver
