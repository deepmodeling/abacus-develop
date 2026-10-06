#include "source_hsolver/hsolver_pw_tddft.h"

#include "source_base/module_device/memory_op.h"
#include "source_base/parallel_device.h"
#include "source_base/timer.h"
#include "source_base/tool_quit.h"
#include "source_basis/module_pw/pw_momentum.h"
#include "source_hsolver/kernels/linear_op.h"

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
    ShiftedHOperator(const HSOperator<T, Device>& op, const T coefficient, const int dim) : op_(op), coefficient_(coefficient), dim_(dim)
    {
    }
    void apply(const T* x, T* y, const int ld, const int nvec) const override
    {
        op_.hpsi(x, y, ld, nvec);
        linear_op<T, Device>().batch(ld, dim_, nvec, y, x, y, T(1), coefficient_, nullptr, nullptr, nullptr);
    }

  private:
    const HSOperator<T, Device>& op_;
    const T coefficient_;
    const int dim_;
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
    ModuleBase::WARNING_QUIT("HSolverPWTDDFT", message);
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
    ModuleBase::WARNING_QUIT("HSolverPWTDDFT", message);
    return PWPreconditioner::kinetic;
}

template <typename T, typename Device>
HSolverPWTDDFT<T, Device>::HSolverPWTDDFT(const ModulePW::PW_Basis_K& basis,
                                          const PWLinearOptions& options,
                                          const diag_comm_info& comm,
                                          std::ostream& log)
    : basis_(basis), comm_(comm), options_(options), algebra_(comm), orthonormal_(comm), log_(log), band_products_(comm, 9)
{
    initialize(log);
    log << " PW RT-TDDFT orthonormalization: " << orth_method_name(options_.orthonormal) << '\n';
}

template <typename T, typename Device>
void HSolverPWTDDFT<T, Device>::initialize(std::ostream& log)
{
    const char* method = nullptr;
    for (const LinearMethodName& entry: linear_method_names)
    {
        if (options_.linear.method == entry.method)
        {
            method = entry.name;
            break;
        }
    }
    if (method == nullptr)
    {
        ModuleBase::WARNING_QUIT("HSolverPWTDDFT", "Unsupported linear solver.");
        return;
    }
    const char* preconditioner = nullptr;
    for (const PWPrecondName& entry: pw_precond_names)
    {
        if (options_.preconditioner == entry.preconditioner)
        {
            preconditioner = entry.name;
            break;
        }
    }
    if (preconditioner == nullptr)
    {
        ModuleBase::WARNING_QUIT("HSolverPWTDDFT", "Unsupported preconditioner.");
        return;
    }
    if (options_.linear.method != LinearMethod::gmres)
    {
        // Ignore the GMRES-only option before creating solver state or reporting effective settings.
        options_.linear.reconstruct = false;
    }
    linear_solver_.reset(new HSolverLinear<T, Device>(options_.linear, comm_));
    std::ostringstream info;
    info << " PW RT-TDDFT linear solver: " << method << "; tolerance: " << std::scientific << std::setprecision(6)
         << linear_solver_->tolerance() << "; maximum iterations: " << options_.linear.max_iterations
         << "; restart: " << options_.linear.restart << '\n'
         << "   Preconditioner: " << preconditioner << "; CN initial guess: " << options_.cn_init
         << "; residual reconstruction: " << options_.linear.reconstruct << '\n';
    log << info.str();
}

template <typename T, typename Device>
void HSolverPWTDDFT<T, Device>::prepare_buffers(const int nbands, const int nbasis)
{
    ModuleBase::timer::start("HSolverPWTDDFT", "prepare_buffers");
    using CtDevice = typename ct::PsiToContainer<Device>::type;
    const ct::DeviceType device = ct::DeviceTypeToEnum<CtDevice>::value;
    const int64_t elements = static_cast<int64_t>(nbands) * nbasis;
    const int64_t size = std::max<int64_t>(1, elements);
    if (rhs_.NumElements() < size || rhs_.data_type() != ct::DataTypeToEnum<T>::value || rhs_.device_type() != device)
    {
        rhs_ = ct::Tensor(ct::DataTypeToEnum<T>::value, device, {size});
        hpsi_ = ct::Tensor(ct::DataTypeToEnum<T>::value, device, {size});
    }
    ModuleBase::timer::end("HSolverPWTDDFT", "prepare_buffers");
}

template <typename T, typename Device>
void HSolverPWTDDFT<T, Device>::update_precond(const int ik,
                                               const int dim,
                                               const T coefficient,
                                               const ModuleBase::Vector3<double>& momentum_shift)
{
    ModuleBase::timer::start("HSolverPWTDDFT", "update_precond");
    using CtDevice = typename ct::PsiToContainer<Device>::type;
    const ct::DeviceType device = ct::DeviceTypeToEnum<CtDevice>::value;
    const int elements = std::max(1, dim);
    if (inverse_kinetic_.NumElements() < elements || inverse_kinetic_.data_type() != ct::DataTypeToEnum<T>::value
        || inverse_kinetic_.device_type() != device)
    {
        inverse_kinetic_ = ct::Tensor(ct::DataTypeToEnum<T>::value, device, {elements});
    }
    std::vector<T> inverse(dim);
    for (int ig = 0; ig < dim; ++ig)
    {
        const double kinetic = options_.kinetic_enabled ? ModulePW::shifted_kinetic(basis_, ik, ig, momentum_shift) : 0.0;
        inverse[ig] = T(1) / (T(1) + coefficient * static_cast<Real>(kinetic));
    }
    if (dim > 0)
    {
        base_device::memory::synchronize_memory_op<T, Device, base_device::DEVICE_CPU>()(inverse_kinetic_.template data<T>(),
                                                                                         inverse.data(),
                                                                                         dim);
    }
    ModuleBase::timer::end("HSolverPWTDDFT", "update_precond");
}

template <typename T, typename Device>
bool HSolverPWTDDFT<T, Device>::tracks_state() const
{
    return options_.preconditioner == PWPreconditioner::kinetic_recycle || options_.linear.reconstruct;
}

template <typename T, typename Device>
void HSolverPWTDDFT<T, Device>::prepare_sequence(int nk, int ld, int bands, double dt, int step, int iteration)
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
void HSolverPWTDDFT<T, Device>::invalidate_basis()
{
    for (KPointState& state: states_)
    {
        state.response.clear();
        state.solve_count = 0;
        state.independent_step = -1;
    }
}

template <typename T, typename Device>
bool HSolverPWTDDFT<T, Device>::require_audit(int step, KPointState* state) const
{
    if (!options_.linear.reconstruct)
    {
        return true;
    }
    constexpr unsigned int audit_period = 16;
    const bool required = state->solve_count % audit_period == 0 || state->independent_step == step;
    ++state->solve_count;
    return required;
}

template <typename T, typename Device>
void HSolverPWTDDFT<T, Device>::retry_kinetic(const LinearOperator<T, Device>& op,
                                              T* current,
                                              const SolveBatch& batch,
                                              LinearSolveResult* result)
{
    const int remaining_iterations = options_.linear.max_iterations - result->iterations;
    const LinearSolveControl control{remaining_iterations, false};
    const T* inverse = inverse_kinetic_.template data<T>();
    const LinearLowRank<T, Device> diagonal(algebra_, inverse, batch.dim);
    const T* rhs = rhs_.template data<T>();
    LinearSolveResult retry = linear_solver_->solve(op, diagonal, batch.ld, batch.bands, batch.dim, current, rhs, nullptr, true, control);
    retry.iterations += result->iterations;
    retry.restarts += result->restarts;
    retry.operator_calls += result->operator_calls;
    retry.operator_columns += result->operator_columns;
    retry.true_checks += result->true_checks;
    retry.reconstruction_fallbacks += result->reconstruction_fallbacks;
    *result = retry;
}

template <typename T, typename Device>
void HSolverPWTDDFT<T, Device>::update_state(KPointState* state,
                                             const SolveDetails& details,
                                             int step,
                                             const SolveBatch& batch,
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
    if (options_.preconditioner == PWPreconditioner::kinetic_recycle && details.projected
        && details.linear.status == LinearSolveStatus::converged && details.linear.reconstruction_fallbacks == 0 && !details.retried)
    {
        const T* seed = projection_.seed();
        const T* residual = projection_.residual();
        const double tolerance = linear_solver_->tolerance();
        state->response.update(algebra_, batch.ld, batch.dim, batch.bands, current, seed, residual, tolerance, &response_workspace_);
    }
}

template <typename T, typename Device>
typename HSolverPWTDDFT<T, Device>::SolveDetails HSolverPWTDDFT<T, Device>::solve_kpoint(const LinearOperator<T, Device>& op,
                                                                                         const T* previous,
                                                                                         T* current,
                                                                                         const SolveBatch& batch,
                                                                                         int step,
                                                                                         int iteration,
                                                                                         KPointState* state)
{
    SolveDetails details;
    const T* rhs = rhs_.template data<T>();
    const bool need_projection = (options_.cn_init && iteration == 1) || options_.preconditioner == PWPreconditioner::kinetic_recycle
                                 || options_.preconditioner == PWPreconditioner::kinetic_subspace;
    details.projected = need_projection && projection_.prepare(algebra_, batch.ld, batch.dim, batch.bands, previous, rhs);
    const T* initial_residual = nullptr;
    if (options_.cn_init && iteration == 1 && details.projected)
    {
        const T* seed = projection_.seed();
        band_products_.copy(batch.ld, batch.dim, batch.bands, seed, current);
        initial_residual = projection_.residual();
        details.cn_initial = true;
    }
    const bool force_check = require_audit(step, state);
    const T* inverse = options_.preconditioner == PWPreconditioner::none ? nullptr : inverse_kinetic_.template data<T>();
    LinearLowRank<T, Device> preconditioner(algebra_, inverse, batch.dim);
    if (options_.preconditioner == PWPreconditioner::kinetic_subspace && details.projected)
    {
        const T* image = projection_.image();
        const LinearSmallLU& factor = projection_.factor();
        preconditioner.prepare_subspace(batch.ld, batch.bands, previous, image, factor, &correction_workspace_);
    }
    else if (options_.preconditioner == PWPreconditioner::kinetic_recycle && state->response.rank() > 0)
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
        && details.linear.iterations < options_.linear.max_iterations)
    {
        retry_kinetic(op, current, batch, &details.linear);
        details.retried = true;
    }
    update_state(state, details, step, batch, current);
    return details;
}

template <typename T, typename Device>
void HSolverPWTDDFT<T, Device>::report_solve(const SolveDetails& details, int ik, int step, int iteration, std::ostream& log) const
{
    const LinearSolveResult& result = details.linear;
    const int evolution_step = step + 1;
    std::ostringstream record;
    record << std::scientific << std::setprecision(6) << " PW RT-TDDFT linear solve: evolution_step=" << evolution_step
           << " scf_iter=" << iteration << " k=" << ik << '\n'
           << "   iterations=" << result.iterations << " restarts=" << result.restarts << " operator_calls=" << result.operator_calls
           << " operator_columns=" << result.operator_columns << " pre_orth_residual=" << result.max_residual
           << " residual_kind=" << (result.reconstructed ? "reconstructed" : "independent") << " true_checks=" << result.true_checks << '\n'
           << "   reconstruction_fallbacks=" << result.reconstruction_fallbacks << " coarse_rank=" << details.coarse_rank
           << " kinetic_retry=" << details.retried << " cn_projected=" << details.projected << " cn_initial=" << details.cn_initial << '\n';
    log << record.str();
}

template <typename T, typename Device>
void HSolverPWTDDFT<T, Device>::solve(HSOperator<T, Device>& op,
                                      const psi::Psi<T, Device>& previous,
                                      psi::Psi<T, Device>* current,
                                      const double dt,
                                      const ModuleBase::Vector3<double>& momentum_shift,
                                      const int istep,
                                      const int iter,
                                      const bool detailed_output,
                                      std::ostream& log)
{
    ModuleBase::timer::start("HSolverPWTDDFT", "solve");
    const int bands = current->get_nbands();
    const int ld = current->get_nbasis();
    const int nk = current->get_nk();
    if (options_.out_stat)
    {
        orth_norms_.resize(nk);
    }
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
        const int dim = current->get_ngk(ik);
        KPointState* state = tracks_state() ? &states_[ik] : nullptr;
        const ShiftedHOperator<T, Device> rhs_op(op, rhs_coefficient, dim);
        const ShiftedHOperator<T, Device> lhs_op(op, coefficient, dim);
        const T* previous_data = previous.get_pointer();
        T* current_data = current->get_pointer();
        rhs_op.apply(previous_data, rhs, ld, bands);
        if (options_.preconditioner != PWPreconditioner::none)
        {
            update_precond(ik, dim, coefficient, momentum_shift);
        }
        const SolveBatch batch{ld, dim, bands};
        const SolveDetails details = solve_kpoint(lhs_op, previous_data, current_data, batch, istep, iter, state);
        if (detailed_output && log.good())
        {
            report_solve(details, ik, istep, iter, log);
        }
        if (details.linear.status != LinearSolveStatus::converged)
        {
            const LinearSolveResult& result = details.linear;
            std::ostringstream message;
            const int evolution_step = istep + 1;
            message << "PW RT-TDDFT linear solve failed at electronic evolution step " << evolution_step << ", SCF iteration " << iter
                    << ", k point " << ik << ", band " << result.failed_band << ", after " << result.iterations
                    << " iterations: " << linear_status_name(result.status) << "; residual = " << result.max_residual;
            ModuleBase::WARNING_QUIT("HSolverPWTDDFT", message.str());
        }
        correct_orbitals(current_data, ld, dim, bands, ik, istep, iter, detailed_output);
    }
    ModuleBase::timer::end("HSolverPWTDDFT", "solve");
}

template <typename T, typename Device>
void HSolverPWTDDFT<T, Device>::cal_band_energy(HSOperator<T, Device>& op, const psi::Psi<T, Device>& current, ModuleBase::matrix* energies)
{
    ModuleBase::timer::start("HSolverPWTDDFT", "cal_band_energy");
    const int nband = current.get_nbands();
    const int ld = current.get_nbasis();
    prepare_buffers(nband, ld);
    std::vector<T> expectations(nband);
    T* hpsi = hpsi_.template data<T>();
    for (int ik = 0; ik < current.get_nk(); ++ik)
    {
        op.update_k(ik);
        current.fix_k(ik);
        op.hpsi(current.get_pointer(), hpsi, ld, nband);
        band_products_.dot(ld, current.get_ngk(ik), nband, current.get_pointer(), hpsi, expectations.data());
        for (int band = 0; band < nband; ++band)
        {
            (*energies)(ik, band) = std::real(expectations[band]);
        }
    }
    ModuleBase::timer::end("HSolverPWTDDFT", "cal_band_energy");
}

template <typename T, typename Device>
void HSolverPWTDDFT<T, Device>::reset_orth_stats()
{
    if (options_.out_stat)
    {
        orth_stats_ = TDOrthStats();
    }
    warned_events_ = 0;
}

template <typename T, typename Device>
void HSolverPWTDDFT<T, Device>::correct_orbitals(T* current, int ld, int dim, int bands, int ik, int istep, int iter, bool detailed_output)
{
    ModuleBase::timer::start("HSolverPWTDDFT", "correct_orbitals");
    const OrthResult result = orthonormal_.apply(current, ld, dim, bands, options_.orthonormal, options_.out_stat);
    record_orth(result, ik, istep, iter);
    if (detailed_output && result.gram_checked && log_.good())
    {
        const int evolution_step = istep + 1;
        std::ostringstream record;
        record << std::scientific << std::setprecision(6) << " PW RT-TDDFT orth: evolution_step=" << evolution_step << " scf_iter=" << iter
               << " k=" << ik << " requested=" << orth_method_name(options_.orthonormal) << " orth_before=" << result.before
               << " orth_after=" << result.after << '\n';
        log_ << record.str();
    }
    ModuleBase::timer::end("HSolverPWTDDFT", "correct_orbitals");
}

template <typename T, typename Device>
void HSolverPWTDDFT<T, Device>::record_orth(const OrthResult& result, int ik, int istep, int iter)
{
    if (result.status == OrthStatus::failed)
    {
        const char* stage = istep == 0 ? "initial check" : "propagation";
        std::ostringstream message;
        const int evolution_step = istep + 1;
        message << std::setprecision(16) << " PW RT-TDDFT orbital validation failed (" << stage << "): evolution_step=" << evolution_step
                << " scf_iter=" << iter << " local_k=" << ik << " requested=" << orth_method_name(options_.orthonormal) << '\n'
                << "   last_attempted=" << orth_method_name(result.actual) << " passes=" << result.passes
                << " fallbacks=" << result.fallbacks << '\n'
                << "   " << orth_failure_name(result.failure) << "; " << result.reason;
        if (result.gram_checked)
        {
            message << '\n' << "   before=" << result.before << " after=" << result.after << " tolerance=" << orth_tolerance<T>();
        }
        // stdout is disabled on non-world-root ranks; a failing pool must still report its error.
        if (comm_.rank == 0)
        {
            std::cerr << message.str() << std::endl;
        }
#ifdef __MPI
        // Let the pool root flush its diagnostic before another rank's exit stops the MPI job.
        double diagnostic_written = 1.0;
        Parallel_Common::bcast_data(&diagnostic_written, 1, comm_.comm, 0);
#endif
        ModuleBase::WARNING_QUIT("HSolverPWTDDFT", message.str());
    }
    report_orth_warning(result, ik, istep, iter);
    if (options_.out_stat)
    {
        orth_norms_[ik] = result.norms;
        if (result.gram_checked)
        {
            orth_stats_.before = std::max(orth_stats_.before, result.before);
            orth_stats_.after = std::max(orth_stats_.after, result.after);
        }
    }
}

template <typename T, typename Device>
void HSolverPWTDDFT<T, Device>::report_orth_warning(const OrthResult& result, int ik, int istep, int iter)
{
    if (comm_.rank != 0)
    {
        return;
    }
    const bool initial_warning = istep == 0 && result.gram_checked && result.after > orth_tolerance<T>();
    unsigned int events = 0;
    if (result.fallbacks > 0)
    {
        events |= 1;
    }
    if (result.rejected > 0)
    {
        events |= 2;
    }
    const unsigned int fresh_events = events & ~warned_events_;
    if (!initial_warning && fresh_events == 0)
    {
        return;
    }
    warned_events_ |= events;
    const int evolution_step = istep + 1;
    std::ostringstream message;
    message << std::setprecision(16) << " PW RT-TDDFT orbital warning: evolution_step=" << evolution_step << " scf_iter=" << iter;
    if (options_.global_k_indices.empty())
    {
        message << " local_k=" << ik;
    }
    else
    {
        message << " global_k=" << options_.global_k_indices.at(ik);
    }
    message << " requested=" << orth_method_name(options_.orthonormal) << '\n'
            << "   before=" << result.before << " after=" << result.after << " tolerance=" << orth_tolerance<T>() << '\n'
            << "   ";
    if (initial_warning)
    {
        message << "Initial orbitals are unchanged; finite orthogonality error exceeds the propagation correction tolerance.";
    }
    else
    {
        message << "last_attempted=" << orth_method_name(result.actual) << ". ";
        if (fresh_events & 1)
        {
            message << "A fallback was attempted. ";
        }
        if (fresh_events & 2)
        {
            message << "A candidate was rejected. ";
        }
        message << "The retained state satisfies the tolerance." << '\n'
                << "   " << result.reason << '\n'
                << "   Further events of these types in this pool are suppressed for this electronic evolution step.";
    }
    const std::string text = message.str();
    std::cerr << text << std::endl;
    if (log_.good())
    {
        log_ << text << std::endl;
    }
}

template <typename T, typename Device>
void HSolverPWTDDFT<T, Device>::check_initial(const psi::Psi<T, Device>& current, int iter)
{
    ModuleBase::timer::start("HSolverPWTDDFT", "check_initial");
    const bool full_gram = options_.orthonormal != OrthMethod::none;
    if (options_.out_stat)
    {
        orth_norms_.resize(current.get_nk());
    }
    for (int ik = 0; ik < current.get_nk(); ++ik)
    {
        current.fix_k(ik);
        const OrthResult result = orthonormal_.inspect(current.get_pointer(),
                                                       current.get_nbasis(),
                                                       current.get_ngk(ik),
                                                       current.get_nbands(),
                                                       full_gram,
                                                       options_.out_stat);
        record_orth(result, ik, 0, iter);
    }
    ModuleBase::timer::end("HSolverPWTDDFT", "check_initial");
}

template <typename T, typename Device>
double HSolverPWTDDFT<T, Device>::wave_electrons(const ModuleBase::matrix& occupations) const
{
    double electrons = 0.0;
    for (std::size_t ik = 0; ik < orth_norms_.size(); ++ik)
    {
        for (std::size_t band = 0; band < orth_norms_[ik].size(); ++band)
        {
            electrons += occupations(ik, band) * orth_norms_[ik][band];
        }
    }
    return electrons;
}

template class HSolverPWTDDFT<std::complex<float>, base_device::DEVICE_CPU>;
template class HSolverPWTDDFT<std::complex<double>, base_device::DEVICE_CPU>;
#if defined(__CUDA) || defined(__ROCM)
template class HSolverPWTDDFT<std::complex<float>, base_device::DEVICE_GPU>;
template class HSolverPWTDDFT<std::complex<double>, base_device::DEVICE_GPU>;
#endif
} // namespace hsolver
