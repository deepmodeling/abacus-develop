#include "source_hsolver/hsolver_pw_tddft.h"

#include "source_base/module_device/memory_op.h"
#include "source_base/parallel_device.h"
#include "source_base/timer.h"
#include "source_base/tool_quit.h"
#include "source_basis/module_pw/pw_momentum.h"
#include "source_hsolver/kernels/linear_op.h"

#include <chrono>
#include <functional>
#include <iomanip>
#include <sstream>

namespace hsolver
{
namespace
{

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
    if (name == "bicgstab")
    {
        return LinearMethod::bicgstab;
    }
    if (name == "cgs")
    {
        return LinearMethod::cgs;
    }
    if (name == "gmres")
    {
        return LinearMethod::gmres;
    }
    ModuleBase::WARNING_QUIT("HSolverPWTDDFT", "Unsupported linear solver: " + name);
    return LinearMethod::bicgstab;
}

PWPreconditioner parse_pw_precond(const std::string& name)
{
    if (name == "none")
    {
        return PWPreconditioner::none;
    }
    if (name == "kinetic")
    {
        return PWPreconditioner::kinetic;
    }
    if (name == "kinetic_recycle")
    {
        return PWPreconditioner::kinetic_recycle;
    }
    if (name == "kinetic_subspace")
    {
        return PWPreconditioner::kinetic_subspace;
    }
    ModuleBase::WARNING_QUIT("HSolverPWTDDFT", "Unsupported preconditioner: " + name);
    return PWPreconditioner::kinetic;
}

template <typename T, typename Device>
HSolverPWTDDFT<T, Device>::HSolverPWTDDFT(const ModulePW::PW_Basis_K& basis,
                                          const PWLinearOptions& options,
                                          const diag_comm_info& comm,
                                          std::ostream& log)
    : basis_(basis), comm_(comm), options_(options), algebra_(comm), band_products_(comm)
{
    initialize(log);
}

template <typename T, typename Device>
void HSolverPWTDDFT<T, Device>::initialize(std::ostream& log)
{
    const char* method = nullptr;
    switch (options_.linear.method)
    {
    case LinearMethod::bicgstab:
        method = "bicgstab";
        break;
    case LinearMethod::cgs:
        method = "cgs";
        break;
    case LinearMethod::gmres:
        method = "gmres";
        break;
    default:
        ModuleBase::WARNING_QUIT("HSolverPWTDDFT", "Unsupported linear solver.");
        return;
    }
    const char* preconditioner = nullptr;
    switch (options_.preconditioner)
    {
    case PWPreconditioner::none:
        preconditioner = "none";
        break;
    case PWPreconditioner::kinetic:
        preconditioner = "kinetic";
        break;
    case PWPreconditioner::kinetic_recycle:
        preconditioner = "kinetic_recycle";
        break;
    case PWPreconditioner::kinetic_subspace:
        preconditioner = "kinetic_subspace";
        break;
    default:
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
    info << "RT-TDDFT linear solver: " << method << "; preconditioner: " << preconditioner << "; tolerance: " << std::setprecision(16)
         << linear_solver_->tolerance() << "; maximum iterations: " << options_.linear.max_iterations
         << "; restart: " << options_.linear.restart << "; CN initial guess: " << options_.cn_init
         << "; residual reconstruction: " << options_.linear.reconstruct << '\n';
    log << info.str();
}

template <typename T, typename Device>
void HSolverPWTDDFT<T, Device>::prepare_buffers(const int nbands, const int nbasis)
{
    ModuleBase::timer::start("HSolverPWTDDFT", "prepare_buffers");
    using CtDevice = typename ct::PsiToContainer<Device>::type;
    const ct::DeviceType device = ct::DeviceTypeToEnum<CtDevice>::value;
    const int64_t size = std::max<int64_t>(1, static_cast<int64_t>(nbands) * nbasis);
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
    if (inverse_kinetic_.NumElements() < std::max(1, dim) || inverse_kinetic_.data_type() != ct::DataTypeToEnum<T>::value
        || inverse_kinetic_.device_type() != device)
    {
        inverse_kinetic_ = ct::Tensor(ct::DataTypeToEnum<T>::value, device, {std::max(1, dim)});
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
void HSolverPWTDDFT<T, Device>::prepare_kpoint(int ik, int dim, KPointState* state)
{
    if (!state)
    {
        return;
    }
    std::size_t signature = static_cast<std::size_t>(dim);
    const auto hash_value
        = [&signature](double value) { signature ^= std::hash<double>()(value) + 0x9e3779b9 + (signature << 6) + (signature >> 2); };
    hash_value(basis_.tpiba);
    if (basis_.kvec_c)
    {
        hash_value(basis_.kvec_c[ik].x);
        hash_value(basis_.kvec_c[ik].y);
        hash_value(basis_.kvec_c[ik].z);
    }
    for (int ig = 0; ig < dim && basis_.gcar; ++ig)
    {
        const ModuleBase::Vector3<double>& g = basis_.getgcar(ik, ig);
        hash_value(g.x);
        hash_value(g.y);
        hash_value(g.z);
    }
    double changed = !state->layout_valid || signature != state->layout ? 1.0 : 0.0;
#ifdef __MPI
    if (comm_.nproc > 1)
    {
        Parallel_Common::reduce_data(&changed, 1, comm_.comm);
    }
#endif
    if (changed > 0)
    {
        state->response.clear();
        state->solve_count = 0;
        state->independent_step = -1;
    }
    state->layout = signature;
    state->layout_valid = true;
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
    LinearSolveOptions fallback = options_.linear;
    fallback.max_iterations -= result->iterations;
    fallback.reconstruct = false;
    HSolverLinear<T, Device> solver(fallback, comm_);
    const T* inverse = inverse_kinetic_.template data<T>();
    const LinearLowRank<T, Device> diagonal(algebra_, inverse, batch.dim);
    const T* rhs = rhs_.template data<T>();
    LinearSolveResult retry = solver.solve(op, diagonal, batch.ld, batch.bands, batch.dim, current, rhs);
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
        preconditioner.prepare_subspace(batch.ld, batch.bands, previous, image, factor);
    }
    else if (options_.preconditioner == PWPreconditioner::kinetic_recycle && state->response.rank() > 0)
    {
        const int rank = state->response.rank();
        const T* directions = state->response.directions();
        const T* images = state->response.images();
        preconditioner.prepare_response(batch.ld, rank, directions, images);
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
void HSolverPWTDDFT<T, Device>::report_solve(const SolveDetails& details,
                                             int ik,
                                             int step,
                                             int iteration,
                                             double elapsed,
                                             std::ostream& log) const
{
    std::vector<double> times(comm_.nproc, 0.0);
    times[comm_.rank] = elapsed;
#ifdef __MPI
    if (comm_.nproc > 1)
    {
        Parallel_Common::reduce_data(times.data(), times.size(), comm_.comm);
    }
#endif
    const double seconds = *std::max_element(times.begin(), times.end());
    const LinearSolveResult& result = details.linear;
    std::ostringstream record;
    record << std::setprecision(10) << "Linear solve: step=" << step << " iter=" << iteration << " k=" << ik
           << " iterations=" << result.iterations << " restarts=" << result.restarts << " operator_calls=" << result.operator_calls
           << " operator_columns=" << result.operator_columns << " residual=" << result.max_residual
           << " residual_kind=" << (result.reconstructed ? "reconstructed" : "independent") << " true_checks=" << result.true_checks
           << " reconstruction_fallbacks=" << result.reconstruction_fallbacks << " coarse_rank=" << details.coarse_rank
           << " kinetic_retry=" << details.retried << " cn_projected=" << details.projected << " cn_initial=" << details.cn_initial
           << " seconds=" << seconds << '\n';
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
        if (detailed_output)
        {
            linear_op<T, Device>().synchronize();
        }
        const std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
        KPointState* state = tracks_state() ? &states_[ik] : nullptr;
        prepare_kpoint(ik, dim, state);
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
        if (detailed_output)
        {
            linear_op<T, Device>().synchronize();
            const double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
            report_solve(details, ik, istep, iter, elapsed, log);
        }
        if (details.linear.status != LinearSolveStatus::converged)
        {
            const LinearSolveResult& result = details.linear;
            std::ostringstream message;
            message << "Linear solve failed at step " << istep << ", k point " << ik << ", band " << result.failed_band << ", after "
                    << result.iterations << " iterations: " << linear_status_name(result.status) << "; residual = " << result.max_residual;
            ModuleBase::WARNING_QUIT("HSolverPWTDDFT", message.str());
        }
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

template class HSolverPWTDDFT<std::complex<float>, base_device::DEVICE_CPU>;
template class HSolverPWTDDFT<std::complex<double>, base_device::DEVICE_CPU>;
#if defined(__CUDA) || defined(__ROCM)
template class HSolverPWTDDFT<std::complex<float>, base_device::DEVICE_GPU>;
template class HSolverPWTDDFT<std::complex<double>, base_device::DEVICE_GPU>;
#endif
} // namespace hsolver
