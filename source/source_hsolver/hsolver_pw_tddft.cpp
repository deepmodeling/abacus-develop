#include "source_hsolver/hsolver_pw_tddft.h"

#include "source_base/module_device/memory_op.h"
#include "source_base/timer.h"
#include "source_base/tool_quit.h"
#include "source_basis/module_pw/pw_momentum.h"
#include "source_hsolver/kernels/linear_op.h"

#include <iomanip>
#include <sstream>

namespace hsolver
{
namespace
{

template <typename T, typename Device>
class DiagonalPreconditioner final : public LinearOperator<T, Device>
{
  public:
    DiagonalPreconditioner(const T* inverse, const int dim) : inverse_(inverse), dim_(dim)
    {
    }
    void apply(const T* x, T* y, const int ld, const int nvec) const override
    {
        linear_op<T, Device>().diagonal(ld, dim_, nvec, inverse_, x, y);
    }

  private:
    const T* inverse_;
    const int dim_;
};

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

template <typename T, typename Device>
HSolverPWTDDFT<T, Device>::HSolverPWTDDFT(const ModulePW::PW_Basis_K& basis,
                                          const std::string& method,
                                          const std::string& preconditioner,
                                          const double tolerance,
                                          const int max_iterations,
                                          const bool kinetic_enabled,
                                          const diag_comm_info& comm,
                                          std::ostream& log)
    : basis_(basis), kinetic_enabled_(kinetic_enabled), kinetic_preconditioner_(preconditioner == "kinetic"), band_products_(comm)
{
    LinearSolveOptions options;
    if (method == "bicgstab")
    {
        options.method = LinearMethod::bicgstab;
    }
    else if (method == "cgs")
    {
        options.method = LinearMethod::cgs;
    }
    else
    {
        ModuleBase::WARNING_QUIT("HSolverPWTDDFT", "Unsupported linear solver: " + method);
    }
    if (preconditioner != "kinetic" && preconditioner != "none")
    {
        ModuleBase::WARNING_QUIT("HSolverPWTDDFT", "Unsupported preconditioner: " + preconditioner);
    }
    options.tolerance = tolerance;
    options.max_iterations = max_iterations;
    linear_solver_.reset(new HSolverLinear<T, Device>(options, comm));
    std::ostringstream info;
    info << "RT-TDDFT linear solver: " << method << "; preconditioner: " << preconditioner << "; tolerance: " << std::setprecision(16)
         << linear_solver_->tolerance() << "; maximum iterations: " << max_iterations;
    log << info.str() << std::endl;
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
        const double kinetic = kinetic_enabled_ ? ModulePW::shifted_kinetic(basis_, ik, ig, momentum_shift) : 0.0;
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
    const int nband = current->get_nbands();
    const int ld = current->get_nbasis();
    prepare_buffers(nband, ld);
    T* rhs = rhs_.template data<T>();
    // H is stored in Rydberg; CN includes the conversion to Hartree.
    const T coefficient(0.0, dt / 4.0);
    for (int ik = 0; ik < current->get_nk(); ++ik)
    {
        op.update_k(ik);
        current->fix_k(ik);
        previous.fix_k(ik);
        const int dim = current->get_ngk(ik);
        const ShiftedHOperator<T, Device> rhs_op(op, -coefficient, dim);
        const ShiftedHOperator<T, Device> lhs_op(op, coefficient, dim);
        rhs_op.apply(previous.get_pointer(), rhs, ld, nband);
        LinearSolveResult result;
        if (kinetic_preconditioner_)
        {
            update_precond(ik, dim, coefficient, momentum_shift);
            const DiagonalPreconditioner<T, Device> preconditioner(inverse_kinetic_.template data<T>(), dim);
            result = linear_solver_->solve(lhs_op, preconditioner, ld, nband, dim, current->get_pointer(), rhs);
        }
        else
        {
            result = linear_solver_->solve(lhs_op, ld, nband, dim, current->get_pointer(), rhs);
        }
        if (detailed_output)
        {
            log << "Linear solve: step=" << istep << " iter=" << iter << " k=" << ik << " iterations=" << result.iterations
                << " restarts=" << result.restarts << " operator_calls=" << result.operator_calls
                << " operator_columns=" << result.operator_columns << " residual=" << result.max_residual << '\n';
        }
        if (result.status != LinearSolveStatus::converged)
        {
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
