#include "source_hsolver/module_rt_pw/exp_action.h"

#include "source_hsolver/module_rt_pw/exp_lanczos.h"
#include "source_hsolver/module_rt_pw/exp_taylor.h"

#include <limits>
#include <stdexcept>

namespace hsolver
{
template <typename T, typename Device>
ExpAction<T, Device>::ExpAction(const ExpOpts& exp_opts, const diag_comm_info& comm) : exp_opts_(exp_opts), algebra_(comm), work_(comm, 1)
{
    using Real = typename T::value_type;
    if (exp_opts_.rel_tol == 0.0)
    {
        exp_opts_.rel_tol = std::max(1e-10, 100.0 * std::numeric_limits<Real>::epsilon());
    }
    if (exp_opts_.taylor_order < 2 || exp_opts_.krylov_max < 2 || !std::isfinite(exp_opts_.rel_tol) || exp_opts_.rel_tol <= 0)
    {
        throw std::invalid_argument("Invalid exponential controls.");
    }
}

template <typename T, typename Device>
std::unique_ptr<ExpAction<T, Device>> make_exp_action(const ExpOpts& exp_opts, const diag_comm_info& comm)
{
    if (exp_opts.method == ExpMethod::taylor)
    {
        return std::unique_ptr<ExpAction<T, Device>>(new ExpTaylor<T, Device>(exp_opts, comm));
    }
    return std::unique_ptr<ExpAction<T, Device>>(new ExpLanczos<T, Device>(exp_opts, comm));
}

template class ExpAction<std::complex<float>, base_device::DEVICE_CPU>;
template class ExpAction<std::complex<double>, base_device::DEVICE_CPU>;
template std::unique_ptr<ExpAction<std::complex<float>, base_device::DEVICE_CPU>> make_exp_action(const ExpOpts&, const diag_comm_info&);
template std::unique_ptr<ExpAction<std::complex<double>, base_device::DEVICE_CPU>> make_exp_action(const ExpOpts&, const diag_comm_info&);
#if defined(__CUDA) || defined(__ROCM)
template class ExpAction<std::complex<float>, base_device::DEVICE_GPU>;
template class ExpAction<std::complex<double>, base_device::DEVICE_GPU>;
template std::unique_ptr<ExpAction<std::complex<float>, base_device::DEVICE_GPU>> make_exp_action(const ExpOpts&, const diag_comm_info&);
template std::unique_ptr<ExpAction<std::complex<double>, base_device::DEVICE_GPU>> make_exp_action(const ExpOpts&, const diag_comm_info&);
#endif
} // namespace hsolver
