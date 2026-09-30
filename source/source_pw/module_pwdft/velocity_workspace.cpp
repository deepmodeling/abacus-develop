#include "source_pw/module_pwdft/velocity_workspace.h"

#include "source_base/kernels/math_kernel_op.h"
#include "source_base/module_device/memory_op.h"
#include "source_base/parallel_device.h"
#include "source_base/timer.h"
#include "source_pw/module_pwdft/kernels/nonlocal_op.h"

#include <algorithm>

namespace hamilt
{
template <typename Real, typename Device>
std::complex<Real>* VelocityWorkspace<Real, Device>::prepare(const int count)
{
    using CtDevice = typename ct::PsiToContainer<Device>::type;
    if (count > capacity_ || coefficients_.data_type() != ct::DataTypeToEnum<Complex>::value
        || coefficients_.device_type() != ct::DeviceTypeToEnum<CtDevice>::value)
    {
        coefficients_
            = ct::Tensor(ct::DataTypeToEnum<Complex>::value, ct::DeviceTypeToEnum<CtDevice>::value, {std::max<int64_t>(1, 8LL * count)});
        host_.resize(4 * count);
        capacity_ = count;
    }
    return coefficients_.template data<Complex>();
}

template <typename Real, typename Device>
void VelocityWorkspace<Real, Device>::contract(const UnitCell& cell,
                                               const pseudopot_cell_vnl& pp,
                                               const int spin,
                                               const int bands,
                                               const Real scale,
                                               const ModulePW::PW_Basis_K& basis,
                                               Complex* buffer)
{
    ModuleBase::timer::start("VelocityWorkspace", "contract");
    const int count = bands * pp.nkb;
    if (count == 0)
    {
        ModuleBase::timer::end("VelocityWorkspace", "contract");
        return;
    }
#ifdef __MPI
    if (basis.poolnproc > 1)
    {
        Parallel_Common::reduce_dev<Complex, Device>(buffer, 4 * count, basis.pool_world, host_.data());
    }
#endif
    Complex* output = buffer + 4 * count;
    base_device::memory::set_memory_op<Complex, Device>()(output, 0, 4 * count);
    for (int component = 0; component < 4; ++component)
    {
        int sum = 0;
        int atom = 0;
        const int stride = component == 0 ? pp.nkb : 3 * pp.nkb;
        const Complex* input = component == 0 ? buffer : buffer + count + (component - 1) * pp.nkb;
        for (int type = 0; type < cell.ntype; ++type)
        {
            nonlocal_pw_op<Real, Device>()(nullptr,
                                           cell.atoms[type].na,
                                           bands,
                                           cell.atoms[type].ncpp.nh,
                                           sum,
                                           atom,
                                           spin,
                                           stride,
                                           pp.deeq.getBound2(),
                                           pp.deeq.getBound3(),
                                           pp.deeq.getBound4(),
                                           pp.template get_deeq_data<Real>(),
                                           output + component * count,
                                           input);
        }
    }
    if (scale != Real(1))
    {
        const Complex factor(scale, 0);
        ModuleBase::scal_op<Real, Device>()(4 * count, &factor, output, 1);
    }
    ModuleBase::timer::end("VelocityWorkspace", "contract");
}
template class VelocityWorkspace<float, base_device::DEVICE_CPU>;
template class VelocityWorkspace<double, base_device::DEVICE_CPU>;
#if defined(__CUDA) || defined(__ROCM)
template class VelocityWorkspace<float, base_device::DEVICE_GPU>;
template class VelocityWorkspace<double, base_device::DEVICE_GPU>;
#endif
} // namespace hamilt
