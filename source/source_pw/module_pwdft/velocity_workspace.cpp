#include "source_pw/module_pwdft/velocity_workspace.h"

#include "source_base/kernels/math_kernel_op.h"
#include "source_base/module_device/memory_op.h"
#include "source_base/parallel_device.h"
#include "source_base/timer.h"
#include "source_pw/module_pwdft/kernels/nonlocal_op.h"

#include <algorithm>
#include <limits>

namespace hamilt
{
template <typename Real, typename Device>
std::complex<Real>* VelocityWorkspace<Real, Device>::prepare(const std::int64_t count)
{
    using CtDevice = typename ct::PsiToContainer<Device>::type;
    if (count > capacity_ || coefficients_.data_type() != ct::DataTypeToEnum<Complex>::value
        || coefficients_.device_type() != ct::DeviceTypeToEnum<CtDevice>::value)
    {
        const std::int64_t workspace_elements = 8 * count;
        const std::int64_t capacity = std::max<std::int64_t>(1, workspace_elements);
        const std::int64_t projection_elements = 4 * count;
        coefficients_ = ct::Tensor(ct::DataTypeToEnum<Complex>::value, ct::DeviceTypeToEnum<CtDevice>::value, {capacity});
        host_.resize(projection_elements);
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
    const std::int64_t count = static_cast<std::int64_t>(bands) * pp.nkb;
    if (count == 0)
    {
        ModuleBase::timer::end("VelocityWorkspace", "contract");
        return;
    }
    const std::int64_t projection_elements = 4 * count;
#ifdef __MPI
    if (basis.poolnproc > 1)
    {
        // The complex reduction wrapper doubles the MPI element count.
        const std::int64_t max_chunk = std::numeric_limits<int>::max() / 2;
        for (std::int64_t offset = 0; offset < projection_elements;)
        {
            const std::int64_t remaining = projection_elements - offset;
            const int chunk_size = static_cast<int>(std::min(max_chunk, remaining));
            Complex* chunk = buffer + offset;
            Complex* host_chunk = host_.data() + offset;
            Parallel_Common::reduce_dev<Complex, Device>(chunk, chunk_size, basis.pool_world, host_chunk);
            offset += chunk_size;
        }
    }
#endif
    Complex* output = buffer + projection_elements;
    base_device::memory::set_memory_op<Complex, Device>()(output, 0, projection_elements);
    for (int component = 0; component < 4; ++component)
    {
        int sum = 0;
        int atom = 0;
        const int stride = component == 0 ? pp.nkb : 3 * pp.nkb;
        const std::int64_t gradient_offset = static_cast<std::int64_t>(component - 1) * pp.nkb;
        const Complex* input = component == 0 ? buffer : buffer + count + gradient_offset;
        Complex* component_output = output + component * count;
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
                                           component_output,
                                           input);
        }
    }
    if (scale != Real(1))
    {
        const Complex factor(scale, 0);
        const std::int64_t max_chunk = std::numeric_limits<int>::max();
        for (std::int64_t offset = 0; offset < projection_elements;)
        {
            const std::int64_t remaining = projection_elements - offset;
            const int chunk_size = static_cast<int>(std::min(max_chunk, remaining));
            Complex* chunk = output + offset;
            ModuleBase::scal_op<Real, Device>()(chunk_size, &factor, chunk, 1);
            offset += chunk_size;
        }
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
