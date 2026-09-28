#include "source_pw/module_pwdft/projector_gradient.h"

#include "source_base/timer.h"
#include "source_pw/module_pwdft/kernels/projector_gradient_op.h"

#include <algorithm>
#include <limits>

namespace hamilt
{
template <typename Real, typename Device>
template <typename Value>
void ProjectorGradient<Real, Device>::reserve(ct::Tensor* tensor, const int64_t count)
{
    using CtDevice = typename ct::PsiToContainer<Device>::type;
    const int64_t capacity = std::max<int64_t>(1, count);
    if (tensor->NumElements() < capacity || tensor->data_type() != ct::DataTypeToEnum<Value>::value
        || tensor->device_type() != ct::DeviceTypeToEnum<CtDevice>::value)
    {
        *tensor = ct::Tensor(ct::DataTypeToEnum<Value>::value, ct::DeviceTypeToEnum<CtDevice>::value, {capacity});
    }
}

template <typename Real, typename Device>
void ProjectorGradient<Real, Device>::calculate(pseudopot_cell_vnl* pp,
                                                const UnitCell& cell,
                                                const ModulePW::PW_Basis_K& basis,
                                                const int ik,
                                                const ModuleBase::Vector3<double>& A,
                                                Complex* output)
{
    ModuleBase::timer::start("ProjectorGradient", "calculate");
    if (pp->nkb == 0)
    {
        ModuleBase::timer::end("ProjectorGradient", "calculate");
        return;
    }
    pp->ensure_grad_table(cell);
    const double dq = pp->table_dq();
    if (version_ != pp->table_version())
    {
        reserve<Real>(&radial_, pp->tab.getSize());
        reserve<Real>(&derivative_, pp->tab_dq.getSize());
        base_device::memory::cast_memory_op<Real, double, Device, base_device::DEVICE_CPU>()(radial_.template data<Real>(),
                                                                                             pp->tab.ptr,
                                                                                             pp->tab.getSize());
        base_device::memory::cast_memory_op<Real, double, Device, base_device::DEVICE_CPU>()(derivative_.template data<Real>(),
                                                                                             pp->tab_dq.ptr,
                                                                                             pp->tab_dq.getSize());
        host_metadata_.clear();
        int atom = 0;
        for (int type = 0; type < cell.ntype; ++type)
        {
            for (int ia = 0; ia < cell.atoms[type].na; ++ia, ++atom)
            {
                for (int p = 0; p < cell.atoms[type].ncpp.nh; ++p)
                {
                    host_metadata_.push_back(type);
                    host_metadata_.push_back(static_cast<int>(pp->indv(type, p)));
                    host_metadata_.push_back(static_cast<int>(pp->nhtol(type, p)));
                    host_metadata_.push_back(static_cast<int>(pp->nhtolm(type, p)));
                    host_metadata_.push_back(atom);
                }
            }
        }
        reserve<int>(&metadata_, host_metadata_.size());
        base_device::memory::synchronize_memory_op<int, Device, base_device::DEVICE_CPU>()(metadata_.template data<int>(),
                                                                                           host_metadata_.data(),
                                                                                           host_metadata_.size());
        version_ = pp->table_version();
    }
    const int npw = basis.npwk[ik];
    if (npw == 0)
    {
        ModuleBase::timer::end("ProjectorGradient", "calculate");
        return;
    }
    host_q_.resize(3 * npw);
    for (int ig = 0; ig < npw; ++ig)
    {
        const ModuleBase::Vector3<double> q = basis.getgpluskcar(ik, ig) + A / cell.tpiba;
        for (int d = 0; d < 3; ++d)
        {
            host_q_[3 * ig + d] = static_cast<Real>(q[d]);
        }
        const Real* v = host_q_.data() + 3 * ig;
        const Real position = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) * static_cast<Real>(cell.tpiba) / static_cast<Real>(dq);
        pp->check_vnl_index(static_cast<double>(position) * (1.0 + 8.0 * std::numeric_limits<Real>::epsilon()), true);
    }
    reserve<Real>(&momentum_, 3LL * npw);
    reserve<Complex>(&structure_, static_cast<int64_t>(cell.nat) * npw);
    base_device::memory::synchronize_memory_op<Real, Device, base_device::DEVICE_CPU>()(momentum_.template data<Real>(),
                                                                                        host_q_.data(),
                                                                                        3 * npw);
    Device* ctx = nullptr;
    pp->psf->get_sk(ctx, ik, &basis, structure_.template data<Complex>());
    projector_gradient_op<Real, Device>()(npw,
                                          basis.npwk_max,
                                          pp->nkb,
                                          pp->tab.getBound2(),
                                          pp->tab.getBound3(),
                                          static_cast<Real>(dq),
                                          static_cast<Real>(cell.tpiba),
                                          metadata_.template data<int>(),
                                          momentum_.template data<Real>(),
                                          radial_.template data<Real>(),
                                          derivative_.template data<Real>(),
                                          structure_.template data<Complex>(),
                                          output);
    ModuleBase::timer::end("ProjectorGradient", "calculate");
}
template class ProjectorGradient<float, base_device::DEVICE_CPU>;
template class ProjectorGradient<double, base_device::DEVICE_CPU>;
#if defined(__CUDA) || defined(__ROCM)
template class ProjectorGradient<float, base_device::DEVICE_GPU>;
template class ProjectorGradient<double, base_device::DEVICE_GPU>;
#endif
} // namespace hamilt
