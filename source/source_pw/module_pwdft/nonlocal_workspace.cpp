#include "source_pw/module_pwdft/nonlocal_workspace.h"

#include "source_base/kernels/math_kernel_op.h"
#include "source_base/parallel_device.h"
#include "source_base/timer.h"
#include "source_pw/module_pwdft/kernels/nonlocal_op.h"

namespace hamilt
{
namespace
{
#ifdef __DSP
template <typename T, typename Device>
using Resize = base_device::memory::resize_memory_op_mt<T, Device>;
template <typename T, typename Device>
using Delete = base_device::memory::delete_memory_op_mt<T, Device>;
template <typename T, typename Device>
using Zero = base_device::memory::set_memory_op_mt<T, Device>;
template <typename T, typename Device>
using Gemm = ModuleBase::gemm_op_mt<T, Device>;
#else
template <typename T, typename Device>
using Resize = base_device::memory::resize_memory_op<T, Device>;
template <typename T, typename Device>
using Delete = base_device::memory::delete_memory_op<T, Device>;
template <typename T, typename Device>
using Zero = base_device::memory::set_memory_op<T, Device>;
template <typename T, typename Device>
using Gemm = ModuleBase::gemm_op<T, Device>;
#endif
}

template <typename T, typename Device>
NonlocalWorkspace<T, Device>::~NonlocalWorkspace()
{
    Delete<T, Device>()(becp_);
    Delete<T, Device>()(ps_);
}

template <typename T, typename Device>
void NonlocalWorkspace<T, Device>::project(const T* vkb, const T* psi, const int npw,
                                         const int ldv, const int ldp, const int nkb, const int bands)
{
    const T one(1, 0);
    const T zero(0, 0);
    if (npw == 0)
    {
        Zero<T, Device>()(becp_, 0, nkb * bands);
    }
    else if (bands == 1)
    {
        ModuleBase::gemv_op<T, Device>()('C', npw, nkb, &one, vkb, ldv, psi, 1, &zero, becp_, 1);
    }
    else
    {
        Gemm<T, Device>()('C', 'N', nkb, bands, npw, &one, vkb, ldv, psi, ldp, &zero, becp_, nkb);
    }
}

template <typename T, typename Device>
void NonlocalWorkspace<T, Device>::contract(const UnitCell& cell, const pseudopot_cell_vnl& pp,
                                          const int spin, const int npol, const int bands)
{
    using Real = typename GetTypeReal<T>::type;
    Zero<T, Device>()(ps_, 0, pp.nkb * bands);
    int sum = 0;
    int atom = 0;
    for (int type = 0; type < cell.ntype; ++type)
    {
        const int nproj = cell.atoms[type].ncpp.nh;
        if (npol == 1)
        {
            nonlocal_pw_op<Real, Device>()(nullptr, cell.atoms[type].na, bands, nproj, sum, atom,
                spin, pp.nkb, pp.deeq.getBound2(), pp.deeq.getBound3(), pp.deeq.getBound4(),
                pp.template get_deeq_data<Real>(), ps_, becp_);
        }
        else
        {
            nonlocal_pw_op<Real, Device>()(nullptr, cell.atoms[type].na, bands, nproj, sum, atom,
                pp.nkb, pp.deeq_nc.getBound2(), pp.deeq_nc.getBound3(), pp.deeq_nc.getBound4(),
                pp.template get_deeq_nc_data<Real>(), ps_, becp_);
        }
    }
}

template <typename T, typename Device>
void NonlocalWorkspace<T, Device>::back_project(const T* vkb, T* hpsi, const int npw,
                                              const int ldv, const int ldp, const int nkb, const int bands)
{
    const T one(1, 0);
    if (bands == 1)
    {
        ModuleBase::gemv_op<T, Device>()('N', npw, nkb, &one, vkb, ldv, ps_, 1, &one, hpsi, 1);
    }
    else
    {
        Gemm<T, Device>()('N', 'T', npw, bands, nkb, &one, vkb, ldv, ps_, bands, &one, hpsi, ldp);
    }
}

template <typename T, typename Device>
void NonlocalWorkspace<T, Device>::apply(const UnitCell& cell, const pseudopot_cell_vnl& pp,
                                       const ModulePW::PW_Basis_K& basis, const int spin,
                                       const int bands, const int nbasis, const int npol, const int npw,
                                       const bool is_first_node, const T* vkb, const T* psi, T* hpsi)
{
    ModuleBase::timer::start("NonlocalWorkspace", "apply");
    if (is_first_node)
    {
        Zero<T, Device>()(hpsi, 0, nbasis * bands / npol);
    }
    if (pp.nkb > 0 && bands > 0)
    {
        const size_t count = static_cast<size_t>(pp.nkb) * bands;
        if (count > becp_capacity_)
        {
            Resize<T, Device>()(becp_, count, "Nonlocal::becp");
            becp_capacity_ = count;
        }
        if (count > ps_capacity_)
        {
            Resize<T, Device>()(ps_, count, "Nonlocal::ps");
            ps_capacity_ = count;
        }
        this->project(vkb, psi, npw, pp.vkbnc, nbasis / npol, pp.nkb, bands);
#ifdef __MPI
        if (basis.poolnproc > 1)
        {
            Parallel_Common::reduce_dev<T, Device>(becp_, pp.nkb * bands, basis.pool_world);
        }
#endif
        if (npw > 0)
        {
            this->contract(cell, pp, spin, npol, bands);
            this->back_project(vkb, hpsi, npw, pp.vkbnc, nbasis / npol, pp.nkb, bands);
        }
    }
    ModuleBase::timer::end("NonlocalWorkspace", "apply");
}

template class NonlocalWorkspace<std::complex<float>, base_device::DEVICE_CPU>;
template class NonlocalWorkspace<std::complex<double>, base_device::DEVICE_CPU>;
#if defined(__CUDA) || defined(__ROCM)
template class NonlocalWorkspace<std::complex<float>, base_device::DEVICE_GPU>;
template class NonlocalWorkspace<std::complex<double>, base_device::DEVICE_GPU>;
#endif
} // namespace hamilt
