#include "source_pw/module_pwdft/nonlocal_workspace.h"

#include "source_base/kernels/math_kernel_op.h"
#include "source_base/parallel_device.h"
#include "source_base/timer.h"
#include "source_pw/module_pwdft/kernels/nonlocal_op.h"

#include <algorithm>
#include <limits>

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
        const size_t count = static_cast<size_t>(nkb) * bands;
        Zero<T, Device>()(becp_, 0, count);
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
    const size_t count = static_cast<size_t>(pp.nkb) * bands;
    Zero<T, Device>()(ps_, 0, count);
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
        const size_t count = static_cast<size_t>(nbasis) * bands / npol;
        Zero<T, Device>()(hpsi, 0, count);
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
            // The complex reduction wrapper doubles the MPI element count.
            const size_t max_chunk = std::numeric_limits<int>::max() / 2;
            for (size_t offset = 0; offset < count;)
            {
                const size_t remaining = count - offset;
                const int chunk_size = static_cast<int>(std::min(max_chunk, remaining));
                T* chunk = becp_ + offset;
                Parallel_Common::reduce_dev<T, Device>(chunk, chunk_size, basis.pool_world);
                offset += chunk_size;
            }
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
