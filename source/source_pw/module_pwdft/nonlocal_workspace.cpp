#include "source_pw/module_pwdft/nonlocal_workspace.h"

#include "source_base/kernels/math_kernel_op.h"
#include "source_base/parallel_device.h"
#include "source_base/timer.h"
#include "source_pw/module_pwdft/kernels/nonlocal_op.h"

#include <algorithm>
#include <limits>
#include <type_traits>
#if defined(__CUDA)
#include "source_pw/module_pwdft/kernels/projector_fold.h"
#endif

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
#if defined(__CUDA)
    delete projector_fold_;
    projector_fold_ = nullptr;
#endif
}

template <typename T, typename Device>
void NonlocalWorkspace<T, Device>::project(const T* vkb, const T* psi, const int npw,
                                         const int ldv, const int ldp, const int nkb, const int bands)
{
    const T one(1, 0);
    const T zero(0, 0);
    if (npw == 0)
    {
        // Empty PW partitions still supply zero coefficients to the collective projector reduction.
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
                                       const ModulePW::PW_Basis_K& basis, const int ik, const int spin,
                                       const int bands, const int nbasis, const int npol, const int npw,
                                       const bool is_first_node, const T* vkb, const T* psi, T* hpsi)
{
    ModuleBase::timer::start("NonlocalWorkspace", "apply");
    if (is_first_node)
    {
        const size_t count = static_cast<size_t>(nbasis) * bands / npol;
        Zero<T, Device>()(hpsi, 0, count);
    }
    // Here bands counts spin-component columns, while nbasis includes all components of one orbital.
    project_contract(cell, pp, basis, ik, spin, bands, nbasis, npol, npw, vkb, psi);
    if (pp.nkb > 0 && bands > 0 && npw > 0)
    {
        bool folded = false;
#if defined(__CUDA)
        if (projector_fold_)
        {
            folded = projector_fold_->add(bands, nbasis / npol,
                                          reinterpret_cast<const std::complex<double>*>(ps_),
                                          reinterpret_cast<std::complex<double>*>(hpsi));
        }
#endif
        if (!folded)
        {
            this->back_project(vkb, hpsi, npw, pp.vkbnc, nbasis / npol, pp.nkb, bands);
        }
    }
    ModuleBase::timer::end("NonlocalWorkspace", "apply");
}


template <typename T, typename Device>
void NonlocalWorkspace<T, Device>::project_contract(const UnitCell& cell, const pseudopot_cell_vnl& pp,
                                                   const ModulePW::PW_Basis_K& basis, int ik, int spin,
                                                   int bands, int nbasis, int npol, int npw,
                                                   const T* vkb, const T* psi)
{
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
        bool folded = false;
#if defined(__CUDA)
        if (!fold_attempted_)
        {
            fold_attempted_ = true;
            if (std::is_same<T, std::complex<double>>::value
                && std::is_same<Device, base_device::DEVICE_GPU>::value
                && basis.poolnproc == 1 && cell.ntype == 1 && vkb != nullptr)
            {
                const auto& atoms = cell.atoms[0];
                if (pp.nkb == atoms.na * atoms.ncpp.nh)
                {
                    projector_fold_ = ProjectorFold::create(atoms.na, atoms.ncpp.nh,
                                                            atoms.taud.data(), atoms.tau.data(), basis).release();
                }
            }
        }
        if (projector_fold_ && fold_ik_ != ik)
        {
            projector_fold_->prepare(ik, reinterpret_cast<const std::complex<double>*>(vkb), pp.vkbnc);
            fold_ik_ = ik;
        }
        if (projector_fold_)
        {
            folded = projector_fold_->project(bands, nbasis / npol,
                                              reinterpret_cast<const std::complex<double>*>(psi),
                                              reinterpret_cast<std::complex<double>*>(becp_));
        }
#endif
        if (!folded)
        {
            this->project(vkb, psi, npw, pp.vkbnc, nbasis / npol, pp.nkb, bands);
        }
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
        this->contract(cell, pp, spin, npol, bands);
    }
}
template class NonlocalWorkspace<std::complex<float>, base_device::DEVICE_CPU>;
template class NonlocalWorkspace<std::complex<double>, base_device::DEVICE_CPU>;
#if defined(__CUDA) || defined(__ROCM)
template class NonlocalWorkspace<std::complex<float>, base_device::DEVICE_GPU>;
template class NonlocalWorkspace<std::complex<double>, base_device::DEVICE_GPU>;
#endif
} // namespace hamilt
