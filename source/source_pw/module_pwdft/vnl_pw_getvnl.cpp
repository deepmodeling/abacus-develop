#include "source_base/global_function.h"
#include "source_base/math_ylmreal.h"
#include "source_base/module_device/device.h"
#include "source_base/timer.h"
#include "source_pw/module_pwdft/kernels/vnl_op.h"
#include "vnl_pw.h"

#include <cmath>
#include <limits>
#include <vector>

/**
 * @brief Calculate velocity-gauge nonlocal projectors at k+G+A(t).
 *
 * The radial and angular projector factors use the kinetic momentum shifted
 * by the Hartree-unit vector potential. The structure factor remains evaluated
 * on the static reciprocal grid because its atom-dependent A phase cancels in
 * each same-atom projector outer product.
 */
template <typename FPTYPE, typename Device>
void pseudopot_cell_vnl::getvnl(Device* ctx,
                                const UnitCell& ucell,
                                const int& ik,
                                const ModuleBase::Vector3<double>& vector_potential,
                                std::complex<FPTYPE>* vkb_in) const
{
    ModuleBase::timer::start("pp_cell_vnl", "getvnl");

    using cal_vnl_op = hamilt::cal_vnl_op<FPTYPE, Device>;
    using resmem_int_op = base_device::memory::resize_memory_op<int, Device>;
    using delmem_int_op = base_device::memory::delete_memory_op<int, Device>;
    using syncmem_int_op = base_device::memory::synchronize_memory_op<int, Device, base_device::DEVICE_CPU>;
    using resmem_var_op = base_device::memory::resize_memory_op<FPTYPE, Device>;
    using delmem_var_op = base_device::memory::delete_memory_op<FPTYPE, Device>;
    using castmem_var_h2d_op = base_device::memory::cast_memory_op<FPTYPE, double, Device, base_device::DEVICE_CPU>;
    using castmem_var_h2h_op = base_device::memory::cast_memory_op<FPTYPE, double, base_device::DEVICE_CPU, base_device::DEVICE_CPU>;
    using resmem_complex_op = base_device::memory::resize_memory_op<std::complex<FPTYPE>, Device>;
    using delmem_complex_op = base_device::memory::delete_memory_op<std::complex<FPTYPE>, Device>;

    if (this->lmaxkb < 0 || this->wfcpw->npwk[ik] == 0)
    {
        ModuleBase::timer::end("pp_cell_vnl", "getvnl");
        return;
    }

    const int ylm_count = (this->lmaxkb + 1) * (this->lmaxkb + 1);
    const int npw = this->wfcpw->npwk[ik];

    int* atom_nh = nullptr;
    int* atom_na = nullptr;
    int* atom_nb = nullptr;
    std::vector<int> host_atom_nh(ucell.ntype);
    std::vector<int> host_atom_na(ucell.ntype);
    std::vector<int> host_atom_nb(ucell.ntype);
    for (int it = 0; it < ucell.ntype; ++it)
    {
        host_atom_nb[it] = ucell.atoms[it].ncpp.nbeta;
        host_atom_nh[it] = ucell.atoms[it].ncpp.nh;
        host_atom_na[it] = ucell.atoms[it].na;
    }

    FPTYPE* vkb_radial = nullptr;
    FPTYPE* shifted_gk_data = nullptr;
    FPTYPE* ylm = nullptr;
    FPTYPE* tab_ptr = this->get_tab_data<FPTYPE>();
    FPTYPE* indv_ptr = this->get_indv_data<FPTYPE>();
    FPTYPE* nhtol_ptr = this->get_nhtol_data<FPTYPE>();
    FPTYPE* nhtolm_ptr = this->get_nhtolm_data<FPTYPE>();
    resmem_var_op()(ylm, ylm_count * npw, "VNL::ylm");
    resmem_var_op()(vkb_radial, this->nhm * npw, "VNL::vkb_radial");

    const ModuleBase::Vector3<double> reduced_vector_potential = vector_potential / ucell.tpiba;
    std::vector<ModuleBase::Vector3<double>> shifted_gk(npw);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int ig = 0; ig < npw; ++ig)
    {
        shifted_gk[ig] = this->wfcpw->getgpluskcar(ik, ig) + reduced_vector_potential;
    }

    // Validate the index in the kernel's precision before obtaining any radial values.
    for (int ig = 0; ig < npw; ++ig)
    {
        const FPTYPE x = static_cast<FPTYPE>(shifted_gk[ig].x);
        const FPTYPE y = static_cast<FPTYPE>(shifted_gk[ig].y);
        const FPTYPE z = static_cast<FPTYPE>(shifted_gk[ig].z);
        const FPTYPE q = std::sqrt(x * x + y * y + z * z) * static_cast<FPTYPE>(ucell.tpiba);
        const FPTYPE index = q / static_cast<FPTYPE>(this->table_dq_);
        // Also cover small host/device differences in the norm and division.
        const double guarded_index = static_cast<double>(index) * (1.0 + 8.0 * std::numeric_limits<FPTYPE>::epsilon());
        this->check_vnl_index(guarded_index, false);
    }

    if (this->use_gpu_)
    {
        resmem_int_op()(atom_nh, ucell.ntype);
        resmem_int_op()(atom_nb, ucell.ntype);
        resmem_int_op()(atom_na, ucell.ntype);
        syncmem_int_op()(atom_nh, host_atom_nh.data(), ucell.ntype);
        syncmem_int_op()(atom_nb, host_atom_nb.data(), ucell.ntype);
        syncmem_int_op()(atom_na, host_atom_na.data(), ucell.ntype);

        resmem_var_op()(shifted_gk_data, npw * 3);
        castmem_var_h2d_op()(shifted_gk_data, reinterpret_cast<double*>(shifted_gk.data()), npw * 3);
    }
    else
    {
        atom_nh = host_atom_nh.data();
        atom_nb = host_atom_nb.data();
        atom_na = host_atom_na.data();
        if (std::is_same<FPTYPE, float>::value)
        {
            resmem_var_op()(shifted_gk_data, npw * 3);
            castmem_var_h2h_op()(shifted_gk_data, reinterpret_cast<double*>(shifted_gk.data()), npw * 3);
        }
        else
        {
            shifted_gk_data = reinterpret_cast<FPTYPE*>(shifted_gk.data());
        }
    }

    ModuleBase::YlmReal::Ylm_Real(ctx, ylm_count, npw, shifted_gk_data, ylm);

    std::complex<FPTYPE>* structure_factor = nullptr;
    resmem_complex_op()(structure_factor, ucell.nat * npw);
    this->psf->get_sk(ctx, ik, this->wfcpw, structure_factor);

    cal_vnl_op()(ctx,
                 ucell.ntype,
                 npw,
                 this->wfcpw->npwk_max,
                 this->nhm,
                 this->tab.getBound2(),
                 this->tab.getBound3(),
                 atom_na,
                 atom_nb,
                 atom_nh,
                 static_cast<FPTYPE>(this->table_dq_),
                 static_cast<FPTYPE>(ucell.tpiba),
                 static_cast<std::complex<FPTYPE>>(ModuleBase::NEG_IMAG_UNIT),
                 shifted_gk_data,
                 ylm,
                 indv_ptr,
                 nhtol_ptr,
                 nhtolm_ptr,
                 tab_ptr,
                 vkb_radial,
                 structure_factor,
                 vkb_in);

    delmem_var_op()(ylm);
    delmem_var_op()(vkb_radial);
    delmem_complex_op()(structure_factor);
    if (this->use_gpu_ || std::is_same<FPTYPE, float>::value)
    {
        delmem_var_op()(shifted_gk_data);
    }
    if (this->use_gpu_)
    {
        delmem_int_op()(atom_nh);
        delmem_int_op()(atom_nb);
        delmem_int_op()(atom_na);
    }

    ModuleBase::timer::end("pp_cell_vnl", "getvnl");
}

// Explicit instantiations for CPU/GPU and float/double precision.
template void pseudopot_cell_vnl::getvnl<float, base_device::DEVICE_CPU>(base_device::DEVICE_CPU*,
                                                                         const UnitCell&,
                                                                         const int&,
                                                                         const ModuleBase::Vector3<double>&,
                                                                         std::complex<float>*) const;
template void pseudopot_cell_vnl::getvnl<double, base_device::DEVICE_CPU>(base_device::DEVICE_CPU*,
                                                                          const UnitCell&,
                                                                          const int&,
                                                                          const ModuleBase::Vector3<double>&,
                                                                          std::complex<double>*) const;
#if defined(__CUDA) || defined(__ROCM)
template void pseudopot_cell_vnl::getvnl<float, base_device::DEVICE_GPU>(base_device::DEVICE_GPU*,
                                                                         const UnitCell&,
                                                                         const int&,
                                                                         const ModuleBase::Vector3<double>&,
                                                                         std::complex<float>*) const;
template void pseudopot_cell_vnl::getvnl<double, base_device::DEVICE_GPU>(base_device::DEVICE_GPU*,
                                                                          const UnitCell&,
                                                                          const int&,
                                                                          const ModuleBase::Vector3<double>&,
                                                                          std::complex<double>*) const;
#endif
