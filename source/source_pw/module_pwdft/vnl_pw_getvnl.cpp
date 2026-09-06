#include "vnl_pw.h"

#include "source_io/module_parameter/parameter.h"
#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "source_base/math_ylmreal.h"
#include "source_base/timer.h"
#include "source_base/module_device/device.h"
#include "source_pw/module_pwdft/kernels/vnl_op.h"

#include <vector>

/**
 * @file vnl_pw_getvnl.cpp
 * @brief Compute the Kleinman-Bylander projector functions (beta functions)
 *        in reciprocal space for a given k-point.
 *
 * This file contains the getvnl() template implementation and its explicit
 * instantiations for CPU/GPU and float/double precision.
 */

/**
 * @brief Calculate beta functions (Kleinman-Bylander projectors) with
 *        structure factor for all atoms in reciprocal space.
 *
 * Workflow:
 * 1. Collect per-atom-type counts (nbeta/nh/na) into host vectors.
 * 2. Build the |G+k| Cartesian vectors for the current k-point.
 * 3. Move data to device when GPU is enabled; reuse host vectors otherwise.
 * 4. Evaluate real spherical harmonics Ylm on all G+k vectors.
 * 5. Build structure factors sk for all atoms.
 * 6. Launch cal_vnl_op kernel to assemble vkb = beta * Ylm * sk * tab.
 *
 * @param ctx     device context (CPU or GPU)
 * @param ucell   unit cell (read-only)
 * @param ik      index of the k-point
 * @param vkb_in  output projector array of shape (nkb, npw)
 */
template <typename FPTYPE, typename Device>
void pseudopot_cell_vnl::getvnl(Device* ctx,
                                const UnitCell& ucell,
                                const int& ik,
                                std::complex<FPTYPE>* vkb_in) const
{
    ModuleBase::TITLE("pseudopot_cell_vnl", "getvnl");
    ModuleBase::timer::start("pp_cell_vnl", "getvnl");

    using cal_vnl_op = hamilt::cal_vnl_op<FPTYPE, Device>;
    using resmem_int_op = base_device::memory::resize_memory_op<int, Device>;
    using delmem_int_op = base_device::memory::delete_memory_op<int, Device>;
    using syncmem_int_op = base_device::memory::synchronize_memory_op<int, Device, base_device::DEVICE_CPU>;
    using resmem_var_op = base_device::memory::resize_memory_op<FPTYPE, Device>;
    using delmem_var_op = base_device::memory::delete_memory_op<FPTYPE, Device>;
    using castmem_var_h2d_op = base_device::memory::cast_memory_op<FPTYPE, double, Device, base_device::DEVICE_CPU>;
    using castmem_var_h2h_op
        = base_device::memory::cast_memory_op<FPTYPE, double, base_device::DEVICE_CPU, base_device::DEVICE_CPU>;
    using resmem_complex_op = base_device::memory::resize_memory_op<std::complex<FPTYPE>, Device>;
    using delmem_complex_op = base_device::memory::delete_memory_op<std::complex<FPTYPE>, Device>;

    if (lmaxkb < 0)
    {
        return;
    }

    const int x1 = (lmaxkb + 1) * (lmaxkb + 1);
    const int npw = this->wfcpw->npwk[ik];

    int* atom_nh = nullptr;
    int* atom_na = nullptr;
    int* atom_nb = nullptr;
    std::vector<int> h_atom_nh(ucell.ntype);
    std::vector<int> h_atom_na(ucell.ntype);
    std::vector<int> h_atom_nb(ucell.ntype);
    for (int it = 0; it < ucell.ntype; it++)
    {
        h_atom_nb[it] = ucell.atoms[it].ncpp.nbeta;
        h_atom_nh[it] = ucell.atoms[it].ncpp.nh;
        h_atom_na[it] = ucell.atoms[it].na;
    }
    // When the internal memory is large enough, it is better to make vkb1 be the number of pseudopot_cell_vnl.
    // We only need to initialize it once as long as the cell is unchanged.
    FPTYPE* vkb1 = nullptr;
    FPTYPE* gk = nullptr;
    FPTYPE* ylm = nullptr;
    FPTYPE* tab_ptr = this->get_tab_data<FPTYPE>();
    FPTYPE* indv_ptr = this->get_indv_data<FPTYPE>();
    FPTYPE* nhtol_ptr = this->get_nhtol_data<FPTYPE>();
    FPTYPE* nhtolm_ptr = this->get_nhtolm_data<FPTYPE>();
    resmem_var_op()(ylm, x1 * npw, "VNL::ylm");
    resmem_var_op()(vkb1, nhm * npw, "VNL::vkb1");

    std::vector<ModuleBase::Vector3<double>> gk_vec(npw);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int ig = 0; ig < npw; ig++)
    {
        gk_vec[ig] = this->wfcpw->getgpluskcar(ik, ig);
    }
    if (this->use_gpu_)
    {
        resmem_int_op()(atom_nh, ucell.ntype);
        resmem_int_op()(atom_nb, ucell.ntype);
        resmem_int_op()(atom_na, ucell.ntype);
        syncmem_int_op()(atom_nh, h_atom_nh.data(), ucell.ntype);
        syncmem_int_op()(atom_nb, h_atom_nb.data(), ucell.ntype);
        syncmem_int_op()(atom_na, h_atom_na.data(), ucell.ntype);

        resmem_var_op()(gk, npw * 3);
        castmem_var_h2d_op()(gk, reinterpret_cast<double*>(gk_vec.data()), npw * 3);
    }
    else
    {
        atom_nh = h_atom_nh.data();
        atom_nb = h_atom_nb.data();
        atom_na = h_atom_na.data();
        if (std::is_same<FPTYPE, float>::value)
        {
            resmem_var_op()(gk, npw * 3);
            castmem_var_h2h_op()(gk, reinterpret_cast<double*>(gk_vec.data()), npw * 3);
        }
        else
        {
            gk = reinterpret_cast<FPTYPE*>(gk_vec.data());
        }
    }

    ModuleBase::YlmReal::Ylm_Real(ctx, x1, npw, gk, ylm);

    std::complex<FPTYPE>* sk = nullptr;
    resmem_complex_op()(sk, ucell.nat * npw);
    this->psf->get_sk(ctx, ik, this->wfcpw, sk);

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
                 static_cast<FPTYPE>(PARAM.globalv.dq),
                 static_cast<FPTYPE>(ucell.tpiba),
                 static_cast<std::complex<FPTYPE>>(ModuleBase::NEG_IMAG_UNIT),
                 gk,
                 ylm,
                 indv_ptr,
                 nhtol_ptr,
                 nhtolm_ptr,
                 tab_ptr,
                 vkb1,
                 sk,
                 vkb_in);

    delmem_var_op()(ylm);
    delmem_var_op()(vkb1);
    delmem_complex_op()(sk);
    if (this->use_gpu_ || std::is_same<FPTYPE, float>::value)
    {
        delmem_var_op()(gk);
    }
    if (this->use_gpu_)
    {
        delmem_int_op()(atom_nh);
        delmem_int_op()(atom_nb);
        delmem_int_op()(atom_na);
    }
    ModuleBase::timer::end("pp_cell_vnl", "getvnl");
} // end subroutine getvnl

// Explicit instantiations for CPU/GPU and float/double precision.
// These must stay in the same translation unit as the template definition.
template void pseudopot_cell_vnl::getvnl<float, base_device::DEVICE_CPU>(base_device::DEVICE_CPU*,
                                                                         const UnitCell&,
                                                                         int const&,
                                                                         std::complex<float>*) const;
template void pseudopot_cell_vnl::getvnl<double, base_device::DEVICE_CPU>(base_device::DEVICE_CPU*,
                                                                          const UnitCell&,
                                                                          int const&,
                                                                          std::complex<double>*) const;
#if defined(__CUDA) || defined(__ROCM)
template void pseudopot_cell_vnl::getvnl<float, base_device::DEVICE_GPU>(base_device::DEVICE_GPU*,
                                                                         const UnitCell&,
                                                                         int const&,
                                                                         std::complex<float>*) const;
template void pseudopot_cell_vnl::getvnl<double, base_device::DEVICE_GPU>(base_device::DEVICE_GPU*,
                                                                          const UnitCell&,
                                                                          int const&,
                                                                          std::complex<double>*) const;
#endif
