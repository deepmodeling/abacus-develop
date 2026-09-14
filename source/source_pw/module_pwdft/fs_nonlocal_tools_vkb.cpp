#include "fs_nonlocal_tools.h"

#include "source_base/math_polyint.h"
#include "source_base/math_ylmreal.h"
#include "source_base/parallel_device.h"
#include "source_base/timer.h"
#include "source_base/tool_title.h"
#include "source_pw/module_pwdft/kernels/force_op.h"
#include "source_io/module_parameter/parameter.h"
#include "nonlocal_maths.hpp"

#include "source_base/parallel_comm.h" // different MPI worlds (POOL_WORLD)

namespace hamilt
{

// cal_vkb
template <typename FPTYPE, typename Device>
void FS_Nonlocal_tools<FPTYPE, Device>::cal_vkb(const int& ik, const int& nbdall)
{
    ModuleBase::TITLE("FS_Nonlocal_tools", "cal_vkb");
    const int npol = this->ucell_->get_npol();
    const int size_becp = nbdall * npol * this->nkb;
    if (this->becp == nullptr)
    {
        resmem_complex_op()(becp, size_becp);
    }

    // prepare math tools
    Nonlocal_maths<FPTYPE, Device> maths(this->nlpp_, this->ucell_);
    const int npw = this->wfc_basis_->npwk[ik];

    std::complex<FPTYPE>* vkb_ptr = this->ppcell_vkb;

    // calculate G+K
    this->g_plus_k = maths.cal_gk(ik, this->wfc_basis_);
    FPTYPE *gk = g_plus_k.data(), *vq_tb = this->nlpp_->tab.ptr;
    // calculate sk
    resmem_complex_op()(hd_sk, this->ucell_->nat * npw);
    this->sf_->get_sk(ctx, ik, this->wfc_basis_, hd_sk);
    std::complex<FPTYPE>* d_sk = this->hd_sk;
    // prepare ylm，size: (lmax+1)^2 * this->max_npw
    const int lmax_ = this->nlpp_->lmaxkb;
    maths.cal_ylm(lmax_, npw, g_plus_k.data(), hd_ylm);
    if (this->device == base_device::GpuDevice)
    {
        syncmem_var_h2d_op()(d_g_plus_k, g_plus_k.data(), g_plus_k.size());
        syncmem_var_h2d_op()(d_vq_tab, this->nlpp_->tab.ptr, this->nlpp_->tab.getSize());
        gk = d_g_plus_k;
        vq_tb = d_vq_tab;
    }
    for (int it = 0; it < this->ucell_->ntype; it++) // loop all elements
    {
        int lenth_vq = this->ucell_->atoms[it].ncpp.nbeta * npw;
        // prepare inputs for calculating vkb，vkb1，vkb2
        // prepare vq and vq', size: nq * this->max_npw
        std::vector<double> vq(lenth_vq); // cal_vq(it, g_plus_k.data(), npw);
        // std::vector<double> vq2(vq.size());

        cal_vq_op()(this->ctx,
                    vq_tb,
                    it,
                    gk,
                    npw,
                    this->nlpp_->tab.getBound2(),
                    this->nlpp_->tab.getBound3(),
                    PARAM.globalv.dq,
                    this->ucell_->atoms[it].ncpp.nbeta,
                    hd_vq);

        // prepare（-i）^l, size: nh
        const int nh = this->ucell_->atoms[it].ncpp.nh;
        std::vector<std::complex<double>> pref = maths.cal_pref(it, nh);
        this->dvkb_indexes.resize(nh * 4);
        maths.cal_dvkb_index(this->ucell_->atoms[it].ncpp.nbeta,
                             this->nlpp_->nhtol.c,
                             this->nlpp_->nhtol.nc,
                             npw,
                             it,
                             0,
                             0,
                             this->dvkb_indexes.data());
        if (this->device == base_device::GpuDevice)
        {
            syncmem_int_h2d_op()(d_dvkb_indexes, dvkb_indexes.data(), nh * 4);
            syncmem_complex_h2d_op()(d_pref_in, pref.data(), nh);
        }

        for (int ia = 0; ia < h_atom_na[it]; ia++)
        {
            // 1. calculate becp
            // 1.a calculate vkb

            if (this->device == base_device::CpuDevice)
            {
                d_pref_in = pref.data();
                d_dvkb_indexes = dvkb_indexes.data();
            }
            cal_vkb_op()(this->ctx, nh, npw, d_dvkb_indexes, hd_vq, hd_ylm, d_sk, d_pref_in, vkb_ptr);

            // 2.b calculate becp = vkb * psi
            vkb_ptr += nh * npw;
            d_sk += npw;
        }
    }
}

// cal_becp
template <typename FPTYPE, typename Device>
void FS_Nonlocal_tools<FPTYPE, Device>::cal_becp(const int& ik,
                                                 const int& npm,
                                                 const std::complex<FPTYPE>* ppsi,
                                                 const int& nbd0)
{
    if (npm == 0)
    {
        return;
    }
    const int npol = this->ucell_->get_npol();
    const int npw = this->wfc_basis_->npwk[ik];
    const char transa = 'C';
    const char transb = 'N';
    const int npm_npol = npm * npol;
    const int index0 = nbd0 * npol * nkb;
    gemm_op()(transa,
              transb,
              this->nkb,
              npm_npol,
              npw,
              &ModuleBase::ONE,
              this->ppcell_vkb,
              npw,
              ppsi,
              this->max_npw,
              &ModuleBase::ZERO,
              this->becp + index0,
              this->nkb);
}

template <typename FPTYPE, typename Device>
void FS_Nonlocal_tools<FPTYPE, Device>::reduce_pool_becp(const int& npm)
{
    const int npol = this->ucell_->get_npol();
    const int size_becp_act = npm * npol * this->nkb;
    // becp calculate is over , now we should broadcast this data.
#ifdef __MPI
    if (GlobalV::NPROC_IN_POOL > 1)
    {
        Parallel_Common::reduce_data(this->becp, size_becp_act, POOL_WORLD);
    }
#endif
}

// template instantiation
template void FS_Nonlocal_tools<double, base_device::DEVICE_CPU>::cal_vkb(const int&, const int&);
template void FS_Nonlocal_tools<double, base_device::DEVICE_CPU>::cal_becp(const int&,
                                                                            const int&,
                                                                            const std::complex<double>*,
                                                                            const int&);
template void FS_Nonlocal_tools<double, base_device::DEVICE_CPU>::reduce_pool_becp(const int&);
#if ((defined __CUDA) || (defined __ROCM))
template void FS_Nonlocal_tools<double, base_device::DEVICE_GPU>::cal_vkb(const int&, const int&);
template void FS_Nonlocal_tools<double, base_device::DEVICE_GPU>::cal_becp(const int&,
                                                                            const int&,
                                                                            const std::complex<double>*,
                                                                            const int&);
template void FS_Nonlocal_tools<double, base_device::DEVICE_GPU>::reduce_pool_becp(const int&);
#endif

} // namespace hamilt
