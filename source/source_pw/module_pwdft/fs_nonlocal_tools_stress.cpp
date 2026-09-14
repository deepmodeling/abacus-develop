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

// cal_dbecp
template <typename FPTYPE, typename Device>
void FS_Nonlocal_tools<FPTYPE, Device>::cal_vkb_deri_s(const int& ik,
                                                       const int& nbdall,
                                                       const int& ipol,
                                                       const int& jpol)
{
    ModuleBase::TITLE("FS_Nonlocal_tools", "cal_vkb_deri_s");
    const int npol = this->ucell_->get_npol();
    const int size_becp = nbdall * npol * this->nkb;
    if (this->dbecp == nullptr)
    {
        resmem_complex_op()(dbecp, size_becp);
    }

    // prepare math tools
    Nonlocal_maths<FPTYPE, Device> maths(this->nlpp_, this->ucell_);

    const int npw = this->wfc_basis_->npwk[ik];
    std::complex<FPTYPE>* vkb_deri_ptr = this->ppcell_vkb;

    if (this->pre_ik_s != ik)
    { // k point has changed, we need to recalculate the g_plus_k
        // this->g_plus_k = maths.cal_gk(ik, this->wfc_basis_); //has been calculated by cal_becp

        const int lmax_ = this->nlpp_->lmaxkb;
        // prepare ylm，size: (lmax+1)^2 * this->max_npw
        // maths.cal_ylm(lmax_, npw, g_plus_k.data(), hd_ylm); //has been calculated by cal_becp
        maths.cal_ylm_deri(lmax_, npw, g_plus_k.data(), hd_ylm_deri);
        this->pre_ik_s = ik;
    }
    FPTYPE *gk = g_plus_k.data(), *vq_tb = this->nlpp_->tab.ptr;
    std::complex<FPTYPE>* d_sk = this->hd_sk;
    if (this->device == base_device::GpuDevice)
    {
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
        cal_vq_deri_op()(this->ctx,
                         vq_tb,
                         it,
                         gk,
                         npw,
                         this->nlpp_->tab.getBound2(),
                         this->nlpp_->tab.getBound3(),
                         PARAM.globalv.dq,
                         this->ucell_->atoms[it].ncpp.nbeta,
                         hd_vq_deri);

        // prepare（-i）^l, size: nh
        const int nh = this->ucell_->atoms[it].ncpp.nh;
        std::vector<std::complex<double>> pref = maths.cal_pref(it, nh);
        // prepare indexes for calculate vkb_deri
        this->dvkb_indexes.resize(nh * 4);
        maths.cal_dvkb_index(this->ucell_->atoms[it].ncpp.nbeta,
                             this->nlpp_->nhtol.c,
                             this->nlpp_->nhtol.nc,
                             npw,
                             it,
                             ipol,
                             jpol,
                             this->dvkb_indexes.data());
        if (this->device == base_device::GpuDevice)
        {
            syncmem_int_h2d_op()(d_dvkb_indexes, dvkb_indexes.data(), nh * 4);
            syncmem_complex_h2d_op()(d_pref_in, pref.data(), nh);
        }
        for (int ia = 0; ia < h_atom_na[it]; ia++)
        {
            // 2. calculate dbecp：
            // 2.a. calculate dbecp_noevc, repeat use the memory of ppcell.vkb

            if (this->device == base_device::CpuDevice)
            {
                d_dvkb_indexes = dvkb_indexes.data();
                d_pref_in = pref.data();
                d_g_plus_k = g_plus_k.data();
            }
            cal_vkb_deri_op()(this->ctx,
                              nh,
                              npw,
                              ipol,
                              jpol,
                              d_dvkb_indexes,
                              hd_vq,
                              hd_vq_deri,
                              hd_ylm,
                              hd_ylm_deri,
                              d_sk,
                              d_pref_in,
                              d_g_plus_k,
                              vkb_deri_ptr);
            d_sk += npw;
            vkb_deri_ptr += nh * npw;
        }
    }
}

// cal_dbecp
template <typename FPTYPE, typename Device>
void FS_Nonlocal_tools<FPTYPE, Device>::cal_dbecp_s(const int& ik,
                                                    const int& npm,
                                                    const std::complex<FPTYPE>* ppsi,
                                                    const int& nbd0)
{
    ModuleBase::TITLE("FS_Nonlocal_tools", "cal_dbecp_s");
    const int npol = this->ucell_->get_npol();
    const int npm_npol = npm * npol;
    const int npw = this->wfc_basis_->npwk[ik];
    std::complex<FPTYPE>* dbecp_ptr = this->dbecp + nbd0 * npol * this->nkb;

    // 2.b calculate dbecp = dbecp_noevc * psi
    const char transa = 'C';
    const char transb = 'N';
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
              dbecp_ptr,
              this->nkb);
}

template <typename FPTYPE, typename Device>
void FS_Nonlocal_tools<FPTYPE, Device>::cal_stress(const int& ik,
                                                   const int& npm,
                                                   const bool& occ,
                                                   const int& ipol,
                                                   const int& jpol,
                                                   FPTYPE* stress,
                                                   const int& nbd0)
{
    const int npol = this->ucell_->get_npol();
    const int index0 = nbd0 * npol * this->nkb;
    // calculate stress for target (ipol, jpol)
    if (npol == 1)
    {
        const int current_spin = this->kv_->isk[ik];
        FPTYPE* d_ekb_ik = nullptr;
        if (d_ekb != nullptr)
        {
            d_ekb_ik = d_ekb + this->nbands * ik;
        }
        FPTYPE* d_wg_ik = d_wk + ik;
        if (occ)
        {
            d_wg_ik = d_wg + this->nbands * ik;
        }
        cal_stress_nl_op()(this->ctx,
                           nondiagonal,
                           ipol,
                           jpol,
                           nkb,
                           npm,
                           this->ntype,
                           current_spin, // uspp only
                           this->nlpp_->deeq.getBound2(),
                           this->nlpp_->deeq.getBound3(),
                           this->nlpp_->deeq.getBound4(),
                           atom_nh,
                           atom_na,
                           d_wg_ik,
                           occ,
                           d_ekb_ik,
                           qq_nt,
                           deeq,
                           becp + index0,
                           dbecp + index0,
                           stress);
    }
    else
    {
        FPTYPE* d_ekb_ik = nullptr;
        if (d_ekb != nullptr)
        {
            d_ekb_ik = d_ekb + this->nbands * ik;
        }
        FPTYPE* d_wg_ik = d_wk + ik;
        if (occ)
        {
            d_wg_ik = d_wg + this->nbands * ik;
        }
        cal_stress_nl_op()(this->ctx,
                           ipol,
                           jpol,
                           nkb,
                           npm,
                           this->ntype,
                           this->nlpp_->deeq_nc.getBound2(),
                           this->nlpp_->deeq_nc.getBound3(),
                           this->nlpp_->deeq_nc.getBound4(),
                           atom_nh,
                           atom_na,
                           d_wg_ik,
                           occ,
                           d_ekb_ik,
                           qq_nt,
                           this->nlpp_->template get_deeq_nc_data<FPTYPE>(),
                           becp + index0,
                           dbecp + index0,
                           stress);
    }
}

// template instantiation
template void FS_Nonlocal_tools<double, base_device::DEVICE_CPU>::cal_vkb_deri_s(const int&,
                                                                                  const int&,
                                                                                  const int&,
                                                                                  const int&);
template void FS_Nonlocal_tools<double, base_device::DEVICE_CPU>::cal_dbecp_s(const int&,
                                                                               const int&,
                                                                               const std::complex<double>*,
                                                                               const int&);
template void FS_Nonlocal_tools<double, base_device::DEVICE_CPU>::cal_stress(const int&,
                                                                              const int&,
                                                                              const bool&,
                                                                              const int&,
                                                                              const int&,
                                                                              double*,
                                                                              const int&);
#if ((defined __CUDA) || (defined __ROCM))
template void FS_Nonlocal_tools<double, base_device::DEVICE_GPU>::cal_vkb_deri_s(const int&,
                                                                                  const int&,
                                                                                  const int&,
                                                                                  const int&);
template void FS_Nonlocal_tools<double, base_device::DEVICE_GPU>::cal_dbecp_s(const int&,
                                                                               const int&,
                                                                               const std::complex<double>*,
                                                                               const int&);
template void FS_Nonlocal_tools<double, base_device::DEVICE_GPU>::cal_stress(const int&,
                                                                              const int&,
                                                                              const bool&,
                                                                              const int&,
                                                                              const int&,
                                                                              double*,
                                                                              const int&);
#endif

} // namespace hamilt
