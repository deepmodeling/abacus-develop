#include "fs_nonlocal_tools.h"

#include "module_base/math_polyint.h"
#include "module_base/math_ylmreal.h"
#include "module_base/memory.h"
#include "module_base/parallel_device.h"
#include "module_base/timer.h"
#include "module_base/tool_title.h"
#include "module_hamilt_pw/hamilt_pwdft/kernels/force_op.h"
#include "module_parameter/parameter.h"
#include "nonlocal_maths.hpp"

namespace hamilt
{

template <typename FPTYPE, typename Device>
FS_Nonlocal_tools<FPTYPE, Device>::FS_Nonlocal_tools(const pseudopot_cell_vnl* nlpp_in,
                                                     const UnitCell* ucell_in,
                                                     const K_Vectors* kv_in,
                                                     const ModulePW::PW_Basis_K* wfc_basis_in,
                                                     const Structure_Factor* sf_in,
                                                     const ModuleBase::matrix& wg,
                                                     const ModuleBase::matrix* p_ekb)
    : nlpp_(nlpp_in), ucell_(ucell_in), kv_(kv_in), wfc_basis_(wfc_basis_in), sf_(sf_in)
{
    // get the device context
    this->device = base_device::get_device_type<Device>(this->ctx);
    this->nkb = nlpp_->nkb;
    this->max_npw = wfc_basis_->npwk_max;
    this->ntype = ucell_->ntype;
    this->nbands = wg.nc;

    // There is a contribution for jh<>ih in US case or multi projectors case
    // Actually, the judge of nondiagonal should be done on every atom type
    this->nondiagonal = (PARAM.globalv.use_uspp || this->nlpp_->multi_proj) ? true : false;

    // allocate memory
    this->allocate_memory(wg, p_ekb);
}

template <typename FPTYPE, typename Device>
FS_Nonlocal_tools<FPTYPE, Device>::~FS_Nonlocal_tools()
{
    // delete memory
    delete_memory();
}

template <typename FPTYPE, typename Device>
void FS_Nonlocal_tools<FPTYPE, Device>::allocate_memory(const ModuleBase::matrix& wg, const ModuleBase::matrix* p_ekb)
{
    // allocate memory

    // prepare the memory of stress and init some variables:
    this->h_atom_nh.resize(this->ntype);
    this->h_atom_na.resize(this->ntype);
    for (int ii = 0; ii < this->ntype; ii++)
    {
        h_atom_nh[ii] = this->ucell_->atoms[ii].ncpp.nh;
        h_atom_na[ii] = this->ucell_->atoms[ii].na;
    }

    this->deeq = this->nlpp_->template get_deeq_data<FPTYPE>();
    this->kvec_c = this->wfc_basis_->template get_kvec_c_data<FPTYPE>();
    this->qq_nt = this->nlpp_->template get_qq_nt_data<FPTYPE>();

    int max_nbeta = 0;
    for (int it = 0; it < this->ntype; it++) // loop all elements
    {
        max_nbeta = std::max(this->ucell_->atoms[it].ncpp.nbeta, max_nbeta);
        this->max_nh = std::max(this->ucell_->atoms[it].ncpp.nh, max_nh);
    }

    // allocate the memory for vkb and vkb_deri.
    if (this->device == base_device::GpuDevice)
    {
        resmem_int_op()(this->ctx, this->d_dvkb_indexes, max_nh * 4);
    }

    resmem_var_op()(this->ctx, this->hd_vq, max_nbeta * max_npw);
    resmem_var_op()(this->ctx, this->hd_vq_deri, max_nbeta * max_npw);
    const int _lmax = this->nlpp_->lmaxkb;
    resmem_var_op()(this->ctx, this->hd_ylm, (_lmax + 1) * (_lmax + 1) * max_npw);
    resmem_var_op()(this->ctx, this->hd_ylm_deri, 3 * (_lmax + 1) * (_lmax + 1) * max_npw);
    const int nks = this->kv_->get_nks();
    resmem_var_op()(this->ctx, d_wk, nks);
    syncmem_var_h2d_op()(this->ctx, this->cpu_ctx, d_wk, this->kv_->wk.data(), nks);

    if (this->device == base_device::GpuDevice)
    {
        resmem_var_op()(this->ctx, d_wg, wg.nr * wg.nc);
        syncmem_var_h2d_op()(this->ctx, this->cpu_ctx, d_wg, wg.c, wg.nr * wg.nc);
        if (p_ekb != nullptr)
        {
            resmem_var_op()(this->ctx, d_ekb, p_ekb->nr * p_ekb->nc);
            syncmem_var_h2d_op()(this->ctx, this->cpu_ctx, d_ekb, p_ekb->c, p_ekb->nr * p_ekb->nc);
        }
        resmem_int_op()(this->ctx, atom_nh, this->ntype);
        resmem_int_op()(this->ctx, atom_na, this->ntype);
        syncmem_int_h2d_op()(this->ctx, this->cpu_ctx, atom_nh, h_atom_nh.data(), this->ntype);
        syncmem_int_h2d_op()(this->ctx, this->cpu_ctx, atom_na, h_atom_na.data(), this->ntype);

        resmem_var_op()(this->ctx, d_g_plus_k, max_npw * 5);
        resmem_var_op()(this->ctx, d_pref, max_nh);
        resmem_var_op()(this->ctx, d_vq_tab, this->nlpp_->tab.getSize());
        resmem_complex_op()(this->ctx, d_pref_in, max_nh);

        this->ppcell_vkb = this->nlpp_->template get_vkb_data<FPTYPE>();
    }
    else
    {
        this->d_wg = wg.c;
        if (p_ekb != nullptr)
        {
            this->d_ekb = p_ekb->c;
        }
        this->atom_nh = h_atom_nh.data();
        this->atom_na = h_atom_na.data();
        this->ppcell_vkb = this->nlpp_->vkb.c;
    }
}

template <typename FPTYPE, typename Device>
void FS_Nonlocal_tools<FPTYPE, Device>::delete_memory()
{
    // delete memory

    delmem_var_op()(this->ctx, hd_vq);
    delmem_var_op()(this->ctx, hd_vq_deri);
    delmem_var_op()(this->ctx, hd_ylm);
    delmem_var_op()(this->ctx, hd_ylm_deri);
    delmem_var_op()(this->ctx, d_wk);

    // delete memory on GPU
    if (this->device == base_device::GpuDevice)
    {
        delmem_var_op()(this->ctx, d_wg);
        delmem_var_op()(this->ctx, d_ekb);
        delmem_int_op()(this->ctx, atom_nh);
        delmem_int_op()(this->ctx, atom_na);
        delmem_var_op()(this->ctx, d_g_plus_k);
        delmem_var_op()(this->ctx, d_pref);
        delmem_var_op()(this->ctx, d_vq_tab);
        delmem_complex_op()(this->ctx, this->d_pref_in);
        delmem_int_op()(this->ctx, d_dvkb_indexes);
    }

    if (becp != nullptr)
    {
        delmem_complex_op()(this->ctx, becp);
        delmem_complex_op()(this->ctx, hd_sk);
    }
    if (dbecp != nullptr)
    {
        delmem_complex_op()(this->ctx, dbecp);
    }
    if (this->pre_ik_f != -1)
    {
        delmem_int_op()(this->ctx, gcar_zero_indexes);
        delmem_complex_op()(this->ctx, vkb_save);
        delmem_var_op()(this->ctx, gcar);
    }
}

// cal_vkb
template <typename FPTYPE, typename Device>
void FS_Nonlocal_tools<FPTYPE, Device>::cal_vkb(const int& ik, const int& nbdall)
{
    ModuleBase::TITLE("FS_Nonlocal_tools", "cal_vkb");
    const int npol = this->ucell_->get_npol();
    const int size_becp = nbdall * npol * this->nkb;
    if (this->becp == nullptr)
    {
        resmem_complex_op()(this->ctx, becp, size_becp);
    }

    // prepare math tools
    Nonlocal_maths<FPTYPE, Device> maths(this->nlpp_, this->ucell_);
    const int npw = this->wfc_basis_->npwk[ik];

    std::complex<FPTYPE>* vkb_ptr = this->ppcell_vkb;

    // calculate G+K
    this->g_plus_k = maths.cal_gk(ik, this->wfc_basis_);
    FPTYPE *gk = g_plus_k.data(), *vq_tb = this->nlpp_->tab.ptr;
    // calculate sk
    resmem_complex_op()(ctx, hd_sk, this->ucell_->nat * npw);
    this->sf_->get_sk(ctx, ik, this->wfc_basis_, hd_sk);
    std::complex<FPTYPE>* d_sk = this->hd_sk;
    // prepare ylm，size: (lmax+1)^2 * this->max_npw
    const int lmax_ = this->nlpp_->lmaxkb;
    maths.cal_ylm(lmax_, npw, g_plus_k.data(), hd_ylm);
    if (this->device == base_device::GpuDevice)
    {
        syncmem_var_h2d_op()(this->ctx, this->cpu_ctx, d_g_plus_k, g_plus_k.data(), g_plus_k.size());
        syncmem_var_h2d_op()(this->ctx, this->cpu_ctx, d_vq_tab, this->nlpp_->tab.ptr, this->nlpp_->tab.getSize());
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
        std::vector<std::complex<double>> pref = maths.cal_pref(it);
        const int nh = pref.size();
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
            syncmem_int_h2d_op()(this->ctx, this->cpu_ctx, d_dvkb_indexes, dvkb_indexes.data(), nh * 4);
            syncmem_complex_h2d_op()(this->ctx, this->cpu_ctx, d_pref_in, pref.data(), nh);
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
    gemm_op()(this->ctx,
              transa,
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
    ModuleBase::timer::tick("FS_Nonlocal_tools", "cal_becp");
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
        Parallel_Common::reduce_dev(this->ctx, this->becp, size_becp_act, POOL_WORLD);
    }
#endif
}

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
        resmem_complex_op()(this->ctx, dbecp, size_becp);
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
        std::vector<std::complex<double>> pref = maths.cal_pref(it);
        int nh = pref.size();
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
            syncmem_int_h2d_op()(this->ctx, this->cpu_ctx, d_dvkb_indexes, dvkb_indexes.data(), nh * 4);
            syncmem_complex_h2d_op()(this->ctx, this->cpu_ctx, d_pref_in, pref.data(), nh);
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
    gemm_op()(this->ctx,
              transa,
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

template <typename FPTYPE, typename Device>
void FS_Nonlocal_tools<FPTYPE, Device>::cal_vkb_deri_f(const int& ik, const int& nbdall, const int& ipol)
{
    ModuleBase::TITLE("FS_Nonlocal_tools", "cal_vkb_deri");
    const int npol = this->ucell_->get_npol();
    const int size_becp = nbdall * npol * this->nkb;
    if (this->dbecp == nullptr)
    {
        resmem_complex_op()(this->ctx, dbecp, 3 * size_becp);
    }

    const std::complex<FPTYPE>* vkb_ptr = this->ppcell_vkb;
    std::complex<FPTYPE>* vkb_deri_ptr = this->ppcell_vkb;

    const int npw = this->wfc_basis_->npwk[ik];
    if (this->pre_ik_f == -1)
    {
        resmem_var_op()(this->ctx, gcar, 3 * this->wfc_basis_->npwk_max);
        resmem_int_op()(this->ctx, gcar_zero_indexes, 3 * this->wfc_basis_->npwk_max);
    }

    if (this->pre_ik_f != ik)
    {
        this->transfer_gcar(npw,
                            this->wfc_basis_->npwk_max,
                            &(this->wfc_basis_->gcar[ik * this->wfc_basis_->npwk_max].x));
    }

    this->save_vkb(ik, ipol);

    const std::complex<double> coeff = ipol == 0 ? ModuleBase::NEG_IMAG_UNIT : ModuleBase::ONE;

    // calculate the vkb_deri for ipol with the memory of ppcell_vkb
    cal_vkb1_nl_op<FPTYPE, Device>()(this->ctx, nkb, npw, npw, npw, ipol, coeff, vkb_ptr, gcar, vkb_deri_ptr);

}

template <typename FPTYPE, typename Device>
void FS_Nonlocal_tools<FPTYPE, Device>::cal_dbecp_f(const int& ik,
                                                    const int& nbdall,
                                                    const int& npm,
                                                    const int& ipol,
                                                    const std::complex<FPTYPE>* ppsi,
                                                    const int& nbd0)
{
    ModuleBase::TITLE("FS_Nonlocal_tools", "cal_dbecp_f");
    const int npol = this->ucell_->get_npol();
    const int size_becp = nbdall * npol * this->nkb;
    const int index0 = nbd0 * npol * this->nkb;
    std::complex<FPTYPE>* dbecp_ptr = this->dbecp + ipol * size_becp + index0;
    std::complex<FPTYPE>* vkb_deri_ptr = this->ppcell_vkb;
    const int npm_npol = npm * npol;
    const int npw = this->wfc_basis_->npwk[ik];

    // do gemm to get dbecp and revert the ppcell_vkb for next ipol
    const char transa = 'C';
    const char transb = 'N';
    gemm_op()(this->ctx,
              transa,
              transb,
              this->nkb,
              npm_npol,
              npw,
              &ModuleBase::ONE,
              vkb_deri_ptr,
              npw,
              ppsi,
              this->max_npw,
              &ModuleBase::ZERO,
              dbecp_ptr,
              this->nkb);
}

// save_vkb
template <typename FPTYPE, typename Device>
void FS_Nonlocal_tools<FPTYPE, Device>::save_vkb(const int& ik, const int& ipol)
{
    const int npw = this->wfc_basis_->npwk[ik];
    if (this->device == base_device::CpuDevice)
    {
        const int gcar_zero_count = this->gcar_zero_indexes[ipol * this->wfc_basis_->npwk_max];
        const int* gcar_zero_ptrs = &this->gcar_zero_indexes[ipol * this->wfc_basis_->npwk_max + 1];
        const std::complex<FPTYPE>* vkb_ptr = this->ppcell_vkb;
        std::complex<FPTYPE>* vkb_save_ptr = this->vkb_save;
        // find the zero indexes to save the vkb values to vkb_save
        for (int ikb = 0; ikb < this->nkb; ++ikb)
        {
            for (int icount = 0; icount < gcar_zero_count; ++icount)
            {
                *vkb_save_ptr = vkb_ptr[gcar_zero_ptrs[icount]];
                ++vkb_save_ptr;
            }
            vkb_ptr += npw;
        }
    }
    else
    {
#if __CUDA || __UT_USE_CUDA || __ROCM || __UT_USE_ROCM
        saveVkbValues<FPTYPE>(this->gcar_zero_indexes,
                              this->ppcell_vkb,
                              this->vkb_save,
                              nkb,
                              this->gcar_zero_counts[ipol],
                              npw,
                              ipol,
                              this->wfc_basis_->npwk_max);
#endif
    }
}

// revert_vkb
template <typename FPTYPE, typename Device>
void FS_Nonlocal_tools<FPTYPE, Device>::revert_vkb(const int& ik, const int& ipol)
{
    const int npw = this->wfc_basis_->npwk[ik];
    const std::complex<FPTYPE> coeff = ipol == 0 ? ModuleBase::NEG_IMAG_UNIT : ModuleBase::ONE;
    if (this->device == base_device::CpuDevice)
    {
        const int gcar_zero_count = this->gcar_zero_indexes[ipol * this->wfc_basis_->npwk_max];
        const int* gcar_zero_ptrs = &this->gcar_zero_indexes[ipol * this->wfc_basis_->npwk_max + 1];
        std::complex<FPTYPE>* vkb_ptr = this->ppcell_vkb;
        const std::complex<FPTYPE>* vkb_save_ptr = this->vkb_save;
        // find the zero indexes to save the vkb values to vkb_save
        for (int ikb = 0; ikb < this->nkb; ++ikb)
        {
            for (int icount = 0; icount < gcar_zero_count; ++icount)
            {
                vkb_ptr[gcar_zero_ptrs[icount]] = *vkb_save_ptr * coeff;
                ++vkb_save_ptr;
            }
            vkb_ptr += npw;
        }
    }
    else
    {
#if __CUDA || __UT_USE_CUDA || __ROCM || __UT_USE_ROCM
        revertVkbValues<FPTYPE>(this->gcar_zero_indexes,
                                this->ppcell_vkb,
                                this->vkb_save,
                                nkb,
                                this->gcar_zero_counts[ipol],
                                npw,
                                ipol,
                                this->wfc_basis_->npwk_max,
                                coeff);
#endif
    }
    this->pre_ik_f = ik;
}

template <typename FPTYPE, typename Device>
void FS_Nonlocal_tools<FPTYPE, Device>::transfer_gcar(const int& npw, const int& npw_max, const FPTYPE* gcar_in)
{
    std::vector<FPTYPE> gcar_tmp(3 * npw_max);
    gcar_tmp.assign(gcar_in, gcar_in + 3 * npw_max);
    std::vector<int> gcar_zero_indexes_tmp(3 * npw_max);

    int* gcar_zero_ptrs[3];
    for (int i = 0; i < 3; i++)
    {
        gcar_zero_ptrs[i] = &gcar_zero_indexes_tmp[i * npw_max];
        gcar_zero_ptrs[i][0] = -1;
        this->gcar_zero_counts[i] = 0;
    }
    for (int ig = 0; ig < npw; ig++)
    {
        // calculate gcar.x , gcar.y/gcar.x, gcar.z/gcar.y
        // if individual gcar is less than 1e-15, we will record the index
        for (int i = 0; i < 3; ++i)
        {
            if (std::abs(gcar_tmp[ig * 3 + i]) < 1e-15)
            {
                ++gcar_zero_counts[i];
                gcar_zero_ptrs[i][gcar_zero_counts[i]] = ig;
            }
        }
        // four cases for the gcar of y and z
        if (gcar_zero_ptrs[0][gcar_zero_counts[0]] == ig && gcar_zero_ptrs[1][gcar_zero_counts[1]] == ig)
        { // x == y == 0, z = z
        }
        else if (gcar_zero_ptrs[0][gcar_zero_counts[0]] != ig && gcar_zero_ptrs[1][gcar_zero_counts[1]] == ig)
        { // x != 0, y == 0, z = z/x
            gcar_tmp[ig * 3 + 2] /= gcar_tmp[ig * 3];
        }
        else if (gcar_zero_ptrs[0][gcar_zero_counts[0]] == ig && gcar_zero_ptrs[1][gcar_zero_counts[1]] != ig)
        { // x == 0, y != 0, y = y, z = z/y
            gcar_tmp[ig * 3 + 2] /= gcar_tmp[ig * 3 + 1];
        }
        else
        { // x != 0, y != 0, y = y/x, z = z/y
            gcar_tmp[ig * 3 + 2] /= gcar_tmp[ig * 3 + 1];
            gcar_tmp[ig * 3 + 1] /= gcar_tmp[ig * 3];
        }
    }
    for (int i = 0; i < 3; ++i)
    { // record the counts to the first element
        gcar_zero_ptrs[i][0] = gcar_zero_counts[i];
    }
    // prepare the memory for vkb_save
    const int max_count = std::max(gcar_zero_counts[0], std::max(gcar_zero_counts[1], gcar_zero_counts[2]));
    resmem_complex_op()(this->ctx, this->vkb_save, this->nkb * max_count);
    // transfer the gcar and gcar_zero_indexes to the device
    syncmem_var_h2d_op()(this->ctx, this->cpu_ctx, gcar, gcar_tmp.data(), 3 * npw_max);
    syncmem_int_h2d_op()(this->ctx, this->cpu_ctx, gcar_zero_indexes, gcar_zero_indexes_tmp.data(), 3 * npw_max);
}

// cal_force
template <typename FPTYPE, typename Device>
void FS_Nonlocal_tools<FPTYPE, Device>::cal_force(const int& ik,
                                                  const int& nbdall,
                                                  const int& npm,
                                                  const bool& occ,
                                                  FPTYPE* force,
                                                  const int& nbd0)
{
    const int current_spin = this->kv_->isk[ik];
    const int force_nc = 3;
    const int npol = this->ucell_->get_npol();
    const int index0 = nbd0 * npol * this->nkb;
    // calculate the force
    if (npol == 1)
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
        
        cal_force_nl_op<FPTYPE, Device>()(this->ctx,
                                          nondiagonal,
                                          npm,
                                          this->ntype,
                                          current_spin,
                                          this->nlpp_->deeq.getBound2(),
                                          this->nlpp_->deeq.getBound3(),
                                          this->nlpp_->deeq.getBound4(),
                                          force_nc,
                                          nbdall,
                                          nkb,
                                          atom_nh,
                                          atom_na,
                                          this->ucell_->tpiba,
                                          d_wg_ik,
                                          occ,
                                          d_ekb_ik,
                                          qq_nt,
                                          deeq,
                                          becp + index0,
                                          dbecp + index0,
                                          force);
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
        cal_force_nl_op<FPTYPE, Device>()(this->ctx,
                                          npm,
                                          this->ntype,
                                          this->nlpp_->deeq_nc.getBound2(),
                                          this->nlpp_->deeq_nc.getBound3(),
                                          this->nlpp_->deeq_nc.getBound4(),
                                          force_nc,
                                          nbdall,
                                          nkb,
                                          atom_nh,
                                          atom_na,
                                          this->ucell_->tpiba,
                                          d_wg_ik,
                                          occ,
                                          d_ekb_ik,
                                          qq_nt,
                                          this->nlpp_->template get_deeq_nc_data<FPTYPE>(),
                                          becp + index0,
                                          dbecp + index0,
                                          force);
    }
}

// template instantiation
template class FS_Nonlocal_tools<double, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class FS_Nonlocal_tools<double, base_device::DEVICE_GPU>;
#endif

} // namespace hamilt
