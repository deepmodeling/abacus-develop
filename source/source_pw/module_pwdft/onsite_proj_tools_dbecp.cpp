#include "onsite_proj_tools.h"

#include "source_base/math_polyint.h"
#include "source_base/math_ylmreal.h"
#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include "source_base/tool_title.h"
#include "source_pw/module_pwdft/kernels/force_op.h"
#include "source_io/module_parameter/parameter.h"
#include "nonlocal_maths.hpp"

namespace hamilt
{

// cal_dbecp
template <typename FPTYPE, typename Device>
void Onsite_Proj_tools<FPTYPE, Device>::cal_dbecp_s(int ik, int npm, int ipol, int jpol)
{
    ModuleBase::TITLE("Onsite_Proj_tools", "cal_dbecp_s");
    ModuleBase::timer::start("Onsite_Proj_tools", "cal_dbecp_s");
    this->current_ik = -1; // reset the current ik, vkb has been reused to save dvkb
    const int npol = this->ucell_->get_npol();
    const int size_becp = this->nbands * npol * this->nkb;
    const int npm_npol = npm * npol;
    if (this->dbecp == nullptr)
    {
        resmem_complex_op()(dbecp, size_becp);
    }

    // prepare math tools
    Nonlocal_maths<FPTYPE, Device> maths(*(this->nhtol), this->lprojmax, this->ucell_);

    const std::complex<FPTYPE>* ppsi = &(this->psi_[0](ik, 0, 0));
    const int npw = this->wfc_basis_->npwk[ik];
    std::complex<FPTYPE>* vkb_deri_ptr = this->ppcell_vkb;

    if (this->pre_ik_s != ik)
    { // k point has changed, we need to recalculate the g_plus_k
        // this->g_plus_k = maths.cal_gk(ik, this->wfc_basis_); //has been calculated by cal_becp

        const int lmax_ = this->lprojmax;
        // prepare ylm，size: (lmax+1)^2 * this->max_npw
        // maths.cal_ylm(lmax_, npw, g_plus_k.data(), hd_ylm); //has been calculated by cal_becp
        maths.cal_ylm_deri(lmax_, npw, g_plus_k.data(), hd_ylm_deri);
        this->pre_ik_s = ik;
    }
    FPTYPE *gk = g_plus_k.data(), *vq_tb = this->tabtpr->ptr;
    std::complex<FPTYPE>* d_sk = this->hd_sk;
    if (this->device == base_device::GpuDevice)
    {
        gk = d_g_plus_k;
        vq_tb = d_vq_tab;
    }

    for (int it = 0; it < this->ucell_->ntype; it++) // loop all elements
    {
        cal_vq_op()(this->ctx,
                    vq_tb,
                    it,
                    gk,
                    npw,
                    this->tabtpr->getBound2(),
                    this->tabtpr->getBound3(),
                    PARAM.globalv.dq,
                    this->nproj[it],
                    hd_vq);
        cal_vq_deri_op()(this->ctx,
                         vq_tb,
                         it,
                         gk,
                         npw,
                         this->tabtpr->getBound2(),
                         this->tabtpr->getBound3(),
                         PARAM.globalv.dq,
                         this->nproj[it],
                         hd_vq_deri);

        // prepare（-i）^l, size: nh
        std::vector<std::complex<double>> pref = maths.cal_pref(it, h_atom_nh[it]);
        int nh = pref.size();
        // prepare indexes for calculate vkb_deri
        this->dvkb_indexes.resize(nh * 4);
        maths.cal_dvkb_index(this->nproj[it],
                             this->nhtol->c,
                             this->nhtol->nc,
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
    // 2.b calculate dbecp = dbecp_noevc * psi
    const char transa = 'C';
    const char transb = 'N';

    gemm_op()(transa,
              transb,
              nkb,
              npm_npol,
              npw,
              &ModuleBase::ONE,
              ppcell_vkb,
              npw,
              ppsi,
              this->max_npw,
              &ModuleBase::ZERO,
              dbecp,
              nkb);
    ModuleBase::timer::end("Onsite_Proj_tools", "cal_dbecp_s");
}

// cal_dbecp_f
// starts from vkb (nkb, ng) table
// it should be again merely the multiplication of matrix (vkb, ng) * (ng, nbands) -> (vkb, nbands)
// the vkb is backed-up, and the memory space is reused for calculate ONE COMPONENT of dbecp
// . the direction of force is indexed by ipol (for stress, there are two, ipol and jpol).
// the dbecp_f is simply the becp multiplied with -i(G+k)_i
template <typename FPTYPE, typename Device>
void Onsite_Proj_tools<FPTYPE, Device>::cal_dbecp_f(int ik, int npm, int ipol)
{
    ModuleBase::TITLE("Onsite_Proj_tools", "cal_dbecp_f");
    ModuleBase::timer::start("Onsite_Proj_tools", "cal_dbecp_f");

    this->current_ik = -1; // reset the current ik, vkb has been reused to save dvkb

    const int npw = this->wfc_basis_->npwk[ik];

    // STAGE1: calculate dvkb_f
    // calculate gcarx, gcary/gcarx and gcarz/gcary, overwrite gcar
    if (this->pre_ik_f == -1) // if it is the very first run, we allocate
    {
        resmem_var_op()(gcar, 3 * this->wfc_basis_->npwk_max);
        resmem_int_op()(gcar_zero_indexes, 3 * this->wfc_basis_->npwk_max);
    }
    // first refresh the value of gcar_zero_indexes, gcar_zero_counts
    if (this->pre_ik_f != ik)
    { // the following lines will cause UNDEFINED BEHAVIOR because memory layout of vector3 instance
      // is assumed to be always contiguous but it is not guaranteed.
        this->transfer_gcar(npw,
                            this->wfc_basis_->npwk_max,
                            &(this->wfc_basis_->gcar[ik * this->wfc_basis_->npwk_max].x));
    }

    // backup vkb values to vkb_save
    this->save_vkb(npw, ipol);
    // for x, the coef is -i, for y and z it is 1
    const std::complex<double> coeff = ipol == 0 ? ModuleBase::NEG_IMAG_UNIT : ModuleBase::ONE;

    const std::complex<FPTYPE>* vkb_ptr = this->ppcell_vkb;
    std::complex<FPTYPE>* vkb_deri_ptr = this->ppcell_vkb;
    // calculate the vkb_deri for ipol with the memory of ppcell_vkb
    cal_vkb1_nl_op<FPTYPE, Device>()(this->ctx, nkb, npw, npw, npw, ipol, coeff, vkb_ptr, gcar, vkb_deri_ptr);

    // ------------------------------------------------------------------------------->8

    // STAGE2: calculate dbecp_f
    // NPOL
    // either 1 or 2, for NSPIN 1, 2 or 4 calculation
    // once NSPIN 4, there are doubled number of pw in each "row" of psi
    // on the other hand, for NSPIN 4 calculation, the number of bands is also doubled
    const int npol = this->ucell_->get_npol();
    const int npm_npol = npm * npol;
    const int size_becp = this->nbands * npol * this->nkb;
    if (this->dbecp == nullptr) // if it is the very first run, we allocate
    {                           // why not judging whether dbecp == nullptr inside resmem_complex_op?
        resmem_complex_op()(dbecp, 3 * size_becp);
    }
    // do gemm to get dbecp and revert the ppcell_vkb for next ipol
    const std::complex<FPTYPE>* ppsi = &(this->psi_[0](ik, 0, 0));
    // move the pointer to corresponding read&write position, according to ipol
    std::complex<FPTYPE>* dbecp_ptr = this->dbecp + ipol * size_becp; // [out]
    const char transa = 'C';
    const char transb = 'N';
    gemm_op()(transa,
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
              nkb);
    this->revert_vkb(npw, ipol);
    this->pre_ik_f = ik;
    ModuleBase::timer::end("Onsite_Proj_tools", "cal_dbecp_f");
}

// save_vkb
template <typename FPTYPE, typename Device>
void Onsite_Proj_tools<FPTYPE, Device>::save_vkb(int npw, int ipol)
{
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
void Onsite_Proj_tools<FPTYPE, Device>::revert_vkb(int npw, int ipol)
{
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
}

template <typename FPTYPE, typename Device>
void Onsite_Proj_tools<FPTYPE, Device>::transfer_gcar(int npw, int npw_max, const FPTYPE* gcar_in)
{
    std::vector<FPTYPE> gcar_tmp(3 * npw_max); // [out], will overwritten this->gcar
    gcar_tmp.assign(gcar_in,
                    gcar_in + 3 * npw_max); // UNDEFINED BEHAVIOR!!! nobody always knows the memory layout of vector3
    std::vector<int> gcar_zero_indexes_tmp(3 * npw_max); // a "checklist"

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
                ++gcar_zero_counts[i]; // num of zeros on each direction
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
    resmem_complex_op()(this->vkb_save, this->nkb * max_count);
    // transfer the gcar and gcar_zero_indexes to the device
    syncmem_var_h2d_op()(gcar, gcar_tmp.data(), 3 * npw_max);
    syncmem_int_h2d_op()(gcar_zero_indexes, gcar_zero_indexes_tmp.data(), 3 * npw_max);
}

// template instantiation
template void Onsite_Proj_tools<double, base_device::DEVICE_CPU>::cal_dbecp_s(int, int, int, int);
template void Onsite_Proj_tools<double, base_device::DEVICE_CPU>::cal_dbecp_f(int, int, int);
template void Onsite_Proj_tools<double, base_device::DEVICE_CPU>::save_vkb(int, int);
template void Onsite_Proj_tools<double, base_device::DEVICE_CPU>::revert_vkb(int, int);
template void Onsite_Proj_tools<double, base_device::DEVICE_CPU>::transfer_gcar(int, int, const double*);
#if ((defined __CUDA) || (defined __ROCM))
template void Onsite_Proj_tools<double, base_device::DEVICE_GPU>::cal_dbecp_s(int, int, int, int);
template void Onsite_Proj_tools<double, base_device::DEVICE_GPU>::cal_dbecp_f(int, int, int);
template void Onsite_Proj_tools<double, base_device::DEVICE_GPU>::save_vkb(int, int);
template void Onsite_Proj_tools<double, base_device::DEVICE_GPU>::revert_vkb(int, int);
template void Onsite_Proj_tools<double, base_device::DEVICE_GPU>::transfer_gcar(int, int, const double*);
#endif

} // namespace hamilt
