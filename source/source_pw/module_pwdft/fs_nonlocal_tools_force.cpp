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

template <typename FPTYPE, typename Device>
void FS_Nonlocal_tools<FPTYPE, Device>::cal_vkb_deri_f(const int& ik, const int& nbdall, const int& ipol)
{
    ModuleBase::TITLE("FS_Nonlocal_tools", "cal_vkb_deri");
    const int npol = this->ucell_->get_npol();
    const int size_becp = nbdall * npol * this->nkb;
    if (this->dbecp == nullptr)
    {
        resmem_complex_op()(dbecp, 3 * size_becp);
    }

    const std::complex<FPTYPE>* vkb_ptr = this->ppcell_vkb;
    std::complex<FPTYPE>* vkb_deri_ptr = this->ppcell_vkb;

    const int npw = this->wfc_basis_->npwk[ik];
    if (this->pre_ik_f == -1)
    {
        resmem_var_op()(gcar, 3 * this->wfc_basis_->npwk_max);
        resmem_int_op()(gcar_zero_indexes, 3 * this->wfc_basis_->npwk_max);
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
    resmem_complex_op()(this->vkb_save, this->nkb * max_count);
    // transfer the gcar and gcar_zero_indexes to the device
    syncmem_var_h2d_op()(gcar, gcar_tmp.data(), 3 * npw_max);
    syncmem_int_h2d_op()(gcar_zero_indexes, gcar_zero_indexes_tmp.data(), 3 * npw_max);
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
template void FS_Nonlocal_tools<double, base_device::DEVICE_CPU>::cal_vkb_deri_f(const int&,
                                                                                  const int&,
                                                                                  const int&);
template void FS_Nonlocal_tools<double, base_device::DEVICE_CPU>::cal_dbecp_f(const int&,
                                                                               const int&,
                                                                               const int&,
                                                                               const int&,
                                                                               const std::complex<double>*,
                                                                               const int&);
template void FS_Nonlocal_tools<double, base_device::DEVICE_CPU>::save_vkb(const int&, const int&);
template void FS_Nonlocal_tools<double, base_device::DEVICE_CPU>::revert_vkb(const int&, const int&);
template void FS_Nonlocal_tools<double, base_device::DEVICE_CPU>::transfer_gcar(const int&,
                                                                                const int&,
                                                                                const double*);
template void FS_Nonlocal_tools<double, base_device::DEVICE_CPU>::cal_force(const int&,
                                                                            const int&,
                                                                            const int&,
                                                                            const bool&,
                                                                            double*,
                                                                            const int&);
#if ((defined __CUDA) || (defined __ROCM))
template void FS_Nonlocal_tools<double, base_device::DEVICE_GPU>::cal_vkb_deri_f(const int&,
                                                                                  const int&,
                                                                                  const int&);
template void FS_Nonlocal_tools<double, base_device::DEVICE_GPU>::cal_dbecp_f(const int&,
                                                                               const int&,
                                                                               const int&,
                                                                               const int&,
                                                                               const std::complex<double>*,
                                                                               const int&);
template void FS_Nonlocal_tools<double, base_device::DEVICE_GPU>::save_vkb(const int&, const int&);
template void FS_Nonlocal_tools<double, base_device::DEVICE_GPU>::revert_vkb(const int&, const int&);
template void FS_Nonlocal_tools<double, base_device::DEVICE_GPU>::transfer_gcar(const int&,
                                                                                const int&,
                                                                                const double*);
template void FS_Nonlocal_tools<double, base_device::DEVICE_GPU>::cal_force(const int&,
                                                                            const int&,
                                                                            const int&,
                                                                            const bool&,
                                                                            double*,
                                                                            const int&);
#endif

} // namespace hamilt
