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
    this->device = base_device::get_device_type(this->ctx);
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
        resmem_int_op()(this->d_dvkb_indexes, max_nh * 4);
    }

    resmem_var_op()(this->hd_vq, max_nbeta * max_npw);
    resmem_var_op()(this->hd_vq_deri, max_nbeta * max_npw);
    const int _lmax = this->nlpp_->lmaxkb;
    resmem_var_op()(this->hd_ylm, (_lmax + 1) * (_lmax + 1) * max_npw);
    resmem_var_op()(this->hd_ylm_deri, 3 * (_lmax + 1) * (_lmax + 1) * max_npw);
    const int nks = this->kv_->get_nks();
    resmem_var_op()(d_wk, nks);
    syncmem_var_h2d_op()(d_wk, this->kv_->wk.data(), nks);

    if (this->device == base_device::GpuDevice)
    {
        resmem_var_op()(d_wg, wg.nr * wg.nc);
        syncmem_var_h2d_op()(d_wg, wg.c, wg.nr * wg.nc);
        if (p_ekb != nullptr)
        {
            resmem_var_op()(d_ekb, p_ekb->nr * p_ekb->nc);
            syncmem_var_h2d_op()(d_ekb, p_ekb->c, p_ekb->nr * p_ekb->nc);
        }
        resmem_int_op()(atom_nh, this->ntype);
        resmem_int_op()(atom_na, this->ntype);
        syncmem_int_h2d_op()(atom_nh, h_atom_nh.data(), this->ntype);
        syncmem_int_h2d_op()(atom_na, h_atom_na.data(), this->ntype);

        resmem_var_op()(d_g_plus_k, max_npw * 5);
        resmem_var_op()(d_pref, max_nh);
        resmem_var_op()(d_vq_tab, this->nlpp_->tab.getSize());
        resmem_complex_op()(d_pref_in, max_nh);

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

    delmem_var_op()(hd_vq);
    delmem_var_op()(hd_vq_deri);
    delmem_var_op()(hd_ylm);
    delmem_var_op()(hd_ylm_deri);
    delmem_var_op()(d_wk);

    // delete memory on GPU
    if (this->device == base_device::GpuDevice)
    {
        delmem_var_op()(d_wg);
        delmem_var_op()(d_ekb);
        delmem_int_op()(atom_nh);
        delmem_int_op()(atom_na);
        delmem_var_op()(d_g_plus_k);
        delmem_var_op()(d_pref);
        delmem_var_op()(d_vq_tab);
        delmem_complex_op()(this->d_pref_in);
        delmem_int_op()(d_dvkb_indexes);
    }

    if (becp != nullptr)
    {
        delmem_complex_op()(becp);
        delmem_complex_op()(hd_sk);
    }
    if (dbecp != nullptr)
    {
        delmem_complex_op()(dbecp);
    }
    if (this->pre_ik_f != -1)
    {
        delmem_int_op()(gcar_zero_indexes);
        delmem_complex_op()(vkb_save);
        delmem_var_op()(gcar);
    }
}

// template instantiation
template class FS_Nonlocal_tools<double, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class FS_Nonlocal_tools<double, base_device::DEVICE_GPU>;
#endif

} // namespace hamilt
