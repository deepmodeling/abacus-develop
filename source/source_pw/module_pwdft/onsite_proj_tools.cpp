#include "onsite_proj_tools.h"

#include "source_base/math_polyint.h"
#include "source_base/math_ylmreal.h"
#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include "source_base/tool_title.h"
#include "source_pw/module_pwdft/kernels/force_op.h"
#include "source_io/module_parameter/parameter.h"
#include "nonlocal_maths.hpp"

#include <numeric>

namespace hamilt
{
template <typename FPTYPE, typename Device>
Onsite_Proj_tools<FPTYPE, Device>::Onsite_Proj_tools(const pseudopot_cell_vnl* nlpp_in,
                                                     const UnitCell* ucell_in,
                                                     const psi::Psi<std::complex<FPTYPE>, Device>* psi_in,
                                                     const K_Vectors* kv_in,
                                                     const ModulePW::PW_Basis_K* wfc_basis_in,
                                                     const Structure_Factor* sf_in,
                                                     const ModuleBase::matrix& wg,
                                                     const ModuleBase::matrix& ekb)
    : nlpp_(nlpp_in), ucell_(ucell_in), psi_(psi_in), kv_(kv_in), wfc_basis_(wfc_basis_in), sf_(sf_in)
{
    // get the device context
    this->device = base_device::get_device_type(this->ctx);

    // seems kvec_c never used...
    this->kvec_c = this->wfc_basis_->template get_kvec_c_data<FPTYPE>();
    // the following is important for calculating the whole contribution to
    // Hamiltonian or force, stress: sum{nk} fnk*sum_{ij}<psi_nk|ai>Dij<aj|psi_nk>
    // among, Dij is deeq.
    // For DFT+U and other projection involved operators, deeq also plays.
    this->deeq = this->nlpp_->template get_deeq_data<FPTYPE>();
    this->deeq_dims[0] = this->nlpp_->deeq.getBound1();
    this->deeq_dims[1] = this->nlpp_->deeq.getBound2();
    this->deeq_dims[2] = this->nlpp_->deeq.getBound3();
    this->deeq_dims[3] = this->nlpp_->deeq.getBound4();
    this->deeq_nc = this->nlpp_->template get_deeq_nc_data<FPTYPE>();
    this->deeq_nc_dims[0] = this->nlpp_->deeq_nc.getBound1();
    this->deeq_nc_dims[1] = this->nlpp_->deeq_nc.getBound2();
    this->deeq_nc_dims[2] = this->nlpp_->deeq_nc.getBound3();
    this->deeq_nc_dims[3] = this->nlpp_->deeq_nc.getBound4();
    // ultrasoft pseudopotential
    this->qq_nt = this->nlpp_->template get_qq_nt_data<FPTYPE>();
    // total number of projectors (all types, all atoms, not m-distinguishive)
    this->nkb = nlpp_->nkb;
    // not clear why do these following...
    this->nbands = psi_->get_nbands();
    this->max_npw = wfc_basis_->npwk_max;
    this->ntype = ucell_->ntype;
    // because the code is needed to reuse, therefore all other parts should be general
    // and not strongly depend on any structure of class pseudopot_cell_vnl, therefore
    // here unpack all needed information.
    this->tabtpr = &(nlpp_->tab);
    this->nhtol = &(nlpp_->nhtol);
    this->lprojmax = nlpp_->lmaxkb;
    // There is a contribution for jh<>ih in US case or multi projectors case
    // Actually, the judge of nondiagonal should be done on every atom type
    this->nondiagonal = (PARAM.globalv.use_uspp || this->nlpp_->multi_proj) ? true : false;

    this->nproj.resize(this->ntype);
    std::vector<int> nch(this->ntype);
    for (int it = 0; it < this->ntype; it++)
    {
        this->nproj[it] = this->ucell_->atoms[it].ncpp.nbeta;
        nch[it] = this->ucell_->atoms[it].ncpp.nh;
    }
    // allocate memory
    this->allocate_memory(wg, ekb, this->nproj, nch);
    this->ppcell_vkb
        = (this->device == base_device::GpuDevice) ? this->nlpp_->template get_vkb_data<FPTYPE>() : this->nlpp_->vkb.c;
}

template <typename FPTYPE, typename Device>
Onsite_Proj_tools<FPTYPE, Device>::Onsite_Proj_tools(
    const std::vector<int>& nproj, // number of projectors for each atom type
    const std::vector<int>& lproj,
    const ModuleBase::realArray& tab, // radials' spherical bessel transform
    const ModuleBase::matrix& nhtol,  // (it, ich) -> l, the ich is (l, m)-distinctive index
    std::complex<FPTYPE>* vkb_buf,    // the vkb buffer
    const UnitCell* ucell_in,
    const psi::Psi<std::complex<FPTYPE>, Device>* psi_in,
    const K_Vectors* kv_in,
    const ModulePW::PW_Basis_K* wfc_basis_in,
    const Structure_Factor* sf_in,
    const ModuleBase::matrix& wg,
    const ModuleBase::matrix& ekb)
{
    // this is a constructor for general case, including vnl, dftu, deltaspin, deepks, etc.
    // what is needed for this kind of constructor?

    // ntype: from unitcell
    // nproj: number of projectors own by each atom type
    // projs: beta function or radial function
    // lproj: angular momentum of projectors
    // rgrid: radial grid
    // deeq: the Dij matrix, Hubbard parameters or other quantities...

    // what are already programmed to be needed?

    // tab: the spherical transform of radial functions, with q = linspace(0, GlobalV::NQX, GlobalV::DQ)
    // nhtol: the (it, ich) -> l, the ich is (l, m)-distinctive index
    // nkb: total # of projectors <- std::accumulate(nproj.begin(), nproj.end(), 0)
    // atom_nh: # of (l, m)-distinctive projectors for each atom type
    // h_atom_nh: counterpart of atom_nh on host
    // max_nh: std::max_element(atom_nh.begin(), atom_nh.end())

    // in conclusion, this constructor needs the following individual information:

    // nproj
    // tab (projs is not needed, should be calculated elsewhere)
    // lproj
    // deeq, with its dims. it will be good to pass the whole realarray

    // what can be built here
    // nhtol
    // nkb
    // atom_nh, h_atom_nh, max_nh
    // deeq_dims

    ucell_ = ucell_in;
    psi_ = psi_in;
    kv_ = kv_in;
    wfc_basis_ = wfc_basis_in;
    sf_ = sf_in;

    this->device = base_device::get_device_type(this->ctx);

    this->kvec_c = this->wfc_basis_->template get_kvec_c_data<FPTYPE>();
    // skip deeq, qq_nt
    this->nbands = psi_->get_nbands();
    this->max_npw = wfc_basis_->npwk_max;
    this->ntype = nproj.size();
    this->tabtpr = &tab;

    this->nhtol = &nhtol;
    this->lprojmax = *std::max_element(lproj.begin(), lproj.end());
    this->nondiagonal = false;

    this->nkb = 0;
    this->h_atom_nh.resize(this->ntype, 0);
    int iproj = 0;
    for (int it = 0; it < this->ntype; it++)
    {
        int nproj_it = nproj[it];
        for (int ip = 0; ip < nproj_it; ip++)
        {
            this->h_atom_nh[it] += 2 * lproj[iproj] + 1;
            this->nkb += (2 * lproj[iproj] + 1) * this->ucell_->atoms[it].na;
            iproj++;
        }
    }
    this->nproj = nproj;
    this->allocate_memory(wg, ekb, nproj, this->h_atom_nh);
    // what is this??? seems it is not on gpu
    this->ppcell_vkb = vkb_buf;
}

template <typename FPTYPE, typename Device>
Onsite_Proj_tools<FPTYPE, Device>::~Onsite_Proj_tools()
{
    // delete memory
    delete_memory();
}

template <typename FPTYPE, typename Device>
void Onsite_Proj_tools<FPTYPE, Device>::allocate_memory(const ModuleBase::matrix& wg,
                                                        const ModuleBase::matrix& ekb,
                                                        const std::vector<int>& nproj,
                                                        const std::vector<int>& nch)
{
    // allocate memory

    // prepare the memory of stress and init some variables:
    this->h_atom_nh.resize(this->ntype);
    this->h_atom_na.resize(this->ntype);
    for (int it = 0; it < this->ntype; it++)
    {
        h_atom_nh[it] = nch[it];
        h_atom_na[it] = this->ucell_->atoms[it].na;
    }

    int nprojmax = 0;
    for (int it = 0; it < this->ntype; it++) // loop all elements
    {
        nprojmax = std::max(nproj[it], nprojmax); // 0000000000000000000000000
        this->max_nh = std::max(h_atom_nh[it], max_nh);
    }

    // allocate the memory for vkb and vkb_deri.
    if (this->device == base_device::GpuDevice)
    {
        resmem_int_op()(this->d_dvkb_indexes, max_nh * 4);
    }

    resmem_var_op()(this->hd_vq, nprojmax * max_npw);
    resmem_var_op()(this->hd_vq_deri, nprojmax * max_npw);
    resmem_var_op()(this->hd_ylm, (lprojmax + 1) * (lprojmax + 1) * max_npw);
    resmem_var_op()(this->hd_ylm_deri, 3 * (lprojmax + 1) * (lprojmax + 1) * max_npw);

    if (this->device == base_device::GpuDevice)
    {
        resmem_var_op()(d_wg, wg.nr * wg.nc);
        resmem_var_op()(d_ekb, ekb.nr * ekb.nc);
        syncmem_var_h2d_op()(d_wg, wg.c, wg.nr * wg.nc);
        syncmem_var_h2d_op()(d_ekb, ekb.c, ekb.nr * ekb.nc);
        resmem_int_op()(atom_nh, this->ntype);
        resmem_int_op()(atom_na, this->ntype);
        syncmem_int_h2d_op()(atom_nh, h_atom_nh.data(), this->ntype);
        syncmem_int_h2d_op()(atom_na, h_atom_na.data(), this->ntype);

        resmem_var_op()(d_g_plus_k, max_npw * 5);
        resmem_var_op()(d_pref, max_nh);
        resmem_var_op()(d_vq_tab, this->tabtpr->getSize());
        resmem_complex_op()(d_pref_in, max_nh);
    }
    else
    {
        this->d_wg = wg.c;
        this->d_ekb = ekb.c;
        this->atom_nh = h_atom_nh.data();
        this->atom_na = h_atom_na.data();
    }
}

template <typename FPTYPE, typename Device>
void Onsite_Proj_tools<FPTYPE, Device>::delete_memory()
{
    // delete memory

    delmem_var_op()(hd_vq);
    delmem_var_op()(hd_vq_deri);
    delmem_var_op()(hd_ylm);
    delmem_var_op()(hd_ylm_deri);

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
template class Onsite_Proj_tools<double, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class Onsite_Proj_tools<double, base_device::DEVICE_GPU>;
#endif

} // namespace hamilt
