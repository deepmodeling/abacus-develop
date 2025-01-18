#include "onsite_proj_tools.h"

#include "module_base/math_polyint.h"
#include "module_base/math_ylmreal.h"
#include "module_base/memory.h"
#include "module_base/timer.h"
#include "module_base/tool_title.h"
#include "module_hamilt_pw/hamilt_pwdft/kernels/force_op.h"
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
    this->device = base_device::get_device_type<Device>(this->ctx);

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

    this->device = base_device::get_device_type<Device>(this->ctx);

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

// cal_becp
// starts from vkb (nkb, ng) table
// it should be merely the multiplication of matrix (vkb, ng) * (ng, nbands) -> (vkb, nbands)
// should be irrelevant with what the matrix is.
// the vkb index generation should be maintained elsewhere.
// vkb already has atomic position information, calculated from the vq and sk
// . the multiplication with sk should be within specific operator
// because the atom selection task is operator-specific.
template <typename FPTYPE, typename Device>
void Onsite_Proj_tools<FPTYPE, Device>::cal_becp(int ik,
                                                 int npm,
                                                 std::complex<FPTYPE>* becp_in,
                                                 const std::complex<FPTYPE>* ppsi_in)
{
    ModuleBase::TITLE("Onsite_Proj_tools", "cal_becp");
    ModuleBase::timer::tick("Onsite_Proj_tools", "cal_becp");

    const int npol = this->ucell_->get_npol();
    const std::complex<FPTYPE>* ppsi = ppsi_in == nullptr ? &(this->psi_[0](ik, 0, 0)) : ppsi_in;
    const int npw = this->wfc_basis_->npwk[ik];
    if (becp_in == nullptr && this->becp == nullptr)
    {
        resmem_complex_op()(becp, this->nbands * npol * this->nkb);
    }
    std::complex<FPTYPE>* becp_tmp = becp_in == nullptr ? this->becp : becp_in;
    const int size_becp_act = npm * npol * this->nkb;
    if (ik != this->current_ik) // different ik, need to recalculate vkb
    {
        const int size_becp = this->nbands * npol * this->nkb;
        if (this->becp == nullptr)
        {
            resmem_complex_op()(becp, size_becp);
        }

        // prepare math tools
        Nonlocal_maths<FPTYPE, Device> maths(*(this->nhtol), this->lprojmax, this->ucell_);

        std::complex<FPTYPE>* vkb_ptr = this->ppcell_vkb;

        // calculate G+K
        this->g_plus_k = maths.cal_gk(ik, this->wfc_basis_);
        FPTYPE *gk = g_plus_k.data(), *vq_tb = this->tabtpr->ptr;
        // vq_tb has dimension (ntype, nproj, GlobalV::NQX)

        // calculate sk
        resmem_complex_op()(hd_sk, this->ucell_->nat * npw);
        this->sf_->get_sk(ctx, ik, this->wfc_basis_, hd_sk);
        std::complex<FPTYPE>* d_sk = this->hd_sk;
        // prepare ylm，size: (lmax+1)^2 * this->max_npw
        const int lmax_ = this->lprojmax;
        maths.cal_ylm(lmax_, npw, g_plus_k.data(), hd_ylm);

        // DEBUG: ONCE YOU CHECK ylm VALUES, YOU UNCOMMENT THE FOLLOW
        // std::vector<ModuleBase::Vector3<double>> qs(npw);
        // for (int ig = 0; ig < npw; ig++)
        // {
        //     qs[ig] = this->wfc_basis_->getgpluskcar(ik, ig);
        // }
        // const int total_lm = (lmax_ + 1) * (lmax_ + 1);
        // ModuleBase::matrix ylmref(total_lm, npw);
        // ModuleBase::YlmReal::Ylm_Real(total_lm, npw, qs.data(), ylmref);
        // std::cout << "Compare the Ylm values of two methods:" << std::endl;
        // int lm = 0;
        // for(int l_ = 0; l_ < lmax_ + 1; l_++)
        // {
        //     for(int m_ = -l_; m_ <= l_; m_++)
        //     {
        //         std::cout << "l = " << l_ << " m = " << m_ << std::endl;
        //         lm = l_ * l_ + l_ + m_;
        //         for(int ig = 0; ig < npw; ig++)
        //         {
        //             std::cout << "[" << ylmref(lm, ig) << " " << hd_ylm[lm * npw + ig] << "]" << " ";
        //         }
        //         std::cout << std::endl;
        //     }
        //     std::cout << std::endl;
        // }
        // ModuleBase::WARNING_QUIT("Onsite_Proj_tools", "cal_becp");

        if (this->device == base_device::GpuDevice)
        {
            syncmem_var_h2d_op()(d_g_plus_k, g_plus_k.data(), g_plus_k.size());
            syncmem_var_h2d_op()(d_vq_tab, this->tabtpr->ptr, this->tabtpr->getSize());
            gk = d_g_plus_k;
            vq_tb = d_vq_tab;
        }

        // int vkb_size = 0;
        for (int it = 0; it < this->ucell_->ntype; it++) // loop all elements
        {
            // interpolate (it, 0..nproj[it], 0..npw) to get hd_vq
            cal_vq_op()(this->ctx,
                        vq_tb, // its data is correct, dimension (ntype, nprojmax, GlobalV::NQX)
                        it,    // but maybe it is (ntype, nprojmax*npol, GlobalV::NQX)
                        gk,
                        npw,
                        this->tabtpr->getBound2(),
                        this->tabtpr->getBound3(),
                        PARAM.globalv.dq,
                        nproj[it],
                        hd_vq); // hd_vq has dimension (nprojmax, npwx), this size will be the largest needed.

            // DEBUG: ONCE YOU CHECK vq VALUES, YOU UNCOMMENT THE FOLLOWING LINE
            // for(int ip = 0; ip < nproj[it]; ip++)
            // {
            //     std::cout << "projector #" << ip << " of atom type " << it << std::endl;
            //     for(int iq = 0; iq < npw; iq++)
            //     {
            //         std::cout << hd_vq[ip * npw + iq] << " ";
            //     }
            //     std::cout << std::endl;
            // }
            // std::cout << std::endl;

            // prepare（-i）^l, size: nh
            std::vector<std::complex<double>> pref = maths.cal_pref(it, h_atom_nh[it]);
            const int nh = pref.size();
            this->dvkb_indexes.resize(nh * 4);
            // print the value of nhtol
            // nhtol->print(std::cout); // as checked, nhtol works as expected
            maths.cal_dvkb_index(nproj[it], this->nhtol->c, this->nhtol->nc, npw, it, 0, 0, this->dvkb_indexes.data());

            if (this->device == base_device::GpuDevice)
            {
                syncmem_int_h2d_op()(d_dvkb_indexes, dvkb_indexes.data(), nh * 4);
                syncmem_complex_h2d_op()(d_pref_in, pref.data(), nh);
            }

            for (int ia = 0; ia < h_atom_na[it]; ia++)
            {
                if (this->device == base_device::CpuDevice)
                {
                    d_pref_in = pref.data();
                    d_dvkb_indexes = dvkb_indexes.data();
                }
                cal_vkb_op()(this->ctx, nh, npw, d_dvkb_indexes, hd_vq, hd_ylm, d_sk, d_pref_in, vkb_ptr);
                vkb_ptr += nh * npw; // vkb_ptr has dimension (nhtot, npwx), this size will be the largest needed.
                d_sk += npw;
                // vkb_size += nh * npw;
            }
        }
        this->current_ik = ik;
    }
    // DEBUG: ONCE YOU CHECK vkb VALUES, YOU UNCOMMENT THE FOLLOWING LINE
    // for(int i = 0; i < vkb_size; i++)
    // {
    //     if (i % npw == 0)
    //     {
    //         std::cout << "The #" << i / npw << " projector" << std::endl;
    //     }
    //     std::cout << this->ppcell_vkb[i] << " ";
    // }
    // std::cout << std::endl;
    // ModuleBase::WARNING_QUIT("Onsite_Proj_tools", "cal_becp");

    // PLAN: seperate the lower and upper into two parts, individually called.
    const char transa = 'C';
    const char transb = 'N';
    int npm_npol = npm * npol;
    gemm_op()(this->ctx,
              transa,
              transb,
              this->nkb,
              npm_npol, // nbands(occ)*npol
              npw,
              &ModuleBase::ONE,
              this->ppcell_vkb,
              npw,
              ppsi,
              this->max_npw,
              &ModuleBase::ZERO,
              becp_tmp,
              this->nkb);

    if (this->device == base_device::GpuDevice)
    {
        std::complex<FPTYPE>* h_becp = nullptr;
        resmem_complex_h_op()(h_becp, size_becp_act);
        syncmem_complex_d2h_op()(h_becp, becp_tmp, size_becp_act);
        Parallel_Reduce::reduce_pool(h_becp, size_becp_act);
        syncmem_complex_h2d_op()(becp_tmp, h_becp, size_becp_act);
        delmem_complex_h_op()(h_becp);
    }
    else
    {
        Parallel_Reduce::reduce_pool(becp_tmp, size_becp_act);
    }
    // DEBUG: ONCE YOU CHECK becp VALUES, YOU UNCOMMENT THE FOLLOWING LINE
    // std::cout << "ik: " << ik << std::endl;
    // for (int i = 0; i < npm_npol*this->nkb; i++)
    // {
    //     std::cout << "becp[" << i << "]: " << becp[i] << std::endl;
    // }
    ModuleBase::timer::tick("Onsite_Proj_tools", "cal_becp");
}

// cal_dbecp
template <typename FPTYPE, typename Device>
void Onsite_Proj_tools<FPTYPE, Device>::cal_dbecp_s(int ik, int npm, int ipol, int jpol)
{
    ModuleBase::TITLE("Onsite_Proj_tools", "cal_dbecp_s");
    ModuleBase::timer::tick("Onsite_Proj_tools", "cal_dbecp_s");
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

    gemm_op()(this->ctx,
              transa,
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
    ModuleBase::timer::tick("Onsite_Proj_tools", "cal_dbecp_s");
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
    ModuleBase::timer::tick("Onsite_Proj_tools", "cal_dbecp_f");

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
              nkb);
    this->revert_vkb(npw, ipol);
    this->pre_ik_f = ik;
    ModuleBase::timer::tick("Onsite_Proj_tools", "cal_dbecp_f");
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

template <typename FPTYPE, typename Device>
void Onsite_Proj_tools<FPTYPE, Device>::cal_force_dftu(int ik,
                                                       int npm,
                                                       FPTYPE* force,
                                                       const int* orbital_corr,
                                                       const std::complex<FPTYPE>* vu,
                                                       const int size_vu,
                                                       const FPTYPE* h_wg)
{
    int* orbital_corr_tmp = nullptr;
    std::complex<FPTYPE>* vu_tmp = nullptr;
#if defined(__CUDA) || defined(__ROCM)
    if (this->device == base_device::GpuDevice)
    {
        resmem_int_op()(orbital_corr_tmp, this->ucell_->ntype);
        syncmem_int_h2d_op()(orbital_corr_tmp, orbital_corr, this->ucell_->ntype);
        resmem_complex_op()(vu_tmp, size_vu);
        syncmem_complex_h2d_op()(vu_tmp, vu, size_vu);
        syncmem_var_h2d_op()(d_wg, h_wg, this->nbands * (ik+1));
    }
    else
#endif
    {
        orbital_corr_tmp = const_cast<int*>(orbital_corr);
        vu_tmp = const_cast<std::complex<FPTYPE>*>(vu);
        d_wg = const_cast<FPTYPE*>(h_wg);
    }
    const int force_nc = 3;
    cal_force_nl_op<FPTYPE, Device>()(this->ctx,
                                      npm,
                                      this->nbands,
                                      this->ntype,
                                      force_nc,
                                      this->nbands,
                                      ik,
                                      nkb,
                                      atom_nh,
                                      atom_na,
                                      this->ucell_->tpiba,
                                      d_wg,
                                      vu_tmp,
                                      orbital_corr_tmp,
                                      becp,
                                      dbecp,
                                      force);
#if defined(__CUDA) || defined(__ROCM)
    if (this->device == base_device::GpuDevice)
    {
        delmem_complex_op()(vu_tmp);
        delmem_int_op()(orbital_corr_tmp);
    }
#endif
}

template <typename FPTYPE, typename Device>
void Onsite_Proj_tools<FPTYPE, Device>::cal_force_dspin(int ik,
                                                        int npm,
                                                        FPTYPE* force,
                                                        const ModuleBase::Vector3<double>* lambda,
                                                        const FPTYPE* h_wg)
{
    std::vector<FPTYPE> lambda_array(this->ucell_->nat * 3);
    for (int iat = 0; iat < this->ucell_->nat; iat++)
    {
        lambda_array[iat * 3] = lambda[iat].x;
        lambda_array[iat * 3 + 1] = lambda[iat].y;
        lambda_array[iat * 3 + 2] = lambda[iat].z;
    }
    FPTYPE* lambda_tmp = nullptr;
#if defined(__CUDA) || defined(__ROCM)
    if (this->device == base_device::GpuDevice)
    {
        resmem_var_op()(lambda_tmp, this->ucell_->nat * 3);
        syncmem_var_h2d_op()(lambda_tmp, lambda_array.data(), this->ucell_->nat * 3);
        syncmem_var_h2d_op()(d_wg, h_wg, this->nbands * (ik+1));
    }
    else
#endif
    {
        lambda_tmp = lambda_array.data();
        d_wg = const_cast<FPTYPE*>(h_wg);
    }
    const int force_nc = 3;
    cal_force_nl_op<FPTYPE, Device>()(this->ctx,
                                      npm,
                                      this->nbands,
                                      this->ntype,
                                      force_nc,
                                      this->nbands,
                                      ik,
                                      nkb,
                                      atom_nh,
                                      atom_na,
                                      this->ucell_->tpiba,
                                      d_wg,
                                      lambda_tmp,
                                      becp,
                                      dbecp,
                                      force);

#if defined(__CUDA) || defined(__ROCM)
    if (this->device == base_device::GpuDevice)
    {
        delmem_var_op()(lambda_tmp);
    }
#endif
}

template <typename FPTYPE, typename Device>
void Onsite_Proj_tools<FPTYPE, Device>::cal_stress_dftu(int ik,
                                                        int npm,
                                                        FPTYPE* stress,
                                                        const int* orbital_corr,
                                                        const std::complex<FPTYPE>* vu,
                                                        const int size_vu,
                                                        const FPTYPE* h_wg)
{
    int* orbital_corr_tmp = nullptr;
    std::complex<FPTYPE>* vu_tmp = nullptr;
#if defined(__CUDA) || defined(__ROCM)
    if (this->device == base_device::GpuDevice)
    {
        resmem_int_op()(orbital_corr_tmp, this->ucell_->ntype);
        syncmem_int_h2d_op()(orbital_corr_tmp, orbital_corr, this->ucell_->ntype);
        resmem_complex_op()(vu_tmp, size_vu);
        syncmem_complex_h2d_op()(vu_tmp, vu, size_vu);
        syncmem_var_h2d_op()(d_wg, h_wg, this->nbands * (ik+1));
    }
    else
#endif
    {
        orbital_corr_tmp = const_cast<int*>(orbital_corr);
        vu_tmp = const_cast<std::complex<FPTYPE>*>(vu);
        d_wg = const_cast<FPTYPE*>(h_wg);
    }
    cal_stress_nl_op()(this->ctx,
                       nkb,
                       npm,
                       this->ntype,
                       this->nbands,
                       ik,
                       atom_nh,
                       atom_na,
                       d_wg,
                       vu_tmp,
                       orbital_corr_tmp,
                       becp,
                       dbecp,
                       stress);
#if defined(__CUDA) || defined(__ROCM)
    if (this->device == base_device::GpuDevice)
    {
        delmem_complex_op()(vu_tmp);
        delmem_int_op()(orbital_corr_tmp);
    }
#endif
}

template <typename FPTYPE, typename Device>
void Onsite_Proj_tools<FPTYPE, Device>::cal_stress_dspin(int ik,
                                                         int npm,
                                                         FPTYPE* stress,
                                                         const ModuleBase::Vector3<double>* lambda,
                                                         const FPTYPE* h_wg)
{
    std::vector<FPTYPE> lambda_array(this->ucell_->nat * 3);
    for (int iat = 0; iat < this->ucell_->nat; iat++)
    {
        lambda_array[iat * 3] = lambda[iat].x;
        lambda_array[iat * 3 + 1] = lambda[iat].y;
        lambda_array[iat * 3 + 2] = lambda[iat].z;
    }
    FPTYPE* lambda_tmp = nullptr;
#if defined(__CUDA) || defined(__ROCM)
    if (this->device == base_device::GpuDevice)
    {
        resmem_var_op()(lambda_tmp, this->ucell_->nat * 3);
        syncmem_var_h2d_op()(lambda_tmp, lambda_array.data(), this->ucell_->nat * 3);
        syncmem_var_h2d_op()(d_wg, h_wg, this->nbands * (ik+1));
    }
    else
#endif
    {
        lambda_tmp = lambda_array.data();
        d_wg = const_cast<FPTYPE*>(h_wg);
    }
    const int force_nc = 3;
    cal_stress_nl_op()(this->ctx,
                       nkb,
                       npm,
                       this->ntype,
                       this->nbands,
                       ik,
                       atom_nh,
                       atom_na,
                       d_wg,
                       lambda_tmp,
                       becp,
                       dbecp,
                       stress);

#if defined(__CUDA) || defined(__ROCM)
    if (this->device == base_device::GpuDevice)
    {
        delmem_var_op()(lambda_tmp);
    }
#endif
}

// template instantiation
template class Onsite_Proj_tools<double, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class Onsite_Proj_tools<double, base_device::DEVICE_GPU>;
#endif

} // namespace hamilt
