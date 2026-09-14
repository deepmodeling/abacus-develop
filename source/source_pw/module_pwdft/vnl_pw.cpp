#include "vnl_pw.h"

#include "source_io/module_parameter/parameter.h"
#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "source_base/output.h"
#include "source_base/memory_recorder.h"
#include "source_base/module_device/device.h"
#include "source_base/timer.h"

#include <cmath>

/**
 * @file vnl_pw.cpp
 * @brief Core resource management for pseudopot_cell_vnl.
 *
 * This file now keeps only the central lifecycle and accessor glue:
 * - ctor / dtor / release_memory(): device buffer cleanup
 * - init(): allocate and wire vkb, tab, tab_at, deeq, qq_nt/qq_so, nhtol etc.
 * - print_vnl(): dump tab to stream
 * - rescale_vnl(): rescale tab/tab_at/qrad when the cell volume changes
 * - get_*_data<T>() template specializations for CPU/GPU typed pointer access
 *
 * Heavy logic (getvnl, init_vnl, qrad, deeq, alpha channel) lives in:
 *   vnl_pw_getvnl.cpp, vnl_pw_init_vnl.cpp, vnl_pw_qrad.cpp,
 *   vnl_pw_deeq.cpp, vnl_pw_alpha.cpp, vnl_pw_grad.cpp
 */


pseudopot_cell_vnl::pseudopot_cell_vnl()
{
    this->use_gpu_ = (PARAM.inp.device == "gpu");
}

pseudopot_cell_vnl::~pseudopot_cell_vnl()
{
}

void pseudopot_cell_vnl::release_memory()
{
    if (this->nhm <= 0 || memory_released) {
        return;
}
    if (this->use_gpu_)
    {
        delmem_sd_op()(this->s_deeq);
        delmem_sd_op()(this->s_nhtol);
        delmem_sd_op()(this->s_nhtolm);
        delmem_sd_op()(this->s_indv);
        delmem_sd_op()(this->s_tab);
        delmem_sd_op()(this->s_qq_nt);
        delmem_cd_op()(this->c_deeq_nc);
        delmem_cd_op()(this->c_vkb);
        delmem_cd_op()(this->c_qq_so);
        delmem_zd_op()(this->z_deeq_nc);
        delmem_zd_op()(this->z_qq_so);
        delmem_dd_op()(this->d_deeq);
        delmem_zd_op()(this->z_vkb);
        delmem_dd_op()(this->d_tab);
        delmem_dd_op()(this->d_indv);
        delmem_dd_op()(this->d_nhtol);
        delmem_dd_op()(this->d_nhtolm);
        delmem_dd_op()(this->d_qq_nt);
    }
    else
    {
        delmem_sh_op()(this->s_deeq);
        delmem_sh_op()(this->s_nhtol);
        delmem_sh_op()(this->s_nhtolm);
        delmem_sh_op()(this->s_indv);
        delmem_sh_op()(this->s_tab);
        delmem_sh_op()(this->s_qq_nt);
        delmem_ch_op()(this->c_deeq_nc);
        delmem_ch_op()(this->c_vkb);
        delmem_ch_op()(this->c_qq_so);
#ifdef __DSP
        if (this->z_vkb != nullptr)
        {
            base_device::memory::delete_memory_op_mt<std::complex<double>, base_device::DEVICE_CPU>()(this->z_vkb);
            this->z_vkb = nullptr;
        }
#endif
        // There's no need to delete double precision pointers while in a CPU environment.
    }
    memory_released = true;
}

//-----------------------------------
// setup lmaxkb, nhm, nkb, lmaxq
// allocate vkb, PARAM.globalv.nqx, tab, tab_at
//-----------------------------------
void pseudopot_cell_vnl::init(const UnitCell& ucell,
                              Structure_Factor* psf_in,
                              const ModulePW::PW_Basis_K* wfc_basis,
                              const bool allocate_vkb)
{
    const int ntype = ucell.ntype;
    ModuleBase::TITLE("pseudopot_cell_vnl", "init");
    ModuleBase::timer::start("ppcell_vnl", "init");

    GlobalV::ofs_running << "\n SETUP NONLOCAL PSEUDOPOTENTIALS IN PLANE WAVE BASIS" << std::endl;

    int it = 0;
    this->wfcpw = wfc_basis;
    this->psf = psf_in;
    //----------------------------------------------------------
    // MEMBER VARIABLE :
    // NAME : lmaxkb(max angular momentum,(see pseudo_h))
    //----------------------------------------------------------
    this->lmaxkb = -1;
    for (it = 0; it < ntype; it++)
    {
        GlobalV::ofs_running << " " << ucell.atoms[it].label << " non-local projectors:" << std::endl;
        for (int ibeta = 0; ibeta < ucell.atoms[it].ncpp.nbeta; ibeta++)
        {
            GlobalV::ofs_running << " projector " << ibeta + 1 << " L=" << ucell.atoms[it].ncpp.lll[ibeta]
                                 << std::endl;
            this->lmaxkb = std::max(this->lmaxkb, ucell.atoms[it].ncpp.lll[ibeta]);
        }
    }

    //----------------------------------------------------------
    // MEMBER VARIABLE :
    // NAME : nhm(max number of different beta functions per atom)
    // NAME : nbetam(max number of beta functions)
    // NAME : nwfcm(max number of wavefunctions)
    //----------------------------------------------------------
    this->nhm = 0;
    this->nbetam = 0;
    int nwfcm = 0;
    for (it = 0; it < ntype; it++)
    {
        this->nhm = std::max(nhm, ucell.atoms[it].ncpp.nh);
        this->nbetam = std::max(nbetam, ucell.atoms[it].ncpp.nbeta);
        nwfcm = std::max(nwfcm, ucell.atoms[it].ncpp.nchi);
    }

    //----------------------------------------------------------
    // MEMBER VARIABLE :
    // NAME : nkb(total number of beta functions, with struct.fact.)
    //----------------------------------------------------------
    this->nkb = 0;
    for (it = 0; it < ntype; it++)
    {
        this->nkb += ucell.atoms[it].ncpp.nh * ucell.atoms[it].na;
    }

    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "TOTAL NUMBER OF NONLOCAL PROJECTORS", nkb);

    if (this->nhm > 0)
    {
        this->indv.create(ntype, this->nhm);
        this->nhtol.create(ntype, this->nhm);
        this->nhtolm.create(ntype, this->nhm);
        this->nhtoj.create(ntype, this->nhm);
        this->deeq.create(PARAM.inp.nspin, ucell.nat, this->nhm, this->nhm);
        this->deeq_nc.create(PARAM.inp.nspin, ucell.nat, this->nhm, this->nhm);
        this->qq_nt.create(ntype, this->nhm, this->nhm);
        this->qq_so.create(ntype, 4, this->nhm, this->nhm);
        if (this->use_gpu_)
        {
            if (PARAM.globalv.has_float_data)
            {
                resmem_sd_op()(s_deeq, PARAM.inp.nspin * ucell.nat * this->nhm * this->nhm);
                resmem_sd_op()(s_nhtol, ntype * this->nhm);
                resmem_sd_op()(s_nhtolm, ntype * this->nhm);
                resmem_sd_op()(s_indv, ntype * this->nhm);
                resmem_sd_op()(s_qq_nt, ntype * this->nhm * this->nhm);
                resmem_cd_op()(c_deeq_nc, PARAM.inp.nspin * ucell.nat * this->nhm * this->nhm);
                resmem_cd_op()(c_qq_so, ntype * 4 * this->nhm * this->nhm);
            }
            if (PARAM.globalv.has_double_data)
            {
                resmem_zd_op()(z_deeq_nc, PARAM.inp.nspin * ucell.nat * this->nhm * this->nhm);
                resmem_zd_op()(z_qq_so, ntype * 4 * this->nhm * this->nhm);
            }
            resmem_dd_op()(d_deeq, PARAM.inp.nspin * ucell.nat * this->nhm * this->nhm);
            resmem_dd_op()(d_indv, ntype * this->nhm);
            resmem_dd_op()(d_nhtol, ntype * this->nhm);
            resmem_dd_op()(d_nhtolm, ntype * this->nhm);
            resmem_dd_op()(d_qq_nt, ntype * this->nhm * this->nhm);
        }
        else
        {
            if (PARAM.globalv.has_float_data)
            {
                resmem_sh_op()(s_deeq,
                               PARAM.inp.nspin * ucell.nat * this->nhm * this->nhm,
                               "VNL::s_deeq");
                resmem_sh_op()(s_nhtol, ntype * this->nhm, "VNL::s_nhtol");
                resmem_sh_op()(s_nhtolm, ntype * this->nhm, "VNL::s_nhtolm");
                resmem_sh_op()(s_indv, ntype * this->nhm, "VNL::s_indv");
                resmem_sh_op()(s_qq_nt, ntype * this->nhm * this->nhm, "VNL::s_qq_nt");
                resmem_ch_op()(c_deeq_nc,
                               PARAM.inp.nspin * ucell.nat * this->nhm * this->nhm,
                               "VNL::c_deeq_nc");
                resmem_ch_op()(c_qq_so, ntype * 4 * this->nhm * this->nhm, "VNL::c_qq_so");
            }
            if (PARAM.globalv.has_double_data)
            {
                this->z_deeq_nc = this->deeq_nc.ptr;
                this->z_qq_so = this->qq_so.ptr;
            }
            this->d_deeq = this->deeq.ptr;
            this->d_indv = this->indv.c;
            this->d_nhtol = this->nhtol.c;
            this->d_nhtolm = this->nhtolm.c;
            this->d_qq_nt = this->qq_nt.ptr;
            // There's no need to delete double precision pointers while in a CPU environment.
        }
        this->dvan.create(ntype, this->nhm, this->nhm);
        this->dvan_so.create(PARAM.inp.nspin, ntype, this->nhm, this->nhm);

        this->ijtoh.create(ntype, this->nhm, this->nhm);
        this->qq_at.create(ucell.nat, this->nhm, this->nhm);
    }
    else
    {
        GlobalV::ofs_running << "\n nhm = 0, not allocate some matrix.";
    }

    // nqxq = ((sqrt(gcutm)+sqrt(xqq[1]*xqq[1]+xqq[2]*xqq[2]+xqq[3]*xqq[3])/
    // dq+4)*cell_factor;
    this->lmaxq = 2 * this->lmaxkb + 1;
    int npwx = this->wfcpw->npwk_max;
    this->vkbnc = npwx;
    if (nkb > 0 && allocate_vkb)
    {
        if (!this->use_gpu_)
    {
        vkb.create(nkb, npwx);
        ModuleBase::Memory::record("VNL::vkb", nkb * npwx * sizeof(std::complex<double>));
        }
        // GPU path: vkb ComplexMatrix is not allocated.
        // Column dimension is stored in vkbnc for gemm/gemv leading dimension.
        // Actual GPU buffers (c_vkb/z_vkb) are allocated below.
    }

    // this->nqx = 10000;        // calculted in allocate_nlpot.f90
    // PARAM.globalv.nqxq = static_cast<int>(((sqrt(INPUT.ecutrho) + qnorm) / PARAM.globalv.dq + 4.0) * cell_factor);

    // mohan update 2021-02-22
    // liuyu update 2023-09-28
    if (nbetam > 0)
    {
        const int nbrx_nc = 2 * nbetam;
        // nbetam: max number of beta functions
        if (PARAM.inp.nspin != 4)
        {
            this->tab.create(ntype, nbetam, PARAM.globalv.nqx);
            ModuleBase::Memory::record("VNL::tab", ntype * nbetam * PARAM.globalv.nqx * sizeof(double));
        }
        else
        {
            this->tab.create(ntype, nbrx_nc, PARAM.globalv.nqx);
            ModuleBase::Memory::record("VNL::tab", ntype * nbrx_nc * PARAM.globalv.nqx * sizeof(double));
        }

        if (lmaxq > 0)
        {
            this->qrad.create(ntype, lmaxq, nbetam * (nbetam + 1) / 2, PARAM.globalv.nqxq);
        }
    }

    // mohan update 2021-02-22
    // liuyu update 2023-09-28
    if (nwfcm > 0)
    {
        int nchix_nc = 2 * nwfcm;
        // nwfcm : max number of atomic wavefunctions per atom
        if (PARAM.inp.nspin != 4)
        {
            this->tab_at.create(ntype, nwfcm, PARAM.globalv.nqx);
            ModuleBase::Memory::record("VNL::tab_at", ntype * nwfcm * PARAM.globalv.nqx * sizeof(double));
        }
        else
        {
            this->tab_at.create(ntype, nchix_nc, PARAM.globalv.nqx);
            ModuleBase::Memory::record("VNL::tab_at", ntype * nchix_nc * PARAM.globalv.nqx * sizeof(double));
        }
    }
    if (this->use_gpu_)
    {
        if (PARAM.globalv.has_float_data)
        {
            resmem_sd_op()(s_tab, this->tab.getSize());
            resmem_cd_op()(c_vkb, nkb * npwx);
        }
        resmem_zd_op()(z_vkb, nkb * npwx);
        resmem_dd_op()(d_tab, this->tab.getSize());
    }
    else
    {
        if (PARAM.globalv.has_float_data)
        {
            resmem_sh_op()(s_tab, this->tab.getSize());
            resmem_ch_op()(c_vkb, nkb * npwx);
        }
#ifdef __DSP
        base_device::memory::resize_memory_op_mt<std::complex<double>, base_device::DEVICE_CPU>()
        (this->z_vkb, this->vkb.size, "VNL::z_vkb");
        // memcpy(this->z_vkb,this->vkb.c,this->vkb.size*16);
#else
        this->z_vkb = this->vkb.c;
#endif
        this->d_tab = this->tab.ptr;
        // There's no need to delete double precision pointers while in a CPU environment.
    }

    ModuleBase::timer::end("ppcell_vnl", "init");
    return;
}

void pseudopot_cell_vnl::print_vnl(std::ofstream& ofs)
{
    output::printr3_d(ofs, " tab : ", tab);
}

// ----------------------------------------------------------------------
// scale the non-local pseudopotential tables
void pseudopot_cell_vnl::rescale_vnl(const double& omega_in)
{
    const double ratio = this->omega_old / omega_in;
    const double sqrt_ratio = std::sqrt(ratio);
    this->omega_old = omega_in;

    for (int i = 0; i < this->tab.getSize(); i++)
    {
        this->tab.ptr[i] *= sqrt_ratio;
    }
    for (int i = 0; i < this->tab_at.getSize(); i++)
    {
        this->tab_at.ptr[i] *= sqrt_ratio;
    }
    for (int i = 0; i < this->qrad.getSize(); i++)
    {
        this->qrad.ptr[i] *= ratio;
    }
    if (this->use_gpu_)
    {
        if (this->s_tab != nullptr)
        {
            castmem_d2s_h2d_op()(this->s_tab, this->tab.ptr, this->tab.getSize());
        }
        // The double table is also used by GPU force and stress paths.
        syncmem_d2d_h2d_op()(this->d_tab, this->tab.ptr, this->tab.getSize());
    }
    else if (this->s_tab != nullptr)
    {
        castmem_d2s_h2h_op()(this->s_tab, this->tab.ptr, this->tab.getSize());
    }
}

template <>
float* pseudopot_cell_vnl::get_nhtol_data() const
{
    return this->s_nhtol;
}
template <>
double* pseudopot_cell_vnl::get_nhtol_data() const
{
    return this->d_nhtol;
}

template <>
float* pseudopot_cell_vnl::get_nhtolm_data() const
{
    return this->s_nhtolm;
}
template <>
double* pseudopot_cell_vnl::get_nhtolm_data() const
{
    return this->d_nhtolm;
}

template <>
float* pseudopot_cell_vnl::get_indv_data() const
{
    return this->s_indv;
}
template <>
double* pseudopot_cell_vnl::get_indv_data() const
{
    return this->d_indv;
}

template <>
float* pseudopot_cell_vnl::get_tab_data() const
{
    return this->s_tab;
}
template <>
double* pseudopot_cell_vnl::get_tab_data() const
{
    return this->d_tab;
}

template <>
float* pseudopot_cell_vnl::get_deeq_data() const
{
    return this->s_deeq;
}
template <>
double* pseudopot_cell_vnl::get_deeq_data() const
{
    return this->d_deeq;
}

template <>
float* pseudopot_cell_vnl::get_qq_nt_data() const
{
    return this->s_qq_nt;
}
template <>
double* pseudopot_cell_vnl::get_qq_nt_data() const
{
    return this->d_qq_nt;
}

template <>
std::complex<float>* pseudopot_cell_vnl::get_vkb_data() const
{
    return this->c_vkb;
}
template <>
std::complex<double>* pseudopot_cell_vnl::get_vkb_data() const
{
    return this->z_vkb;
}

template <>
std::complex<float>* pseudopot_cell_vnl::get_deeq_nc_data() const
{
    return this->c_deeq_nc;
}
template <>
std::complex<double>* pseudopot_cell_vnl::get_deeq_nc_data() const
{
    return this->z_deeq_nc;
}

template <>
std::complex<float>* pseudopot_cell_vnl::get_qq_so_data() const
{
    return this->c_qq_so;
}
template <>
std::complex<double>* pseudopot_cell_vnl::get_qq_so_data() const
{
    return this->z_qq_so;
}
