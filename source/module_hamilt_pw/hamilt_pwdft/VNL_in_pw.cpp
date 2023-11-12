#include "VNL_in_pw.h"

#include "module_base/clebsch_gordan_coeff.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/math_integral.h"
#include "module_base/math_polyint.h"
#include "module_base/math_sphbes.h"
#include "module_base/math_ylmreal.h"
#include "module_base/memory.h"
#include "module_base/timer.h"
#include "module_basis/module_ao/ORB_gen_tables.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_hamilt_pw/hamilt_pwdft/kernels/vnl_op.h"
#include "module_hamilt_pw/hamilt_pwdft/wavefunc.h"
#include "module_psi/kernels/device.h"

pseudopot_cell_vnl::pseudopot_cell_vnl()
{
}

pseudopot_cell_vnl::~pseudopot_cell_vnl()
{
    if(GlobalV::use_paw) return;
    delete[] indv_ijkb0;
    if (GlobalV::device_flag == "gpu")
    {
        if (GlobalV::precision_flag == "single")
        {
            delmem_sd_op()(gpu_ctx, this->s_deeq);
            delmem_sd_op()(gpu_ctx, this->s_nhtol);
            delmem_sd_op()(gpu_ctx, this->s_nhtolm);
            delmem_sd_op()(gpu_ctx, this->s_indv);
            delmem_sd_op()(gpu_ctx, this->s_tab);
            delmem_sd_op()(gpu_ctx, this->s_qq_nt);
            delmem_cd_op()(gpu_ctx, this->c_deeq_nc);
            delmem_cd_op()(gpu_ctx, this->c_vkb);
            delmem_cd_op()(gpu_ctx, this->c_qq_so);
        }
        else
        {
            delmem_zd_op()(gpu_ctx, this->z_deeq_nc);
            delmem_zd_op()(gpu_ctx, this->z_qq_so);
        }
        delmem_dd_op()(gpu_ctx, this->d_deeq);
        delmem_zd_op()(gpu_ctx, this->z_vkb);
        delmem_dd_op()(gpu_ctx, this->d_tab);
        delmem_dd_op()(gpu_ctx, this->d_indv);
        delmem_dd_op()(gpu_ctx, this->d_nhtol);
        delmem_dd_op()(gpu_ctx, this->d_nhtolm);
        delmem_dd_op()(gpu_ctx, this->d_qq_nt);
    }
    else
    {
        if (GlobalV::precision_flag == "single")
        {
            delmem_sh_op()(cpu_ctx, this->s_deeq);
            delmem_sh_op()(cpu_ctx, this->s_nhtol);
            delmem_sh_op()(cpu_ctx, this->s_nhtolm);
            delmem_sh_op()(cpu_ctx, this->s_indv);
            delmem_sh_op()(cpu_ctx, this->s_tab);
            delmem_sh_op()(cpu_ctx, this->s_qq_nt);
            delmem_ch_op()(cpu_ctx, this->c_deeq_nc);
            delmem_ch_op()(cpu_ctx, this->c_vkb);
            delmem_ch_op()(cpu_ctx, this->c_qq_so);
        }
        // There's no need to delete double precision pointers while in a CPU environment.
    }
}

//-----------------------------------
// setup lmaxkb, nhm, nkb, lmaxq
// allocate vkb, GlobalV::NQX, tab, tab_at
//-----------------------------------
void pseudopot_cell_vnl::init(const int ntype,
                              Structure_Factor* psf_in,
                              const ModulePW::PW_Basis_K* wfc_basis,
                              const bool allocate_vkb)
{
    if(GlobalV::use_paw) return;

    ModuleBase::TITLE("pseudopot_cell_vnl", "init");
    ModuleBase::timer::tick("ppcell_vnl", "init");

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
        GlobalV::ofs_running << " " << GlobalC::ucell.atoms[it].label << " non-local projectors:" << std::endl;
        for (int ibeta = 0; ibeta < GlobalC::ucell.atoms[it].ncpp.nbeta; ibeta++)
        {
            GlobalV::ofs_running << " projector " << ibeta + 1 << " L=" << GlobalC::ucell.atoms[it].ncpp.lll[ibeta]
                                 << std::endl;
            this->lmaxkb = std::max(this->lmaxkb, GlobalC::ucell.atoms[it].ncpp.lll[ibeta]);
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
        this->nhm = std::max(nhm, GlobalC::ucell.atoms[it].ncpp.nh);
        this->nbetam = std::max(nbetam, GlobalC::ucell.atoms[it].ncpp.nbeta);
        nwfcm = std::max(nwfcm, GlobalC::ucell.atoms[it].ncpp.nchi);
    }

    //----------------------------------------------------------
    // MEMBER VARIABLE :
    // NAME : nkb(total number of beta functions, with struct.fact.)
    //----------------------------------------------------------
    this->nkb = 0;
    for (it = 0; it < ntype; it++)
    {
        this->nkb += GlobalC::ucell.atoms[it].ncpp.nh * GlobalC::ucell.atoms[it].na;
    }

    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "TOTAL NUMBER OF NONLOCAL PROJECTORS", nkb);

    if (this->nhm > 0)
    {
        this->indv.create(ntype, this->nhm);
        this->nhtol.create(ntype, this->nhm);
        this->nhtolm.create(ntype, this->nhm);
        this->nhtoj.create(ntype, this->nhm);
        this->deeq.create(GlobalV::NSPIN, GlobalC::ucell.nat, this->nhm, this->nhm);
        this->deeq_nc.create(GlobalV::NSPIN, GlobalC::ucell.nat, this->nhm, this->nhm);
        this->qq_nt.create(ntype, this->nhm, this->nhm);
        this->qq_so.create(ntype, 4, this->nhm, this->nhm);
        if (GlobalV::device_flag == "gpu")
        {
            if (GlobalV::precision_flag == "single")
            {
                resmem_sd_op()(gpu_ctx, s_deeq, GlobalV::NSPIN * GlobalC::ucell.nat * this->nhm * this->nhm);
                resmem_sd_op()(gpu_ctx, s_nhtol, ntype * this->nhm);
                resmem_sd_op()(gpu_ctx, s_nhtolm, ntype * this->nhm);
                resmem_sd_op()(gpu_ctx, s_indv, ntype * this->nhm);
                resmem_sd_op()(gpu_ctx, s_qq_nt, ntype * this->nhm * this->nhm);
                resmem_cd_op()(gpu_ctx, c_deeq_nc, GlobalV::NSPIN * GlobalC::ucell.nat * this->nhm * this->nhm);
                resmem_cd_op()(gpu_ctx, c_qq_so, ntype * 4 * this->nhm * this->nhm);
            }
            else
            {
                resmem_zd_op()(gpu_ctx, z_deeq_nc, GlobalV::NSPIN * GlobalC::ucell.nat * this->nhm * this->nhm);
                resmem_zd_op()(gpu_ctx, z_qq_so, ntype * 4 * this->nhm * this->nhm);
            }
            resmem_dd_op()(gpu_ctx, d_deeq, GlobalV::NSPIN * GlobalC::ucell.nat * this->nhm * this->nhm);
            resmem_dd_op()(gpu_ctx, d_indv, ntype * this->nhm);
            resmem_dd_op()(gpu_ctx, d_nhtol, ntype * this->nhm);
            resmem_dd_op()(gpu_ctx, d_nhtolm, ntype * this->nhm);
            resmem_dd_op()(gpu_ctx, d_qq_nt, ntype * this->nhm * this->nhm);
        }
        else
        {
            if (GlobalV::precision_flag == "single")
            {
                resmem_sh_op()(cpu_ctx,
                               s_deeq,
                               GlobalV::NSPIN * GlobalC::ucell.nat * this->nhm * this->nhm,
                               "VNL::s_deeq");
                resmem_sh_op()(cpu_ctx, s_nhtol, ntype * this->nhm, "VNL::s_nhtol");
                resmem_sh_op()(cpu_ctx, s_nhtolm, ntype * this->nhm, "VNL::s_nhtolm");
                resmem_sh_op()(cpu_ctx, s_indv, ntype * this->nhm, "VNL::s_indv");
                resmem_sh_op()(cpu_ctx, s_qq_nt, ntype * this->nhm * this->nhm, "VNL::s_qq_nt");
                resmem_ch_op()(cpu_ctx,
                               c_deeq_nc,
                               GlobalV::NSPIN * GlobalC::ucell.nat * this->nhm * this->nhm,
                               "VNL::c_deeq_nc");
                resmem_ch_op()(cpu_ctx, c_qq_so, ntype * 4 * this->nhm * this->nhm, "VNL::c_qq_so");
            }
            else
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
        this->dvan_so.create(GlobalV::NSPIN, ntype, this->nhm, this->nhm);

        this->ijtoh.create(ntype, this->nhm, this->nhm);
        this->qq_at.create(GlobalC::ucell.nat, this->nhm, this->nhm);
    }
    else
    {
        GlobalV::ofs_running << "\n nhm = 0, not allocate some matrix.";
    }

    // nqxq = ((sqrt(gcutm)+sqrt(xqq[1]*xqq[1]+xqq[2]*xqq[2]+xqq[3]*xqq[3])/
    // dq+4)*cell_factor;
    this->lmaxq = 2 * this->lmaxkb + 1;
    int npwx = this->wfcpw->npwk_max;
    if (nkb > 0 && allocate_vkb)
    {
        vkb.create(nkb, npwx);
        ModuleBase::Memory::record("VNL::vkb", nkb * npwx * sizeof(double));
    }

    // this->nqx = 10000;		// calculted in allocate_nlpot.f90
    GlobalV::NQX = static_cast<int>((sqrt(INPUT.ecutwfc) / GlobalV::DQ + 4.0) * cell_factor);
    GlobalV::NQXQ = static_cast<int>((sqrt(INPUT.ecutrho) / GlobalV::DQ + 4.0) * cell_factor);
    // GlobalV::NQXQ = static_cast<int>(((sqrt(INPUT.ecutrho) + qnorm) / GlobalV::DQ + 4.0) * cell_factor);

    // mohan update 2021-02-22
    // liuyu update 2023-09-28
    if (nbetam > 0)
    {
        const int nbrx_nc = 2 * nbetam;
        // nbetam: max number of beta functions
        if (GlobalV::NSPIN != 4)
        {
            this->tab.create(ntype, nbetam, GlobalV::NQX);
            ModuleBase::Memory::record("VNL::tab", ntype * nbetam * GlobalV::NQX * sizeof(double));
        }
        else
        {
            this->tab.create(ntype, nbrx_nc, GlobalV::NQX);
            ModuleBase::Memory::record("VNL::tab", ntype * nbrx_nc * GlobalV::NQX * sizeof(double));
        }

        if (lmaxq > 0)
        {
            this->qrad.create(ntype, lmaxq, nbetam * (nbetam + 1) / 2, GlobalV::NQXQ);
        }
    }

    // mohan update 2021-02-22
    // liuyu update 2023-09-28
    if (nwfcm > 0)
    {
        int nchix_nc = 2 * nwfcm;
        // nwfcm : max number of atomic wavefunctions per atom
        if (GlobalV::NSPIN != 4)
        {
            this->tab_at.create(ntype, nwfcm, GlobalV::NQX);
            ModuleBase::Memory::record("VNL::tab_at", ntype * nwfcm * GlobalV::NQX * sizeof(double));
        }
        else
        {
            this->tab_at.create(ntype, nchix_nc, GlobalV::NQX);
            ModuleBase::Memory::record("VNL::tab_at", ntype * nchix_nc * GlobalV::NQX * sizeof(double));
        }
    }
    if (GlobalV::device_flag == "gpu")
    {
        if (GlobalV::precision_flag == "single")
        {
            resmem_sd_op()(gpu_ctx, s_tab, this->tab.getSize());
            resmem_cd_op()(gpu_ctx, c_vkb, nkb * npwx);
        }
        resmem_zd_op()(gpu_ctx, z_vkb, nkb * npwx);
        resmem_dd_op()(gpu_ctx, d_tab, this->tab.getSize());
    }
    else
    {
        if (GlobalV::precision_flag == "single")
        {
            resmem_sh_op()(cpu_ctx, s_tab, this->tab.getSize());
            resmem_ch_op()(cpu_ctx, c_vkb, nkb * npwx);
        }
        this->z_vkb = this->vkb.c;
        this->d_tab = this->tab.ptr;
        // There's no need to delete double precision pointers while in a CPU environment.
    }

    ModuleBase::timer::tick("ppcell_vnl", "init");
    return;
}

//----------------------------------------------------------
// Calculates beta functions (Kleinman-Bylander projectors),
// with structure factor, for all atoms, in reciprocal space
//----------------------------------------------------------
void pseudopot_cell_vnl::getvnl(const int& ik, ModuleBase::ComplexMatrix& vkb_in) const
{
    if(GlobalV::use_paw) return;
    if (GlobalV::test_pp)
        ModuleBase::TITLE("pseudopot_cell_vnl", "getvnl");
    ModuleBase::timer::tick("pp_cell_vnl", "getvnl");

    if (lmaxkb < 0)
    {
        return;
    }

    const int npw = this->wfcpw->npwk[ik];

    // When the internal memory is large enough, it is better to make vkb1 be the number of pseudopot_cell_vnl.
    // We only need to initialize it once as long as the cell is unchanged.
    ModuleBase::matrix vkb1(nhm, npw);
    double* vq = new double[npw];
    const int x1 = (lmaxkb + 1) * (lmaxkb + 1);

    ModuleBase::matrix ylm(x1, npw);
    ModuleBase::Memory::record("VNL::ylm", x1 * npw * sizeof(double));
    ModuleBase::Vector3<double>* gk = new ModuleBase::Vector3<double>[npw];
    for (int ig = 0; ig < npw; ig++)
    {
        gk[ig] = this->wfcpw->getgpluskcar(ik, ig);
    }

    ModuleBase::YlmReal::Ylm_Real(cpu_ctx, x1, npw, reinterpret_cast<double*>(gk), ylm.c);

    using Device = psi::DEVICE_CPU;
    Device* ctx = {};
    using resmem_complex_op = psi::memory::resize_memory_op<std::complex<double>, Device>;
    using delmem_complex_op = psi::memory::delete_memory_op<std::complex<double>, Device>;
    std::complex<double>* sk = nullptr;
    resmem_complex_op()(ctx, sk, GlobalC::ucell.nat * npw, "VNL::sk");
    this->psf->get_sk(ctx, ik, this->wfcpw, sk);

    int jkb = 0, iat = 0;
    for (int it = 0; it < GlobalC::ucell.ntype; it++)
    {
        // calculate beta in G-space using an interpolation table
        const int nbeta = GlobalC::ucell.atoms[it].ncpp.nbeta;
        const int nh = GlobalC::ucell.atoms[it].ncpp.nh;

        if (GlobalV::test_pp > 1)
            ModuleBase::GlobalFunc::OUT("nbeta", nbeta);

        for (int nb = 0; nb < nbeta; nb++)
        {
            if (GlobalV::test_pp > 1)
                ModuleBase::GlobalFunc::OUT("ib", nb);
            for (int ig = 0; ig < npw; ig++)
            {
                const double gnorm = gk[ig].norm() * GlobalC::ucell.tpiba;

                vq[ig] = ModuleBase::PolyInt::Polynomial_Interpolation(this->tab,
                                                                       it,
                                                                       nb,
                                                                       GlobalV::NQX,
                                                                       GlobalV::DQ,
                                                                       gnorm);
            }

            // add spherical harmonic part
            for (int ih = 0; ih < nh; ih++)
            {
                if (nb == this->indv(it, ih))
                {
                    const int lm = static_cast<int>(nhtolm(it, ih));
                    for (int ig = 0; ig < npw; ig++)
                    {
                        vkb1(ih, ig) = ylm(lm, ig) * vq[ig];
                    }
                }
            } // end ih
        }     // end nbeta

        // vkb1 contains all betas including angular part for type nt
        // now add the structure factor and factor (-i)^l
        for (int ia = 0; ia < GlobalC::ucell.atoms[it].na; ia++)
        {
            for (int ih = 0; ih < nh; ih++)
            {
                std::complex<double> pref = pow(ModuleBase::NEG_IMAG_UNIT, nhtol(it, ih)); //?
                std::complex<double>* pvkb = &vkb_in(jkb, 0);
                for (int ig = 0; ig < npw; ig++)
                {
                    pvkb[ig] = vkb1(ih, ig) * sk[iat * npw + ig] * pref;
                }
                ++jkb;
            } // end ih
            iat++;
        } // end ia
    }     // enddo

    delete[] gk;
    delete[] vq;
    delmem_complex_op()(ctx, sk);
    ModuleBase::timer::tick("pp_cell_vnl", "getvnl");

    return;
} // end subroutine getvnl

template <typename FPTYPE, typename Device>
void pseudopot_cell_vnl::getvnl(Device* ctx, const int& ik, std::complex<FPTYPE>* vkb_in) const
{
    if(GlobalV::use_paw) return;
    if (GlobalV::test_pp)
        ModuleBase::TITLE("pseudopot_cell_vnl", "getvnl");
    ModuleBase::timer::tick("pp_cell_vnl", "getvnl");

    using cal_vnl_op = hamilt::cal_vnl_op<FPTYPE, Device>;
    using resmem_int_op = psi::memory::resize_memory_op<int, Device>;
    using delmem_int_op = psi::memory::delete_memory_op<int, Device>;
    using syncmem_int_op = psi::memory::synchronize_memory_op<int, Device, psi::DEVICE_CPU>;
    using resmem_var_op = psi::memory::resize_memory_op<FPTYPE, Device>;
    using delmem_var_op = psi::memory::delete_memory_op<FPTYPE, Device>;
    using castmem_var_h2d_op = psi::memory::cast_memory_op<FPTYPE, double, Device, psi::DEVICE_CPU>;
    using castmem_var_h2h_op = psi::memory::cast_memory_op<FPTYPE, double, psi::DEVICE_CPU, psi::DEVICE_CPU>;
    using resmem_complex_op = psi::memory::resize_memory_op<std::complex<FPTYPE>, Device>;
    using delmem_complex_op = psi::memory::delete_memory_op<std::complex<FPTYPE>, Device>;

    if (lmaxkb < 0)
    {
        return;
    }

    const int x1 = (lmaxkb + 1) * (lmaxkb + 1);
    const int npw = this->wfcpw->npwk[ik];

    int *atom_nh = nullptr, *atom_na = nullptr, *atom_nb = nullptr, *h_atom_nh = new int[GlobalC::ucell.ntype],
        *h_atom_na = new int[GlobalC::ucell.ntype], *h_atom_nb = new int[GlobalC::ucell.ntype];
    for (int it = 0; it < GlobalC::ucell.ntype; it++)
    {
        h_atom_nb[it] = GlobalC::ucell.atoms[it].ncpp.nbeta;
        h_atom_nh[it] = GlobalC::ucell.atoms[it].ncpp.nh;
        h_atom_na[it] = GlobalC::ucell.atoms[it].na;
    }
    // When the internal memory is large enough, it is better to make vkb1 be the number of pseudopot_cell_vnl.
    // We only need to initialize it once as long as the cell is unchanged.
    FPTYPE *vkb1 = nullptr, *gk = nullptr, *ylm = nullptr, *_tab = this->get_tab_data<FPTYPE>(),
           *_indv = this->get_indv_data<FPTYPE>(), *_nhtol = this->get_nhtol_data<FPTYPE>(),
           *_nhtolm = this->get_nhtolm_data<FPTYPE>();
    resmem_var_op()(ctx, ylm, x1 * npw, "VNL::ylm");
    resmem_var_op()(ctx, vkb1, nhm * npw, "VNL::vkb1");

    ModuleBase::Vector3<double>* _gk = new ModuleBase::Vector3<double>[npw];
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096 / sizeof(FPTYPE))
#endif
    for (int ig = 0; ig < npw; ig++)
    {
        _gk[ig] = this->wfcpw->getgpluskcar(ik, ig);
    }
    if (GlobalV::device_flag == "gpu")
    {
        resmem_int_op()(ctx, atom_nh, GlobalC::ucell.ntype);
        resmem_int_op()(ctx, atom_nb, GlobalC::ucell.ntype);
        resmem_int_op()(ctx, atom_na, GlobalC::ucell.ntype);
        syncmem_int_op()(ctx, cpu_ctx, atom_nh, h_atom_nh, GlobalC::ucell.ntype);
        syncmem_int_op()(ctx, cpu_ctx, atom_nb, h_atom_nb, GlobalC::ucell.ntype);
        syncmem_int_op()(ctx, cpu_ctx, atom_na, h_atom_na, GlobalC::ucell.ntype);

        resmem_var_op()(ctx, gk, npw * 3);
        castmem_var_h2d_op()(ctx, cpu_ctx, gk, reinterpret_cast<double*>(_gk), npw * 3);
    }
    else
    {
        atom_nh = h_atom_nh;
        atom_nb = h_atom_nb;
        atom_na = h_atom_na;
        if (GlobalV::precision_flag == "single")
        {
            resmem_var_op()(ctx, gk, npw * 3);
            castmem_var_h2h_op()(cpu_ctx, cpu_ctx, gk, reinterpret_cast<double*>(_gk), npw * 3);
        }
        else
        {
            gk = reinterpret_cast<FPTYPE*>(_gk);
        }
    }

    ModuleBase::YlmReal::Ylm_Real(ctx, x1, npw, gk, ylm);

    std::complex<FPTYPE>* sk = nullptr;
    resmem_complex_op()(ctx, sk, GlobalC::ucell.nat * npw);
    this->psf->get_sk(ctx, ik, this->wfcpw, sk);

    cal_vnl_op()(ctx,
                 GlobalC::ucell.ntype,
                 npw,
                 this->wfcpw->npwk_max,
                 this->nhm,
                 GlobalV::NQX,
                 this->tab.getBound2(),
                 this->tab.getBound3(),
                 atom_na,
                 atom_nb,
                 atom_nh,
                 static_cast<FPTYPE>(GlobalV::DQ),
                 static_cast<FPTYPE>(GlobalC::ucell.tpiba),
                 static_cast<std::complex<FPTYPE>>(ModuleBase::NEG_IMAG_UNIT),
                 gk,
                 ylm,
                 _indv,
                 _nhtol,
                 _nhtolm,
                 _tab,
                 vkb1,
                 sk,
                 vkb_in);

    delete[] _gk;
    delete[] h_atom_nh;
    delete[] h_atom_na;
    delete[] h_atom_nb;
    delmem_var_op()(ctx, ylm);
    delmem_var_op()(ctx, vkb1);
    delmem_complex_op()(ctx, sk);
    if (psi::device::get_device_type<Device>(ctx) == psi::GpuDevice)
    {
        delmem_var_op()(ctx, gk);
        delmem_int_op()(ctx, atom_nh);
        delmem_int_op()(ctx, atom_nb);
        delmem_int_op()(ctx, atom_na);
    }
    ModuleBase::timer::tick("pp_cell_vnl", "getvnl");
} // end subroutine getvnl

void pseudopot_cell_vnl::init_vnl(UnitCell& cell, const ModulePW::PW_Basis* rho_basis)
{
    if(GlobalV::use_paw) return;
    ModuleBase::TITLE("pseudopot_cell_vnl", "init_vnl");
    ModuleBase::timer::tick("ppcell_vnl", "init_vnl");

    // from init_us_1
    //    a) For each non vanderbilt pseudopotential it computes the D and
    //       the betar in the same form of the Vanderbilt pseudopotential.
    //    b) It computes the indices indv which establish the correspondence
    //       nh <-> beta in the atom
    //    c) It computes the indices nhtol which establish the correspondence
    //       nh <-> angular momentum of the beta function
    //    d) It computes the indices nhtolm which establish the correspondence
    //       nh <-> combined (l,m) index for the beta function.

    // the following prevents an out-of-bound error: upf(nt)%nqlc=2*lmax+1
    // but in some versions of the PP files lmax is not set to the maximum
    // l of the beta functions but includes the l of the local potential
    for (int it = 0; it < cell.ntype; it++)
    {
        if (cell.atoms[it].ncpp.tvanp)
        {
            cell.atoms[it].ncpp.nqlc = std::min(cell.atoms[it].ncpp.nqlc, lmaxq);
            if (cell.atoms[it].ncpp.nqlc < 0)
                cell.atoms[it].ncpp.nqlc = 0;
        }
    }

    // In the spin-orbit case we need the unitary matrix u which rotates the
    // real spherical harmonics and yields the complex ones.
    soc.fcoef.create(cell.ntype, this->nhm, this->nhm);
    if (GlobalV::LSPINORB)
    {
        soc.rot_ylm(this->lmaxkb);
    }

    // For each pseudopotential we initialize the indices nhtol, nhtolm,
    // nhtoj, indv, and if the pseudopotential is of KB type we initialize
    // the atomic D terms
    this->dvan.zero_out();
    this->dvan_so.zero_out(); // added by zhengdy-soc
    delete[] indv_ijkb0;
    this->indv_ijkb0 = new int[GlobalC::ucell.nat];
    int ijkb0 = 0;
    for (int it = 0; it < cell.ntype; it++)
    {
        int BetaIndex = 0;
        const int Nprojectors = cell.atoms[it].ncpp.nh;
        for (int ib = 0; ib < cell.atoms[it].ncpp.nbeta; ib++)
        {
            const int l = cell.atoms[it].ncpp.lll[ib];
            const double j = cell.atoms[it].ncpp.jjj[ib];
            for (int m = 0; m < 2 * l + 1; m++)
            {
                this->nhtol(it, BetaIndex) = l;
                this->nhtolm(it, BetaIndex) = l * l + m;
                this->nhtoj(it, BetaIndex) = j;
                this->indv(it, BetaIndex) = ib;
                ++BetaIndex;
            }
        }

        // ijtoh map augmentation channel indexes ih and jh to composite
        // "triangular" index ijh
        for (int ih1 = 0; ih1 < nhm; ih1++)
        {
            for (int ih2 = ih1; ih2 < nhm; ih2++)
            {
                this->ijtoh(it, ih1, ih2) = -1;
                this->ijtoh(it, ih2, ih1) = -1;
            }
        }
        int ijv = 0;
        for (int ih1 = 0; ih1 < Nprojectors; ih1++)
        {
            for (int ih2 = ih1; ih2 < Nprojectors; ih2++)
            {
                ijv++; // liuyu I'm not sure from 0 or 1
                this->ijtoh(it, ih1, ih2) = ijv;
                this->ijtoh(it, ih2, ih1) = ijv;
            }
        }

        // ijkb0 points to the first beta "in the solid" for atom ia
        // i.e. ijkb0,.. ijkb0+nh(ityp(ia))-1 are the nh beta functions of
        // atom ia in the global list of beta functions (ijkb0=0 for ia=0)
        for (int ia = 0; ia < cell.nat; ia++)
        {
            if (it == cell.iat2it[ia])
            {
                this->indv_ijkb0[ia] = ijkb0;
                ijkb0 += Nprojectors;
            }
        }

        // From now on the only difference between KB and US pseudopotentials
        // is in the presence of the q and Q functions.
        // Here we initialize the D of the solid
        if (cell.atoms[it].ncpp.has_so)
        {
            for (int ip = 0; ip < Nprojectors; ip++)
            {
                const int l1 = this->nhtol(it, ip);
                const double j1 = this->nhtoj(it, ip);
                const int m1 = this->nhtolm(it, ip) - l1 * l1;
                // const int v1 = static_cast<int>( indv(it, ip ) );
                for (int ip2 = 0; ip2 < Nprojectors; ip2++)
                {
                    const int l2 = this->nhtol(it, ip2);
                    const double j2 = this->nhtoj(it, ip2);
                    const int m2 = this->nhtolm(it, ip2) - l2 * l2;
                    // const int v2 = static_cast<int>( indv(it, ip2 ) );
                    if (l1 == l2 && fabs(j1 - j2) < 1e-7)
                    {
                        for (int is1 = 0; is1 < 2; is1++)
                        {
                            for (int is2 = 0; is2 < 2; is2++)
                            {
                                soc.set_fcoef(l1, l2, is1, is2, m1, m2, j1, j2, it, ip, ip2);
                            }
                        }
                    }
                }
            }
            //
            //   and calculate the bare coefficients
            //
            for (int ip = 0; ip < Nprojectors; ++ip)
            {
                const int ir = static_cast<int>(indv(it, ip));
                for(int ip2=0; ip2<Nprojectors; ++ip2)
				{
					const int is = static_cast<int>( indv(it, ip2) );
					int ijs =0;
					for(int is1=0;is1<2;++is1)
					{
						for(int is2=0;is2<2;++is2)
						{
							this->dvan_so(ijs,it,ip,ip2) = cell.atoms[it].ncpp.dion(ir, is) * soc.fcoef(it,is1,is2,ip,ip2);
							++ijs;
							if(ir != is) soc.fcoef(it,is1,is2,ip,ip2) = std::complex<double>(0.0,0.0);
						}
					}
				}
            }
        }
        else
            for (int ip = 0; ip < Nprojectors; ip++)
            {
                for (int ip2 = 0; ip2 < Nprojectors; ip2++)
                {
                    if (this->nhtol(it, ip) == nhtol(it, ip2) && this->nhtolm(it, ip) == nhtolm(it, ip2))
                    {
                        const int ir = static_cast<int>(indv(it, ip));
                        const int is = static_cast<int>(indv(it, ip2));
                        if (GlobalV::LSPINORB)
                        {
                            this->dvan_so(0, it, ip, ip2) = cell.atoms[it].ncpp.dion(ir, is);
                            this->dvan_so(3, it, ip, ip2) = cell.atoms[it].ncpp.dion(ir, is);
                        }
                        else
                        {
                            this->dvan(it, ip, ip2) = cell.atoms[it].ncpp.dion(ir, is);
                        }
                    }
                }
            }
    }

    // e) It computes the coefficients c_{LM}^{nm} which relates the
    //    spherical harmonics in the Q expansion
    // f) It computes the interpolation table "qrad" for Q(G)
    // g) It computes the qq terms which define the S matrix.

    // compute Clebsch-Gordan coefficients
    if (GlobalV::use_uspp)
    {
        ModuleBase::Clebsch_Gordan::clebsch_gordan(lmaxkb + 1, ap, lpx, lpl);
    }

    // here for the US types we compute the Fourier transform of the Q functions.
    if (lmaxq > 0)
    {
        this->compute_qrad(cell);
    }

    // compute the qq coefficients by integrating the Q.
    // The qq are the g=0 components of Q
    if (rho_basis->ig_gge0 >= 0)
    {
        ModuleBase::matrix ylmk0(lmaxq * lmaxq, 1);
        const double qnorm = 0.0; // only G=0 term
        std::complex<double> qgm(0.0, 0.0);
        ModuleBase::YlmReal::Ylm_Real(lmaxq * lmaxq, 1, &(rho_basis->gcar[rho_basis->ig_gge0]), ylmk0);
        for (int it = 0; it < cell.ntype; it++)
        {
            Atom_pseudo* upf = &cell.atoms[it].ncpp;
            if (upf->tvanp)
            {
                if (upf->has_so)
                {
                    for (int ih = 0; ih < upf->nh; ih++)
                    {
                        for (int jh = 0; jh < upf->nh; jh++)
                        {
                            this->radial_fft_q(1, ih, jh, it, &qnorm, ylmk0, &qgm);
                            this->qq_nt(it, ih, jh) = cell.omega * qgm.real();
                            for (int kh = 0; kh < upf->nh; kh++)
                            {
                                for (int lh = 0; lh < upf->nh; lh++)
                                {
                                    int ijs = 0;
                                    for (int is1 = 0; is1 < 2; ++is1)
                                    {
                                        for (int is2 = 0; is2 < 2; ++is2)
                                        {
                                            for (int is = 0; is < 2; is++)
                                            {
                                                this->qq_so(it, ijs, kh, lh) += cell.omega * qgm.real()
                                                                                * soc.fcoef(it, is1, is, kh, ih)
                                                                                * soc.fcoef(it, is, is2, jh, lh);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                else
                {
                    for (int ih = 0; ih < upf->nh; ih++)
                    {
                        for (int jh = ih; jh < upf->nh; jh++)
                        {
                            this->radial_fft_q(1, ih, jh, it, &qnorm, ylmk0, &qgm);
                            if (GlobalV::LSPINORB)
                            {
                                this->qq_so(it, 0, ih, jh) = cell.omega * qgm.real();
                                this->qq_so(it, 0, jh, ih) = this->qq_so(it, 0, ih, jh);
                                this->qq_so(it, 3, ih, jh) = this->qq_so(it, 0, ih, jh);
                                this->qq_so(it, 3, jh, ih) = this->qq_so(it, 0, ih, jh);
                            }
                            this->qq_nt(it, ih, jh) = cell.omega * qgm.real();
                            this->qq_nt(it, jh, ih) = cell.omega * qgm.real();
                        }
                    }
                }
            }
        }
    }

#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, this->qq_nt.ptr, this->qq_nt.getSize(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, this->qq_so.ptr, this->qq_so.getSize(), MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
#endif

    // set the atomic specific qq_at matrices
    for (int ia = 0; ia < cell.nat; ia++)
    {
        int it = cell.iat2it[ia];
        for (int ih = 0; ih < nhm; ih++)
        {
            for (int jh = ih; jh < nhm; jh++)
            {
                this->qq_at(ia, ih, jh) = qq_nt(it, ih, jh);
                this->qq_at(ia, jh, ih) = qq_nt(it, jh, ih);
            }
        }
    }

    // h) It fills the interpolation table for the beta functions
    /**********************************************************
    // He Lixin: this block is used for non-local potential
    // fill the interpolation table tab
    ************************************************************/

    const double pref = ModuleBase::FOUR_PI / sqrt(cell.omega);
    this->tab.zero_out();
    GlobalV::ofs_running << "\n Init Non-Local PseudoPotential table : ";
    for (int it = 0; it < cell.ntype; it++)
    {
        const int nbeta = cell.atoms[it].ncpp.nbeta;
        int kkbeta = cell.atoms[it].ncpp.kkbeta;

        // mohan modify 2008-3-31
        // mohan add kkbeta>0 2009-2-27
        if ((kkbeta % 2 == 0) && kkbeta > 0)
        {
            kkbeta--;
        }

        double* jl = new double[kkbeta];
        double* aux = new double[kkbeta];

        for (int ib = 0; ib < nbeta; ib++)
        {
            const int l = cell.atoms[it].ncpp.lll[ib];
            for (int iq = 0; iq < GlobalV::NQX; iq++)
            {
                const double q = iq * GlobalV::DQ;
                ModuleBase::Sphbes::Spherical_Bessel(kkbeta, cell.atoms[it].ncpp.r, q, l, jl);

                for (int ir = 0; ir < kkbeta; ir++)
                {
                    aux[ir] = cell.atoms[it].ncpp.betar(ib, ir) * jl[ir] * cell.atoms[it].ncpp.r[ir];
                }
                double vqint;
                ModuleBase::Integral::Simpson_Integral(kkbeta, aux, cell.atoms[it].ncpp.rab, vqint);
                this->tab(it, ib, iq) = vqint * pref;
            }
        }
        delete[] aux;
        delete[] jl;
    }
    if (GlobalV::device_flag == "gpu") {
        if (GlobalV::precision_flag == "single") {
            castmem_d2s_h2d_op()(gpu_ctx, cpu_ctx, this->s_indv, this->indv.c, this->indv.nr * this->indv.nc);
            castmem_d2s_h2d_op()(gpu_ctx, cpu_ctx, this->s_nhtol, this->nhtol.c, this->nhtol.nr * this->nhtol.nc);
            castmem_d2s_h2d_op()(gpu_ctx, cpu_ctx, this->s_nhtolm, this->nhtolm.c, this->nhtolm.nr * this->nhtolm.nc);
            castmem_d2s_h2d_op()(gpu_ctx, cpu_ctx, this->s_tab, this->tab.ptr, this->tab.getSize());
            castmem_d2s_h2d_op()(gpu_ctx, cpu_ctx, this->s_qq_nt, this->qq_nt.ptr, this->qq_nt.getSize());
            castmem_z2c_h2d_op()(gpu_ctx, cpu_ctx, this->c_qq_so, this->qq_so.ptr, this->qq_so.getSize());
        }
        else
        {
            syncmem_z2z_h2d_op()(gpu_ctx, cpu_ctx, this->z_qq_so, this->qq_so.ptr, this->qq_so.getSize());
        }
        // Even when the single precision flag is enabled,
        // these variables are utilized in the Force/Stress calculation as well.
        // modified by denghuilu at 2023-05-15
        syncmem_d2d_h2d_op()(gpu_ctx, cpu_ctx, this->d_indv, this->indv.c, this->indv.nr * this->indv.nc);
        syncmem_d2d_h2d_op()(gpu_ctx, cpu_ctx, this->d_nhtol, this->nhtol.c, this->nhtol.nr * this->nhtol.nc);
        syncmem_d2d_h2d_op()(gpu_ctx, cpu_ctx, this->d_nhtolm, this->nhtolm.c, this->nhtolm.nr * this->nhtolm.nc);
        syncmem_d2d_h2d_op()(gpu_ctx, cpu_ctx, this->d_tab, this->tab.ptr, this->tab.getSize());
        syncmem_d2d_h2d_op()(gpu_ctx, cpu_ctx, this->d_qq_nt, this->qq_nt.ptr, this->qq_nt.getSize());
    }
    else {
        if (GlobalV::precision_flag == "single") {
            castmem_d2s_h2h_op()(cpu_ctx, cpu_ctx, this->s_indv, this->indv.c, this->indv.nr * this->indv.nc);
            castmem_d2s_h2h_op()(cpu_ctx, cpu_ctx, this->s_nhtol, this->nhtol.c, this->nhtol.nr * this->nhtol.nc);
            castmem_d2s_h2h_op()(cpu_ctx, cpu_ctx, this->s_nhtolm, this->nhtolm.c, this->nhtolm.nr * this->nhtolm.nc);
            castmem_d2s_h2h_op()(cpu_ctx, cpu_ctx, this->s_tab, this->tab.ptr, this->tab.getSize());
            castmem_d2s_h2h_op()(cpu_ctx, cpu_ctx, this->s_qq_nt, this->qq_nt.ptr, this->qq_nt.getSize());
            castmem_z2c_h2h_op()(cpu_ctx, cpu_ctx, this->c_qq_so, this->qq_so.ptr, this->qq_so.getSize());
        }
        // There's no need to synchronize double precision pointers while in a CPU environment.
    }
	ModuleBase::timer::tick("ppcell_vnl","init_vnl");
	GlobalV::ofs_running << "\n Init Non-Local-Pseudopotential done." << std::endl;
	return;
}

void pseudopot_cell_vnl::compute_qrad(UnitCell& cell)
{
    const double pref = ModuleBase::FOUR_PI / cell.omega;

    for (int it = 0; it < cell.ntype; it++)
    {
        Atom_pseudo* upf = &cell.atoms[it].ncpp;
        if (upf->tvanp)
        {
            const int nbeta = upf->nbeta;
            int kkbeta = upf->kkbeta;
            if ((kkbeta % 2 == 0) && kkbeta > 0)
            {
                kkbeta--;
            }
            double* aux = new double[kkbeta];
            double* besr = new double[kkbeta];

            for (int l = 0; l < upf->nqlc; l++)
            {
                for (int iq = 0; iq < GlobalV::NQXQ; iq++)
                {
                    const double q = iq * GlobalV::DQ;
                    // here we compute the spherical bessel function for each q_i
                    ModuleBase::Sphbes::Spherical_Bessel(kkbeta, upf->r, q, l, besr);
                    for (int nb = 0; nb < nbeta; nb++)
                    {
                        // the Q are symmetric with respect to indices nb and mb
                        for (int mb = nb; mb < nbeta; mb++)
                        {
                            const int ijv = mb * (mb + 1) / 2 + nb;
                            if ((l >= std::abs(upf->lll[nb] - upf->lll[mb])) && (l <= (upf->lll[nb] + upf->lll[mb]))
                                && ((l + upf->lll[nb] + upf->lll[mb]) % 2 == 0))
                            {
                                for (int ir = 0; ir < kkbeta; ir++)
                                {
                                    aux[ir] = besr[ir] * upf->qfuncl(l, ijv, ir);
                                }
                                // then we integrate with all the Q functions
                                double vqint;
                                ModuleBase::Integral::Simpson_Integral(kkbeta, aux, upf->rab, vqint);
                                qrad(it, l, ijv, iq) = vqint * pref;
                            }
                        }
                    }
                }
            }
            delete[] aux;
            delete[] besr;
        }
    }
}

void pseudopot_cell_vnl::radial_fft_q(const int ng,
                                      const int ih,
                                      const int jh,
                                      const int itype,
                                      const double* qnorm,
                                      const ModuleBase::matrix ylm,
                                      std::complex<double>* qg) const
{
    // computes the indices which correspond to ih,jh
    const int nb = indv(itype, ih);
    const int mb = indv(itype, jh);
    assert(nb < nbetam);
    assert(mb < nbetam);
    int ijv = 0;
    if (nb >= mb)
    {
        ijv = nb * (nb + 1) / 2 + mb;
    }
    else
    {
        ijv = mb * (mb + 1) / 2 + nb;
    }
    const int ivl = nhtolm(itype, ih);
    const int jvl = nhtolm(itype, jh);

    for (int ig = 0; ig < ng; ig++)
    {
        qg[ig] = {0, 0};
    }
    // makes the sum over the non zero LM
    int l = -1;
    std::complex<double> pref(0.0, 0.0);
    for (int lm = 0; lm < this->lpx(ivl, jvl); lm++)
    {
        int lp = this->lpl(ivl, jvl, lm);
        assert(lp >= 0);
        assert(lp < 49);
        if (lp == 0)
        {
            l = 0;
        }
        else if (lp < 4)
        {
            l = 1;
        }
        else if (lp < 9)
        {
            l = 2;
        }
        else if (lp < 16)
        {
            l = 3;
        }
        else if (lp < 25)
        {
            l = 4;
        }
        else if (lp < 36)
        {
            l = 5;
        }
        else
        {
            l = 6;
        }
        pref = pow(ModuleBase::NEG_IMAG_UNIT, l) * this->ap(lp, ivl, jvl);

        double qm1 = -1.0; // any number smaller than qnorm
        double work = 0.0;
        for (int ig = 0; ig < ng; ig++)
        {
            // calculate quantites depending on the module of G only when needed
            if (std::abs(qnorm[ig] - qm1) > 1e-6)
            {
                work = ModuleBase::PolyInt::Polynomial_Interpolation(this->qrad,
                                                                     itype,
                                                                     l,
                                                                     ijv,
                                                                     GlobalV::NQXQ,
                                                                     GlobalV::DQ,
                                                                     qnorm[ig]);
                qm1 = qnorm[ig];
            }
            qg[ig] += pref * work * ylm(lp, ig);
        }
    }
}

template <typename FPTYPE, typename Device>
void pseudopot_cell_vnl::radial_fft_q(Device* ctx,
                                      const int ng,
                                      const int ih,
                                      const int jh,
                                      const int itype,
                                      const FPTYPE* qnorm,
                                      const FPTYPE* ylm,
                                      std::complex<FPTYPE>* qg) const
{
    using setmem_complex_op = psi::memory::set_memory_op<std::complex<FPTYPE>, Device>;

    // computes the indices which correspond to ih,jh
    const int nb = indv(itype, ih);
    const int mb = indv(itype, jh);
    assert(nb < nbetam);
    assert(mb < nbetam);
    int ijv = 0;
    if (nb >= mb)
    {
        ijv = nb * (nb + 1) / 2 + mb;
    }
    else
    {
        ijv = mb * (mb + 1) / 2 + nb;
    }
    const int ivl = nhtolm(itype, ih);
    const int jvl = nhtolm(itype, jh);

    setmem_complex_op()(ctx, qg, 0, ng);

    const double* qnorm_double = reinterpret_cast<const double*>(qnorm);

    // makes the sum over the non zero LM
    int l = -1;
    std::complex<FPTYPE> pref(0.0, 0.0);
    for (int lm = 0; lm < this->lpx(ivl, jvl); lm++)
    {
        int lp = this->lpl(ivl, jvl, lm);
        assert(lp >= 0);
        assert(lp < 49);
        if (lp == 0)
        {
            l = 0;
        }
        else if (lp < 4)
        {
            l = 1;
        }
        else if (lp < 9)
        {
            l = 2;
        }
        else if (lp < 16)
        {
            l = 3;
        }
        else if (lp < 25)
        {
            l = 4;
        }
        else if (lp < 36)
        {
            l = 5;
        }
        else
        {
            l = 6;
        }
        pref = static_cast<std::complex<FPTYPE>>(pow(ModuleBase::NEG_IMAG_UNIT, l) * this->ap(lp, ivl, jvl));

        double qm1 = -1.0; // any number smaller than qnorm
        double work = 0.0;
        for (int ig = 0; ig < ng; ig++)
        {
            if (std::abs(qnorm_double[ig] - qm1) > 1e-6)
            {
                work = ModuleBase::PolyInt::Polynomial_Interpolation(this->qrad,
                                                                     itype,
                                                                     l,
                                                                     ijv,
                                                                     GlobalV::NQXQ,
                                                                     GlobalV::DQ,
                                                                     qnorm_double[ig]);
                qm1 = qnorm_double[ig];
            }
            qg[ig] += pref * static_cast<FPTYPE>(work) * ylm[lp * ng + ig];
        }
    }
}

#ifdef __LCAO
std::complex<double> pseudopot_cell_vnl::Cal_C(int alpha, int lu, int mu, int L, int M) // pengfei Li  2018-3-23
{
	std::complex<double> cf;
	if(alpha == 0)
	{
		cf = -sqrt(4*ModuleBase::PI/3)*CG(lu,mu,1,1,L,M);
	}
	else if(alpha == 1)
	{
		cf = -sqrt(4*ModuleBase::PI/3)*CG(lu,mu,1,2,L,M);
	}
	else if(alpha == 2)
	{
		cf = sqrt(4*ModuleBase::PI/3)*CG(lu,mu,1,0,L,M);
	}
	else
	{
		ModuleBase::WARNING_QUIT("pseudopot_cell_vnl_alpha", "alpha must be 0~2");
	}
	
	return cf;
}

double pseudopot_cell_vnl::CG(int l1, int m1, int l2, int m2, int L, int M)      // pengfei Li 2018-3-23
{
	int dim = L*L+M;
	int dim1 = l1*l1+m1;
	int dim2 = l2*l2+m2;
	
	//double A = MGT.Gaunt_Coefficients(dim1, dim2, dim);
	
	return MGT.Gaunt_Coefficients(dim1, dim2, dim);
}

// void pseudopot_cell_vnl::getvnl_alpha(const int &ik)           // pengfei Li  2018-3-23
// {
// 	if(GlobalV::test_pp) ModuleBase::TITLE("pseudopot_cell_vnl","getvnl_alpha");
// 	ModuleBase::timer::tick("pp_cell_vnl","getvnl_alpha");

// 	if(lmaxkb < 0) 
// 	{
// 		return;
// 	}
	
// 	const int npw = this->wfcpw->npwk[ik];
// 	int ig, ia, nb, ih, lu, mu;

// 	double *vq = new double[npw];
// 	const int x1= (lmaxkb + 2)*(lmaxkb + 2);

// 	ModuleBase::matrix ylm(x1, npw);
// 	ModuleBase::Vector3<double> *gk = new ModuleBase::Vector3<double>[npw];
// 	for (ig = 0;ig < npw;ig++) 
// 	{
// 		gk[ig] = this->wfcpw->getgpluskcar(ik,ig);
// 	}

// 	vkb1_alpha = new std::complex<double>**[3];
// 	for(int i=0; i<3; i++)
// 	{
// 		vkb1_alpha[i] = new std::complex<double>*[nhm];
// 		for(int j=0; j<nhm; j++)
// 		{
// 			vkb1_alpha[i][j] = new std::complex<double>[npw];
// 		}
// 	}	
	
// 	vkb_alpha = new std::complex<double>**[3];
// 	for(int i=0; i<3; i++)
// 	{
// 		vkb_alpha[i] = new std::complex<double>*[nkb];
// 		for(int j=0; j<nkb; j++)
// 		{
// 			vkb_alpha[i][j] = new std::complex<double>[this->wfcpw->npwk_max];
// 		}
// 	}
	
// 	ModuleBase::YlmReal::Ylm_Real(x1, npw, gk, ylm);

// 	MGT.init_Gaunt_CH( lmaxkb + 2 );
// 	MGT.init_Gaunt( lmaxkb + 2 );
	
// 	int jkb = 0;
// 	for(int it = 0;it < GlobalC::ucell.ntype;it++)
// 	{
// 		if(GlobalV::test_pp>1) ModuleBase::GlobalFunc::OUT("it",it);
// 		// calculate beta in G-space using an interpolation table
// 		const int nbeta = GlobalC::ucell.atoms[it].ncpp.nbeta;
// 		const int nh = GlobalC::ucell.atoms[it].ncpp.nh;

// 		if(GlobalV::test_pp>1) ModuleBase::GlobalFunc::OUT("nbeta",nbeta);

// 		for(int i=0; i<3; i++)
// 			for(int j=0; j<nhm; j++)
// 			{
// 				ModuleBase::GlobalFunc::ZEROS(vkb1_alpha[i][j], npw);
// 			}
			
// 		for (ih = 0;ih < nh; ih++)
// 		{
// 			lu = static_cast<int>( nhtol(it, ih));
// 			mu = static_cast<int>( nhtolm(it, ih)) - lu * lu;
// 			nb = static_cast<int>( indv(it, ih));
			
// 			for (int L= std::abs(lu - 1); L<= (lu + 1); L++)
// 			{
// 				for (ig = 0;ig < npw;ig++)
// 				{
// 					const double gnorm = gk[ig].norm() * GlobalC::ucell.tpiba;
// 					vq [ig] = ModuleBase::PolyInt::Polynomial_Interpolation(
// 							this->tab_alpha, it, nb, L, GlobalV::NQX, GlobalV::DQ, gnorm);
					
// 					for (int M=0; M<2*L+1; M++)
// 					{
// 						int lm = L*L + M;
// 						for (int alpha=0; alpha<3; alpha++)
// 						{
// 							std::complex<double> c = Cal_C(alpha,lu, mu, L, M);
// 							/*if(alpha == 0)
// 							{
// 								std::cout<<"lu= "<<lu<<"  mu= "<<mu<<"  L= "<<L<<"  M= "<<M<<" alpha = "<<alpha<<"  "<<c<<std::endl;
// 							}*/
// 							vkb1_alpha[alpha][ih][ig] += c * vq[ig] * ylm(lm, ig) * pow( ModuleBase::NEG_IMAG_UNIT, L);
// 						}	
// 					}
// 				}
// 			}
// 		} // end nbeta

// 		for (ia=0; ia<GlobalC::ucell.atoms[it].na; ia++)
// 		{
// 			std::complex<double> *sk = this->psf->get_sk(ik, it, ia,this->wfcpw);
// 			for (ih = 0;ih < nh;ih++)
// 			{
// 				for (ig = 0;ig < npw;ig++)
// 				{
// 					for(int alpha=0; alpha<3; alpha++)
// 					{
// 						vkb_alpha[alpha][jkb][ig] = vkb1_alpha[alpha][ih][ig] * sk [ig];
// 					}
// 				}
// 				++jkb;
// 			} // end ih
// 			delete [] sk;
// 		} // end ia
// 	} // enddo

// 	delete [] gk;
// 	delete [] vq;
// 	ModuleBase::timer::tick("pp_cell_vnl","getvnl_alpha");
// 	return;
// } 
#endif

void pseudopot_cell_vnl::init_vnl_alpha(void)          // pengfei Li 2018-3-23
{
	if(GlobalV::test_pp) ModuleBase::TITLE("pseudopot_cell_vnl","init_vnl_alpha");
	ModuleBase::timer::tick("ppcell_vnl","init_vnl_alpha");

	for(int it=0;it<GlobalC::ucell.ntype;it++)
	{
		int BetaIndex=0;
		//const int Nprojectors = GlobalC::ucell.atoms[it].nh;
		for (int ib=0; ib<GlobalC::ucell.atoms[it].ncpp.nbeta; ib++)
		{
			const int l = GlobalC::ucell.atoms[it].ncpp.lll [ib];
			for(int m=0; m<2*l+1; m++)
			{
				this->nhtol(it,BetaIndex) = l;
				this->nhtolm(it,BetaIndex) = l*l + m;
				this->indv(it,BetaIndex) = ib;
				++BetaIndex;
			}
		}
	} 


	// max number of beta functions
	const int nbrx = 10;

	const double pref = ModuleBase::FOUR_PI / sqrt(GlobalC::ucell.omega);
	this->tab_alpha.create(GlobalC::ucell.ntype, nbrx, lmaxkb+2, GlobalV::NQX);
	this->tab_alpha.zero_out();
	GlobalV::ofs_running<<"\n Init Non-Local PseudoPotential table( including L index) : ";
	for (int it = 0;it < GlobalC::ucell.ntype;it++)  
	{
		const int nbeta = GlobalC::ucell.atoms[it].ncpp.nbeta;
		int kkbeta = GlobalC::ucell.atoms[it].ncpp.kkbeta;

		//mohan modify 2008-3-31
		//mohan add kkbeta>0 2009-2-27
		if ( (kkbeta%2 == 0) && kkbeta>0 )
		{
			kkbeta--;
		}

		double *jl = new double[kkbeta];
		double *aux  = new double[kkbeta];

		for (int ib = 0;ib < nbeta;ib++)
		{
			for (int L = 0; L <= lmaxkb+1; L++)
			{
				for (int iq = 0; iq < GlobalV::NQX; iq++)
				{
					const double q = iq * GlobalV::DQ;
					ModuleBase::Sphbes::Spherical_Bessel(kkbeta, GlobalC::ucell.atoms[it].ncpp.r, q, L, jl);
					
					for (int ir = 0;ir < kkbeta;ir++)
					{
						aux[ir] = GlobalC::ucell.atoms[it].ncpp.betar(ib, ir) * jl[ir] * 
								  GlobalC::ucell.atoms[it].ncpp.r[ir] * GlobalC::ucell.atoms[it].ncpp.r[ir];
					}
					double vqint;
					ModuleBase::Integral::Simpson_Integral(kkbeta, aux, GlobalC::ucell.atoms[it].ncpp.rab, vqint);
					this->tab_alpha(it, ib, L, iq) = vqint * pref;
				}
			}
		} 
		delete[] aux;
		delete[] jl;
	}
	ModuleBase::timer::tick("ppcell_vnl","init_vnl_alpha");
	GlobalV::ofs_running << "\n Init Non-Local-Pseudopotential done(including L)." << std::endl;
	return;
}



void pseudopot_cell_vnl::print_vnl(std::ofstream &ofs)
{
	output::printr3_d(ofs, " tab : ", tab);
}

// ----------------------------------------------------------------------
void pseudopot_cell_vnl::cal_effective_D(const ModuleBase::matrix& veff,
                                         const ModulePW::PW_Basis* rho_basis,
                                         UnitCell& cell)
{
    if(GlobalV::use_paw) return;
    ModuleBase::TITLE("pseudopot_cell_vnl", "cal_effective_D");

    /*
	recalculate effective coefficient matrix for non-local pseudo-potential
	1. assign to each atom from element;
	2. extend to each spin when nspin larger than 1
	3. rotate to effective matrix when spin-orbital coupling is used
	*/

    if (!GlobalV::use_uspp)
    {
        for (int iat = 0; iat < cell.nat; iat++)
        {
            const int it = cell.iat2it[iat];
            const int nht = cell.atoms[it].ncpp.nh;
            // nht: number of beta functions per atom type
            for (int is = 0; is < GlobalV::NSPIN; is++)
            {
                for (int ih = 0; ih < nht; ih++)
                {
                    for (int jh = ih; jh < nht; jh++)
                    {
                        if (GlobalV::LSPINORB)
                        {
                            this->deeq_nc(is, iat, ih, jh) = this->dvan_so(is, it, ih, jh);
                            this->deeq_nc(is, iat, jh, ih) = this->dvan_so(is, it, jh, ih);
                        }
                        else if (GlobalV::NSPIN == 4)
                        {
                            if (is == 0)
                            {
                                this->deeq_nc(is, iat, ih, jh) = this->dvan(it, ih, jh);
                                this->deeq_nc(is, iat, jh, ih) = this->dvan(it, ih, jh);
                            }
                            else if (is == 1)
                            {
                                this->deeq_nc(is, iat, ih, jh) = std::complex<double>(0.0, 0.0);
                                this->deeq_nc(is, iat, jh, ih) = std::complex<double>(0.0, 0.0);
                            }
                            else if (is == 2)
                            {
                                this->deeq_nc(is, iat, ih, jh) = std::complex<double>(0.0, 0.0);
                                this->deeq_nc(is, iat, jh, ih) = std::complex<double>(0.0, 0.0);
                            }
                            else if (is == 3)
                            {
                                this->deeq_nc(is, iat, ih, jh) = this->dvan(it, ih, jh);
                                this->deeq_nc(is, iat, jh, ih) = this->dvan(it, ih, jh);
                            }
                        }
                        else
                        {
                            this->deeq(is, iat, ih, jh) = this->dvan(it, ih, jh);
                            this->deeq(is, iat, jh, ih) = this->dvan(it, ih, jh);
                            // in most of pseudopotential files, number of projections of one orbital is only one,
                            // which lead to diagonal matrix of dion
                            // when number larger than 1, non-diagonal dion should be calculated.
                            if (ih != jh && std::fabs(this->deeq(is, iat, ih, jh)) > 0.0)
                            {
                                this->multi_proj = true;
                            }
                        }
                    }
                }
            }
        }
    }
    else
    {
        newq(veff, rho_basis, cell);

        for (int iat = 0; iat < cell.nat; iat++)
        {
            int it = cell.iat2it[iat];
            if (GlobalV::NONCOLIN)
            {
                if (cell.atoms[it].ncpp.has_so)
                {
                    this->newd_so(iat, cell);
                }
                else
                {
                    this->newd_nc(iat, cell);
                }
            }
            else
            {
                for (int is = 0; is < GlobalV::NSPIN; is++)
                {
                    for (int ih = 0; ih < cell.atoms[it].ncpp.nh; ih++)
                    {
                        for (int jh = ih; jh < cell.atoms[it].ncpp.nh; jh++)
                        {
                            deeq(is, iat, ih, jh) += this->dvan(it, ih, jh);
                            deeq(is, iat, jh, ih) = deeq(is, iat, ih, jh);
                        }
                    }
                }
            }
        }
    }
    if (GlobalV::device_flag == "gpu") {
        if (GlobalV::precision_flag == "single") {
            castmem_d2s_h2d_op()(gpu_ctx,
                                 cpu_ctx,
                                 this->s_deeq,
                                 this->deeq.ptr,
                                 GlobalV::NSPIN * cell.nat * this->nhm * this->nhm);
            castmem_z2c_h2d_op()(gpu_ctx,
                                 cpu_ctx,
                                 this->c_deeq_nc,
                                 this->deeq_nc.ptr,
                                 GlobalV::NSPIN * cell.nat * this->nhm * this->nhm);
        }
        else {
            syncmem_z2z_h2d_op()(gpu_ctx,
                                 cpu_ctx,
                                 this->z_deeq_nc,
                                 this->deeq_nc.ptr,
                                 GlobalV::NSPIN * cell.nat * this->nhm * this->nhm);
        }
        syncmem_d2d_h2d_op()(gpu_ctx,
                             cpu_ctx,
                             this->d_deeq,
                             this->deeq.ptr,
                             GlobalV::NSPIN * cell.nat * this->nhm * this->nhm);
    }
    else {
        if (GlobalV::precision_flag == "single") {
            castmem_d2s_h2h_op()(cpu_ctx,
                                 cpu_ctx,
                                 this->s_deeq,
                                 this->deeq.ptr,
                                 GlobalV::NSPIN * cell.nat * this->nhm * this->nhm);
            castmem_z2c_h2h_op()(cpu_ctx,
                                 cpu_ctx,
                                 this->c_deeq_nc,
                                 this->deeq_nc.ptr,
                                 GlobalV::NSPIN * cell.nat * this->nhm * this->nhm);
        }
        // There's no need to synchronize double precision pointers while in a CPU environment.
    }
}

void pseudopot_cell_vnl::newq(const ModuleBase::matrix& veff, const ModulePW::PW_Basis* rho_basis, UnitCell& cell)
{
    ModuleBase::TITLE("pseudopot_cell_vnl", "newq");

    const std::complex<double> ci_tpi = ModuleBase::IMAG_UNIT * ModuleBase::TWO_PI;
    double fact = 1.0;
    if (rho_basis->gamma_only)
    {
        fact = 2.0;
    }

    const int npw = rho_basis->npw;
    ModuleBase::matrix ylmk0(lmaxq * lmaxq, npw);
    ModuleBase::YlmReal::Ylm_Real(lmaxq * lmaxq, npw, rho_basis->gcar, ylmk0);

    double* qnorm = new double[npw];
    for (int ig = 0; ig < npw; ig++)
    {
        qnorm[ig] = rho_basis->gcar[ig].norm() * cell.tpiba;
    }

    // fourier transform of the total effective potential
    ModuleBase::ComplexMatrix vaux(GlobalV::NSPIN, npw);
    for (int is = 0; is < GlobalV::NSPIN; is++)
    {
        rho_basis->real2recip(&veff(is, 0), &vaux(is, 0));
    }

    for (int it = 0; it < cell.ntype; it++)
    {
        Atom_pseudo* upf = &cell.atoms[it].ncpp;
        if (upf->tvanp)
        {
            // nij = max number of (ih,jh) pairs per atom type
            int nij = upf->nh * (upf->nh + 1) / 2;
            ModuleBase::ComplexMatrix qg(nij, npw);

            // Compute and store Q(G) for this atomic species
            // (without structure factor)
            int ijh = 0;
            for (int ih = 0; ih < upf->nh; ih++)
            {
                for (int jh = ih; jh < upf->nh; jh++)
                {
                    radial_fft_q(npw, ih, jh, it, qnorm, ylmk0, &qg(ijh, 0));
                    ijh++;
                }
            }

            // Compute and store V(G) times the structure factor e^(-iG*tau)
            const int natom = cell.atoms[it].na;
            ModuleBase::ComplexMatrix aux(natom, npw);
            ModuleBase::matrix deeaux(natom, nij);
            for (int is = 0; is < GlobalV::NSPIN; is++)
            {
                for (int ia = 0; ia < natom; ia++)
                {
                    const ModuleBase::Vector3<double> tau = cell.atoms[it].tau[ia];
                    for (int ig = 0; ig < npw; ig++)
                    {
                        const ModuleBase::Vector3<double> g = rho_basis->gcar[ig];
                        const std::complex<double> phase = ci_tpi * (g * tau);
                        aux(ia, ig) = vaux(is, ig) * exp(phase);
                    }
                }
                // here we compute the integral Q*V for all atoms of this kind
                const char transa = 'C', transb = 'N';
                const double zero = 0.0;
                const int complex_npw = 2 * npw;
                double* qg_ptr = reinterpret_cast<double*>(qg.c);
                double* aux_ptr = reinterpret_cast<double*>(aux.c);

                dgemm_(&transa,
                       &transb,
                       &nij,
                       &natom,
                       &complex_npw,
                       &fact,
                       qg_ptr,
                       &complex_npw,
                       aux_ptr,
                       &complex_npw,
                       &zero,
                       deeaux.c,
                       &nij);
                // I'm not sure if this is correct for gamma_only
                if (rho_basis->gamma_only && rho_basis->ig_gge0 >= 0)
                {
                    const double neg = -1.0;
                    dger_(&nij, &natom, &neg, qg_ptr, &complex_npw, aux_ptr, &complex_npw, deeaux.c, &nij);
                }

                for (int ia = 0; ia < natom; ia++)
                {
                    int ijh = 0;
                    const int iat = cell.itia2iat(it, ia);
                    for (int ih = 0; ih < upf->nh; ih++)
                    {
                        for (int jh = ih; jh < upf->nh; jh++)
                        {
                            deeq(is, iat, ih, jh) = cell.omega * deeaux(ia, ijh);
                            if (jh > ih)
                            {
                                deeq(is, iat, jh, ih) = deeq(is, iat, ih, jh);
                            }
                            ijh++;
                        }
                    }
                }
            }
        }
    }

#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, deeq.ptr, deeq.getSize(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

    delete[] qnorm;
}

void pseudopot_cell_vnl::newd_so(const int& iat, UnitCell& cell)
{
    ModuleBase::TITLE("pseudopot_cell_vnl", "newd_so");

    const int it = cell.iat2it[iat];
    Atom_pseudo* upf = &cell.atoms[it].ncpp;
    int ijs = 0;
    for (int is1 = 0; is1 < 2; is1++)
    {
        for (int is2 = 0; is2 < 2; is2++)
        {
            for (int ih = 0; ih < upf->nh; ih++)
            {
                for (int jh = 0; jh < upf->nh; jh++)
                {
                    deeq_nc(ijs, iat, ih, jh) = dvan_so(ijs, it, ih, jh);

                    for (int kh = 0; kh < upf->nh; kh++)
                    {
                        for (int lh = 0; lh < upf->nh; lh++)
                        {
                            if (GlobalV::DOMAG)
                            {
                                deeq_nc(ijs, iat, ih, jh)
                                    += deeq(0, iat, kh, lh)
                                           * (soc.fcoef(it, is1, 0, ih, kh) * soc.fcoef(it, 0, is2, lh, jh)
                                              + soc.fcoef(it, is1, 1, ih, kh) * soc.fcoef(it, 1, is2, lh, jh))
                                       + deeq(1, iat, kh, lh)
                                             * (soc.fcoef(it, is1, 0, ih, kh) * soc.fcoef(it, 1, is2, lh, jh)
                                                + soc.fcoef(it, is1, 1, ih, kh) * soc.fcoef(it, 0, is2, lh, jh))
                                       + ModuleBase::NEG_IMAG_UNIT * deeq(2, iat, kh, lh)
                                             * (soc.fcoef(it, is1, 0, ih, kh) * soc.fcoef(it, 1, is2, lh, jh)
                                                - soc.fcoef(it, is1, 1, ih, kh) * soc.fcoef(it, 0, is2, lh, jh))
                                       + deeq(3, iat, kh, lh)
                                             * (soc.fcoef(it, is1, 0, ih, kh) * soc.fcoef(it, 0, is2, lh, jh)
                                                - soc.fcoef(it, is1, 1, ih, kh) * soc.fcoef(it, 1, is2, lh, jh));
                            }
                            else
                            {
                                deeq_nc(ijs, iat, ih, jh)
                                    += deeq(0, iat, kh, lh)
                                       * (soc.fcoef(it, is1, 0, ih, kh) * soc.fcoef(it, 0, is2, lh, jh)
                                          + soc.fcoef(it, is1, 1, ih, kh) * soc.fcoef(it, 1, is2, lh, jh));
                            }
                        }
                    }
                }
            }
            ijs++;
        }
    }
}

void pseudopot_cell_vnl::newd_nc(const int& iat, UnitCell& cell)
{
    ModuleBase::TITLE("pseudopot_cell_vnl", "newd_nc");

    const int it = cell.iat2it[iat];
    Atom_pseudo* upf = &cell.atoms[it].ncpp;

    for (int ih = 0; ih < upf->nh; ih++)
    {
        for (int jh = 0; jh < upf->nh; jh++)
        {
            if (GlobalV::LSPINORB)
            {
                deeq_nc(0, iat, ih, jh) = dvan_so(0, it, ih, jh) + deeq(0, iat, ih, jh) + deeq(3, iat, ih, jh);
                deeq_nc(3, iat, ih, jh) = dvan_so(3, it, ih, jh) + deeq(0, iat, ih, jh) - deeq(3, iat, ih, jh);
            }
            else
            {
                deeq_nc(0, iat, ih, jh) = dvan(it, ih, jh) + deeq(0, iat, ih, jh) + deeq(3, iat, ih, jh);
                deeq_nc(3, iat, ih, jh) = dvan(it, ih, jh) + deeq(0, iat, ih, jh) - deeq(3, iat, ih, jh);
            }
            deeq_nc(1, iat, ih, jh) = deeq(1, iat, ih, jh) + ModuleBase::NEG_IMAG_UNIT * deeq(2, iat, ih, jh);
            deeq_nc(2, iat, ih, jh) = deeq(1, iat, ih, jh) + ModuleBase::IMAG_UNIT * deeq(2, iat, ih, jh);
        }
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

    template void pseudopot_cell_vnl::getvnl<float, psi::DEVICE_CPU>(psi::DEVICE_CPU*,
                                                                     int const&,
                                                                     std::complex<float>*) const;
    template void pseudopot_cell_vnl::getvnl<double, psi::DEVICE_CPU>(psi::DEVICE_CPU*,
                                                                      int const&,
                                                                      std::complex<double>*) const;
#if defined(__CUDA) || defined(__ROCM)
    template void pseudopot_cell_vnl::getvnl<float, psi::DEVICE_GPU>(psi::DEVICE_GPU*,
                                                                     int const&,
                                                                     std::complex<float>*) const;
    template void pseudopot_cell_vnl::getvnl<double, psi::DEVICE_GPU>(psi::DEVICE_GPU*,
                                                                      int const&,
                                                                      std::complex<double>*) const;
#endif

    template void pseudopot_cell_vnl::radial_fft_q<float, psi::DEVICE_CPU>(psi::DEVICE_CPU*,
                                                                           const int,
                                                                           const int,
                                                                           const int,
                                                                           const int,
                                                                           const float*,
                                                                           const float*,
                                                                           std::complex<float>*) const;
    template void pseudopot_cell_vnl::radial_fft_q<double, psi::DEVICE_CPU>(psi::DEVICE_CPU*,
                                                                            const int,
                                                                            const int,
                                                                            const int,
                                                                            const int,
                                                                            const double*,
                                                                            const double*,
                                                                            std::complex<double>*) const;
#if defined(__CUDA) || defined(__ROCM)
    template void pseudopot_cell_vnl::radial_fft_q<float, psi::DEVICE_GPU>(psi::DEVICE_GPU*,
                                                                           const int,
                                                                           const int,
                                                                           const int,
                                                                           const int,
                                                                           const float*,
                                                                           const float*,
                                                                           std::complex<float>*) const;
    template void pseudopot_cell_vnl::radial_fft_q<double, psi::DEVICE_GPU>(psi::DEVICE_GPU*,
                                                                            const int,
                                                                            const int,
                                                                            const int,
                                                                            const int,
                                                                            const double*,
                                                                            const double*,
                                                                            std::complex<double>*) const;
#endif
