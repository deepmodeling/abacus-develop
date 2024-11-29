#include "wavefunc.h"

#include "module_base/memory.h"
#include "module_base/timer.h"
#include "module_hamilt_lcao/hamilt_lcaodft/wavefunc_in_pw.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_hamilt_pw/hamilt_pwdft/hamilt_pw.h"
#include "module_hsolver/diago_iter_assist.h"
#include "module_hsolver/kernels/math_kernel_op.h"
#include "module_io/read_wfc_pw.h"
#include "module_io/winput.h"
#include "module_parameter/parameter.h"
#include "module_psi/psi.h"

wavefunc::wavefunc()
{
}

wavefunc::~wavefunc()
{
    if (PARAM.inp.test_deconstructor)
    {
        std::cout << " ~wavefunc()" << std::endl;
    }
    if (this->irindex != nullptr)
    {
        delete[] this->irindex;
        this->irindex = nullptr;
    }
}

psi::Psi<std::complex<double>>* wavefunc::allocate(const int nkstot, const int nks, const int* ngk, const int npwx_in)
{
    ModuleBase::TITLE("wavefunc", "allocate");

    this->npwx = npwx_in;
    this->nkstot = nkstot;
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "npwx", npwx);

    assert(npwx > 0);
    assert(nks > 0);

    // allocate for kinetic energy

    // if use spin orbital, do not double nks but double allocate evc and wanf2.
    int prefactor = 1;
    if (PARAM.inp.nspin == 4)
    {
        prefactor = PARAM.globalv.npol; // added by zhengdy-soc
    }

    const int nks2 = nks;

    psi::Psi<std::complex<double>>* psi_out = nullptr;
    if (PARAM.inp.calculation == "nscf" && PARAM.inp.mem_saver == 1)
    {
        // initial psi rather than evc
        psi_out = new psi::Psi<std::complex<double>>(1, PARAM.inp.nbands, npwx * PARAM.globalv.npol, ngk);
        if (PARAM.inp.basis_type == "lcao_in_pw")
        {
            wanf2[0].create(PARAM.globalv.nlocal, npwx * PARAM.globalv.npol);

            // WARNING: put the sizeof() be the first to avoid the overflow of the multiplication of int
            const size_t memory_cost
                = sizeof(std::complex<double>) * PARAM.globalv.nlocal * (PARAM.globalv.npol * npwx);

            std::cout << " Memory for wanf2 (MB): " << static_cast<double>(memory_cost) / 1024.0 / 1024.0 << std::endl;
            ModuleBase::Memory::record("WF::wanf2", memory_cost);
        }

        // WARNING: put the sizeof() be the first to avoid the overflow of the multiplication of int
        const size_t memory_cost = sizeof(std::complex<double>) * PARAM.inp.nbands * (PARAM.globalv.npol * npwx);

        std::cout << " MEMORY FOR PSI (MB)  : " << static_cast<double>(memory_cost) / 1024.0 / 1024.0 << std::endl;
        ModuleBase::Memory::record("Psi_PW", memory_cost);
    }
    else if (PARAM.inp.basis_type != "pw")
    {
        if ((PARAM.inp.basis_type == "lcao" || PARAM.inp.basis_type == "lcao_in_pw") || winput::out_spillage == 2)
        { // for lcao_in_pw
            if (this->wanf2 != nullptr)
            {
                delete[] this->wanf2;
            }
            this->wanf2 = new ModuleBase::ComplexMatrix[nks2];

            for (int ik = 0; ik < nks2; ik++)
            {
                this->wanf2[ik].create(PARAM.globalv.nlocal, npwx * PARAM.globalv.npol);
            }

            // WARNING: put the sizeof() be the first to avoid the overflow of the multiplication of int
            const size_t memory_cost
                = sizeof(std::complex<double>) * nks2 * PARAM.globalv.nlocal * (npwx * PARAM.globalv.npol);

            std::cout << " Memory for wanf2 (MB): " << static_cast<double>(memory_cost) / 1024.0 / 1024.0 << std::endl;
            ModuleBase::Memory::record("WF::wanf2", memory_cost);
        }
    }
    else
    {
        // initial psi rather than evc
        psi_out = new psi::Psi<std::complex<double>>(nks2, PARAM.inp.nbands, npwx * PARAM.globalv.npol, ngk);

        // WARNING: put the sizeof() be the first to avoid the overflow of the multiplication of int
        const size_t memory_cost = sizeof(std::complex<double>) * nks2 * PARAM.inp.nbands * (PARAM.globalv.npol * npwx);

        std::cout << " MEMORY FOR PSI (MB)  : " << static_cast<double>(memory_cost) / 1024.0 / 1024.0 << std::endl;
        ModuleBase::Memory::record("Psi_PW", memory_cost);
    }
    return psi_out;

    // showMemStats();
}

//===================================================================
// This routine computes an estimate of the start_ wavefunctions
// from superposition of atomic wavefunctions or random wave functions.
//===================================================================
void wavefunc::wfcinit(psi::Psi<std::complex<double>>* psi_in, ModulePW::PW_Basis_K* wfc_basis)
{
    ModuleBase::TITLE("wavefunc", "wfcinit");
    ModuleBase::timer::tick("wavefunc", "wfcinit");
    if (PARAM.inp.basis_type == "pw")
    {
        if (this->irindex != nullptr)
        {
            delete[] this->irindex;
        }
        this->irindex = new int[wfc_basis->fftnxy];
        wfc_basis->getfftixy2is(this->irindex);
#if defined(__CUDA) || defined(__ROCM)
        if (PARAM.inp.device == "gpu")
        {
            wfc_basis->get_ig2ixyz_k();
        }
#endif
    }
    ModuleBase::timer::tick("wavefunc", "wfcinit");
    return;
}

int wavefunc::get_starting_nw() const
{
    if (PARAM.inp.init_wfc == "file")
    {
        return PARAM.inp.nbands;
    }
    else if (PARAM.inp.init_wfc.substr(0, 6) == "atomic")
    {
        if (GlobalC::ucell.natomwfc >= PARAM.inp.nbands)
        {
            if (PARAM.inp.test_wf)
            {
                GlobalV::ofs_running << " Start wave functions are all pseudo atomic wave functions." << std::endl;
            }
        }
        else
        {
            if (PARAM.inp.test_wf)
            {
                GlobalV::ofs_running << " Start wave functions are atomic + "
                                     << PARAM.inp.nbands - GlobalC::ucell.natomwfc << " random wave functions."
                                     << std::endl;
            }
        }
        return std::max(GlobalC::ucell.natomwfc, PARAM.inp.nbands);
    }
    else if (PARAM.inp.init_wfc == "random")
    {
        if (PARAM.inp.test_wf)
        {
            GlobalV::ofs_running << " Start wave functions are all random." << std::endl;
        }
        return PARAM.inp.nbands;
    }
    else
    {
        throw std::runtime_error("wavefunc::get_starting_nw. Don't know what to do! Please Check source code! "
                                 + ModuleBase::GlobalFunc::TO_STRING(__FILE__) + " line "
                                 + ModuleBase::GlobalFunc::TO_STRING(__LINE__)); // Peize Lin change 2019-05-01
        // ModuleBase::WARNING_QUIT("get_starting_nw","Don't know what to do! Please Check source code!");
    }
}

namespace hamilt
{

template <>
void diago_PAO_in_pw_k2(const base_device::DEVICE_CPU* ctx,
                        const int& ik,
                        psi::Psi<std::complex<float>, base_device::DEVICE_CPU>& wvf,
                        ModulePW::PW_Basis_K* wfc_basis,
                        wavefunc* p_wf,
                        const ModuleBase::realArray& tab_at,
                        const int& lmaxkb,
                        hamilt::Hamilt<std::complex<float>, base_device::DEVICE_CPU>* phm_in)
{
    // TODO float func
}

template <>
void diago_PAO_in_pw_k2(const base_device::DEVICE_CPU* ctx,
                        const int& ik,
                        psi::Psi<std::complex<double>, base_device::DEVICE_CPU>& wvf,
                        ModulePW::PW_Basis_K* wfc_basis,
                        wavefunc* p_wf,
                        const ModuleBase::realArray& tab_at,
                        const int& lmaxkb,
                        hamilt::Hamilt<std::complex<double>, base_device::DEVICE_CPU>* phm_in)
{
    ModuleBase::TITLE("wavefunc", "diago_PAO_in_pw_k2");

    const int nbasis = wvf.get_nbasis();
    const int nbands = wvf.get_nbands();
    const int current_nbasis = wfc_basis->npwk[ik];

    if (PARAM.inp.init_wfc == "file")
    {
        ModuleBase::ComplexMatrix wfcatom(nbands, nbasis);
        std::stringstream filename;
        int ik_tot = K_Vectors::get_ik_global(ik, p_wf->nkstot);
        filename << PARAM.globalv.global_readin_dir << "WAVEFUNC" << ik_tot + 1 << ".dat";
        ModuleIO::read_wfc_pw(filename.str(), wfc_basis, ik, p_wf->nkstot, wfcatom);

        if (PARAM.inp.ks_solver == "cg")
        {
            std::vector<double> etfile(nbands, 0.0);
            if (phm_in != nullptr)
            {
                hsolver::DiagoIterAssist<std::complex<double>>::diagH_subspace_init(phm_in,
                                                                                    wfcatom.c,
                                                                                    wfcatom.nr,
                                                                                    wfcatom.nc,
                                                                                    wvf,
                                                                                    etfile.data());
                return;
            }
            else
            {
                ModuleBase::WARNING_QUIT("wavefunc", "Hamiltonian does not exist!");
            }
        }

        assert(nbands <= wfcatom.nr);
        for (int ib = 0; ib < nbands; ib++)
        {
            for (int ig = 0; ig < nbasis; ig++)
            {
                wvf(ib, ig) = wfcatom(ib, ig);
            }
        }
    }
    else if (PARAM.inp.init_wfc == "random"
             || (PARAM.inp.init_wfc.substr(0, 6) == "atomic" && GlobalC::ucell.natomwfc == 0))
    {
        p_wf->random(wvf.get_pointer(), 0, nbands, ik, wfc_basis);

        if (PARAM.inp.ks_solver == "cg")
        {
            std::vector<double> etrandom(nbands, 0.0);

            if (phm_in != nullptr)
            {
                hsolver::DiagoIterAssist<std::complex<double>>::diagH_subspace(phm_in, wvf, wvf, etrandom.data());
                return;
            }
            else
            {
                ModuleBase::WARNING_QUIT("wavefunc", "Hamiltonian does not exist!");
            }
        }
    }
    else if (PARAM.inp.init_wfc.substr(0, 6) == "atomic")
    {
        const int starting_nw = p_wf->get_starting_nw();
        if (starting_nw == 0)
        {
            return;
        }
        assert(starting_nw > 0);

        ModuleBase::ComplexMatrix wfcatom(starting_nw, nbasis);

        if (PARAM.inp.test_wf)
        {
            ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "starting_nw", starting_nw);
        }

        p_wf->atomic_wfc(ik,
                         current_nbasis,
                         GlobalC::ucell.lmax_ppwf,
                         lmaxkb,
                         wfc_basis,
                         wfcatom,
                         tab_at,
                         PARAM.globalv.nqx,
                         PARAM.globalv.dq);

        if (PARAM.inp.init_wfc == "atomic+random"
            && starting_nw == GlobalC::ucell.natomwfc) // added by qianrui 2021-5-16
        {
            p_wf->atomicrandom(wfcatom, 0, starting_nw, ik, wfc_basis);
        }

        //====================================================
        // If not enough atomic wfc are available, complete
        // with random wfcs
        //====================================================
        p_wf->random(wfcatom.c, GlobalC::ucell.natomwfc, nbands, ik, wfc_basis);

        // (7) Diago with cg method.
        // if(GlobalV::DIAGO_TYPE == "cg") xiaohui modify 2013-09-02
        if (PARAM.inp.ks_solver == "cg") // xiaohui add 2013-09-02
        {
            std::vector<double> etatom(starting_nw, 0.0);
            if (phm_in != nullptr)
            {
                hsolver::DiagoIterAssist<std::complex<double>>::diagH_subspace_init(phm_in,
                                                                                    wfcatom.c,
                                                                                    wfcatom.nr,
                                                                                    wfcatom.nc,
                                                                                    wvf,
                                                                                    etatom.data());
                return;
            }
            else
            {
                ModuleBase::WARNING_QUIT("wavefunc", "Hamiltonian does not exist!");
            }
        }

        assert(nbands <= wfcatom.nr);
        for (int ib = 0; ib < nbands; ib++)
        {
            for (int ig = 0; ig < nbasis; ig++)
            {
                wvf(ib, ig) = wfcatom(ib, ig);
            }
        }
    }
}

#if ((defined __CUDA) || (defined __ROCM))

template <>
void diago_PAO_in_pw_k2(const base_device::DEVICE_GPU* ctx,
                        const int& ik,
                        psi::Psi<std::complex<float>, base_device::DEVICE_GPU>& wvf,
                        ModulePW::PW_Basis_K* wfc_basis,
                        wavefunc* p_wf,
                        const ModuleBase::realArray& tab_at,
                        const int& lmaxkb,
                        hamilt::Hamilt<std::complex<float>, base_device::DEVICE_GPU>* phm_in)
{
    // TODO float 
}

template <>
void diago_PAO_in_pw_k2(const base_device::DEVICE_GPU* ctx,
                        const int& ik,
                        psi::Psi<std::complex<double>, base_device::DEVICE_GPU>& wvf,
                        ModulePW::PW_Basis_K* wfc_basis,
                        wavefunc* p_wf,
                        const ModuleBase::realArray& tab_at,
                        const int& lmaxkb,
                        hamilt::Hamilt<std::complex<double>, base_device::DEVICE_GPU>* phm_in)
{
    ModuleBase::TITLE("wavefunc", "diago_PAO_in_pw_k2");

    const int nbasis = wvf.get_nbasis();
    const int nbands = wvf.get_nbands();
    const int current_nbasis = wfc_basis->npwk[ik];
    int starting_nw = nbands;

    ModuleBase::ComplexMatrix wfcatom(nbands, nbasis);
    if (PARAM.inp.init_wfc == "file")
    {
        std::stringstream filename;
        int ik_tot = K_Vectors::get_ik_global(ik, p_wf->nkstot);
        filename << PARAM.globalv.global_readin_dir << "WAVEFUNC" << ik_tot + 1 << ".dat";
        ModuleIO::read_wfc_pw(filename.str(), wfc_basis, ik, p_wf->nkstot, wfcatom);
    }

    starting_nw = p_wf->get_starting_nw();
    if (starting_nw == 0)
        return;
    assert(starting_nw > 0);
    wfcatom.create(starting_nw, nbasis); // added by zhengdy-soc
    if (PARAM.inp.test_wf)
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "starting_nw", starting_nw);
    if (PARAM.inp.init_wfc.substr(0, 6) == "atomic")
    {
        p_wf->atomic_wfc(ik,
                         current_nbasis,
                         GlobalC::ucell.lmax_ppwf,
                         lmaxkb,
                         wfc_basis,
                         wfcatom,
                         tab_at,
                         PARAM.globalv.nqx,
                         PARAM.globalv.dq);
        if (PARAM.inp.init_wfc == "atomic+random"
            && starting_nw == GlobalC::ucell.natomwfc) // added by qianrui 2021-5-16
        {
            p_wf->atomicrandom(wfcatom, 0, starting_nw, ik, wfc_basis);
        }

        //====================================================
        // If not enough atomic wfc are available, complete
        // with random wfcs
        //====================================================
        p_wf->random(wfcatom.c, GlobalC::ucell.natomwfc, nbands, ik, wfc_basis);
    }
    else if (PARAM.inp.init_wfc == "random")
    {
        p_wf->random(wfcatom.c, 0, nbands, ik, wfc_basis);
    }

    std::complex<double>* z_wfcatom = nullptr;
    if (PARAM.inp.ks_solver != "bpcg")
    {
        // store wfcatom on the GPU
        resmem_zd_op()(gpu_ctx, z_wfcatom, wfcatom.nr * wfcatom.nc);
        syncmem_z2z_h2d_op()(gpu_ctx, cpu_ctx, z_wfcatom, wfcatom.c, wfcatom.nr * wfcatom.nc);
    }
    if (PARAM.inp.ks_solver == "cg") // xiaohui add 2013-09-02
    {
        // (7) Diago with cg method.
        if (phm_in != nullptr)
        {
            std::vector<double> etatom(starting_nw, 0.0);
            hsolver::DiagoIterAssist<std::complex<double>, base_device::DEVICE_GPU>::diagH_subspace_init(phm_in,
                                                                                                         z_wfcatom,
                                                                                                         wfcatom.nr,
                                                                                                         wfcatom.nc,
                                                                                                         wvf,
                                                                                                         etatom.data());
        }
        else
        {
            // this diagonalization method is obsoleted now
            // GlobalC::hm.diagH_subspace(ik ,starting_nw, nbands, wfcatom, wfcatom, etatom.data());
        }
    }
    else if (PARAM.inp.ks_solver == "dav" || PARAM.inp.ks_solver == "dav_subspace")
    {
        assert(nbands <= wfcatom.nr);
        // replace by haozhihan 2022-11-23
        hsolver::matrixSetToAnother<std::complex<double>, base_device::DEVICE_GPU>()(gpu_ctx,
                                                                                     nbands,
                                                                                     z_wfcatom,
                                                                                     wfcatom.nc,
                                                                                     &wvf(0, 0),
                                                                                     nbasis);
    }
    else if (PARAM.inp.ks_solver == "bpcg")
    {
        syncmem_z2z_h2d_op()(gpu_ctx, cpu_ctx, &wvf(0, 0), wfcatom.c, wfcatom.nr * wfcatom.nc);
    }

    if (PARAM.inp.ks_solver != "bpcg")
    {
        delmem_zd_op()(gpu_ctx, z_wfcatom);
    }
}

#endif

} // namespace hamilt

// LiuXh add a new function here,
// 20180515
void wavefunc::init_after_vc(const int nks)
{
    ModuleBase::TITLE("wavefunc", "init");
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "npwx", this->npwx);

    assert(this->npwx > 0);
    assert(nks > 0);
    assert(PARAM.inp.nbands > 0);

    const int nks2 = nks;
    const int nbasis = this->npwx * PARAM.globalv.npol;

    if ((PARAM.inp.basis_type == "lcao" || PARAM.inp.basis_type == "lcao_in_pw") || winput::out_spillage == 2)
    {
        if (wanf2 != nullptr)
        {
            delete[] wanf2;
        }
        this->wanf2 = new ModuleBase::ComplexMatrix[nks2];
        for (int ik = 0; ik < nks2; ik++)
        {
            this->wanf2[ik].create(PARAM.globalv.nlocal, nbasis);
        }
    }

    if (PARAM.inp.test_wf)
    {
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "psi allocation", "Done");
        if (PARAM.inp.basis_type == "lcao" || PARAM.inp.basis_type == "lcao_in_pw")
        {
            ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "wanf2 allocation", "Done");
        }
    }

    return;
}
