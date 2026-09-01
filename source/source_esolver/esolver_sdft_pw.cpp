#include "esolver_sdft_pw.h"

#include "source_base/global_variable.h"
#include "source_base/memory_recorder.h"
#include "source_estate/module_charge/symm_rho.h"
#include "source_hsolver/diago_iter_assist.h"
#include "source_hsolver/diago_params.h"
#include "source_io/module_parameter/parameter.h"
#include "source_pw/module_stodft/sto_dos.h"
#include "source_pw/module_stodft/sto_elecond.h"
#include "source_pw/module_stodft/sto_forces.h"
#include "source_pw/module_stodft/sto_stress_pw.h"

#include <algorithm>
#include <fstream>

namespace ModuleESolver
{

template <typename T, typename Device>
ESolver_SDFT_PW<T, Device>::ESolver_SDFT_PW()
{
    this->classname = "ESolver_SDFT_PW";
    this->basisname = "PW";
}

template <typename T, typename Device>
ESolver_SDFT_PW<T, Device>::~ESolver_SDFT_PW()
{
    //****************************************************
    // do not add any codes in this deconstructor funcion
    //****************************************************
}

template <typename T, typename Device>
void ESolver_SDFT_PW<T, Device>::before_all_runners(BaseCell& basecell, const Input_para& inp)
{
    basecell.require_kind(BaseCell::Kind::unit_cell, __FUNCTION__);
    UnitCell& ucell = static_cast<UnitCell&>(basecell);

    // 1) initialize parameters from int Input class
    this->nche_sto = inp.nche_sto;
    this->method_sto = inp.method_sto;
    this->stoche.init(inp.nche_sto, inp.method_sto, inp.emax_sto, inp.emin_sto);

    // 2) run "before_all_runners" in ESolver_KS
    ESolver_KS_PW<T, Device>::before_all_runners(ucell, inp);

    // 3) initialize the stochastic wave functions
    this->stowf.init(&this->kv, this->pw_wfc->npwk_max);
    if (inp.nbands_sto != 0)
    {
        if (inp.initsto_ecut < inp.ecutwfc)
        {
            this->stowf.init_sto_orbitals(inp.seed_sto);
        }
        else
        {
            this->stowf.init_sto_orbitals_Ecut(inp.seed_sto, this->kv, *this->pw_wfc, inp.initsto_ecut);
        }
    }
    else
    {
        this->stowf.init_com_orbitals();
    }
    if (this->method_sto == 2)
    {
        this->stowf.allocate_chiallorder(this->nche_sto);
    }
    this->stowf.sync_chi0();

    // 4) allocate spaces for \sqrt(f(H))|chi> and |\tilde{chi}>
    size_t size = stowf.chi0->size();
    this->stowf.shchi = new psi::Psi<T, Device>(this->kv.get_nks(),
                                                this->stowf.nchip_max,
                                                this->pw_wfc->npwk_max,
                                                this->kv.ngk,
                                                true);
    ModuleBase::Memory::record("SDFT::shchi", size * sizeof(T));

    if (inp.nbands > 0)
    {
        this->stowf.chiortho = new psi::Psi<T, Device>(this->kv.get_nks(),
                                                       this->stowf.nchip_max,
                                                       this->pw_wfc->npwk_max,
                                                       this->kv.ngk,
                                                       true);
        ModuleBase::Memory::record("SDFT::chiortho", size * sizeof(T));
    }

    return;
}

template <typename T, typename Device>
void ESolver_SDFT_PW<T, Device>::before_scf(UnitCell& ucell, const int istep)
{
    ModuleBase::TITLE("ESolver_SDFT_PW", "before_scf");
    ModuleBase::timer::start("ESolver_SDFT_PW", "before_scf");

    ESolver_KS_PW<T, Device>::before_scf(ucell, istep);
    delete reinterpret_cast<hamilt::HamiltPW<double>*>(this->p_hamilt);
    this->p_hamilt = new hamilt::HamiltSdftPW<T, Device>(this->pelec->pot,
                                                         this->pw_wfc,
                                                         &this->kv,
                                                         &this->ppcell,
                                                         &ucell,
                                                         PARAM.globalv.npol,
                                                         &this->stoche.emin_sto,
                                                         &this->stoche.emax_sto);
    this->p_hamilt_sto = static_cast<hamilt::HamiltSdftPW<T, Device>*>(this->p_hamilt);

    if (istep > 0 && this->inp_->nbands_sto != 0 && this->inp_->initsto_freq > 0 && istep % this->inp_->initsto_freq == 0)
    {
        this->stowf.update_sto_orbitals(this->inp_->seed_sto);
    }

    ModuleBase::timer::end("ESolver_SDFT_PW", "before_scf");
}

template <typename T, typename Device>
void ESolver_SDFT_PW<T, Device>::iter_finish(UnitCell& ucell, const int istep, int& iter, bool& conv_esolver)
{
    // call iter_finish() of ESolver_KS
    ESolver_KS::iter_finish(ucell, istep, iter, conv_esolver);
}

template <typename T, typename Device>
void ESolver_SDFT_PW<T, Device>::after_scf(UnitCell& ucell, const int istep, const bool conv_esolver)
{
    ModuleBase::TITLE("ESolver_SDFT_PW", "after_scf");
    ModuleBase::timer::start("ESolver_SDFT_PW", "after_scf");

    // 1) call after_scf() of ESolver_KS_PW
    ESolver_KS_PW<T, Device>::after_scf(ucell, istep, conv_esolver);

    ModuleBase::timer::end("ESolver_SDFT_PW", "after_scf");
}

template <typename T, typename Device>
void ESolver_SDFT_PW<T, Device>::hamilt2rho_single(UnitCell& ucell, int istep, int iter, double ethr)
{
    ModuleBase::TITLE("ESolver_SDFT_PW", "hamilt2rho");
    ModuleBase::timer::start("ESolver_SDFT_PW", "hamilt2rho");

    // reset energy
    this->pelec->f_en.eband = 0.0;
    this->pelec->f_en.demet = 0.0;

    // setup diagonalization parameters for SDFT
    hsolver::setup_diago_params_sdft<T, Device>(istep, iter, ethr, *this->inp_);

    bool skip_charge = this->inp_->calculation == "nscf" ? true : false;

    // hsolver only exists in this function
    hsolver::HSolverPW_SDFT<T, Device> hsolver_pw_sdft_obj(&this->kv,
                                                           this->pw_wfc,
                                                           this->stowf,
                                                           this->stoche,
                                                           this->p_hamilt_sto,
                                                           this->inp_->calculation,
                                                           this->inp_->basis_type,
                                                           this->inp_->ks_solver,
                                                           PARAM.globalv.use_uspp,
                                                           this->inp_->nspin,
                                                           hsolver::DiagoIterAssist<T, Device>::SCF_ITER,
                                                           hsolver::DiagoIterAssist<T, Device>::PW_DIAG_NMAX,
                                                           hsolver::DiagoIterAssist<T, Device>::PW_DIAG_THR,
                                                           hsolver::DiagoIterAssist<T, Device>::need_subspace,
                                                           this->inp_->nbands,
                                                           this->inp_->diago_smooth_ethr,
                                                           this->inp_->pw_diag_ndim,
                                                           this->inp_->diag_subspace,
                                                           this->inp_->nb2d,
                                                           PARAM.globalv.ks_run,
                                                           PARAM.globalv.all_ks_run,
                                                           this->inp_->bndpar);

    hsolver_pw_sdft_obj.solve(ucell,
                              static_cast<hamilt::Hamilt<T, Device>*>(this->p_hamilt),
                              *this->stp.template get_psi_t<T, Device>(),
                              this->stp.psi_cpu[0],
                              this->pelec,
                              this->pw_wfc,
                              this->stowf,
                              istep,
                              iter,
                              GlobalV::ofs_running,
                              skip_charge);

    // set_diagethr need it
    this->esolver_KS_ne = hsolver_pw_sdft_obj.stoiter.KS_ne;

    if (PARAM.globalv.ks_run)
    {
        Symmetry_rho::symmetrize_rho(this->inp_->nspin, this->chr, this->pw_rho, ucell.symm);
        this->pelec->f_en.deband = this->pelec->cal_delta_eband(ucell);
    }
    else
    {
#ifdef __MPI
        if (ModuleSymmetry::Symmetry::symm_flag == 1)
        {
            MPI_Barrier(MPI_COMM_WORLD);
        }
#endif
    }
#ifdef __MPI
    MPI_Bcast(&(this->pelec->f_en.deband), 1, MPI_DOUBLE, 0, BP_WORLD);
#endif
    ModuleBase::timer::end("ESolver_SDFT_PW", "hamilt2rho");
}

template <typename T, typename Device>
double ESolver_SDFT_PW<T, Device>::cal_energy()
{
    return this->pelec->f_en.etot;
}

template <typename T, typename Device>
void ESolver_SDFT_PW<T, Device>::cal_force(BaseCell& basecell, ModuleBase::matrix& force)
{
    basecell.require_kind(BaseCell::Kind::unit_cell, __FUNCTION__);
    UnitCell& ucell = static_cast<UnitCell&>(basecell);

    Sto_Forces<double, Device> ff(ucell.nat);

    ff.cal_stoforce(force,
                    *this->pelec,
                    this->pw_rho,
                    &ucell.symm,
                    &this->sf,
                    &this->kv,
                    this->pw_wfc,
                    this->locpp,
                    this->ppcell,
                    ucell,
                    *this->stp.template get_psi_t<T, Device>(),
                    this->stowf);
}

template <typename T, typename Device>
void ESolver_SDFT_PW<T, Device>::cal_stress(BaseCell& basecell, ModuleBase::matrix& stress)
{
    basecell.require_kind(BaseCell::Kind::unit_cell, __FUNCTION__);
    UnitCell& ucell = static_cast<UnitCell&>(basecell);

    Sto_Stress_PW<double, Device> ss;
    ss.cal_stress(stress,
                  *this->pelec,
                  this->pw_rho,
                  &ucell.symm,
                  &this->sf,
                  &this->kv,
                  this->pw_wfc,
                  *this->stp.template get_psi_t<T, Device>(),
                  this->stowf,
                  &this->chr,
                  &this->locpp,
                  &this->ppcell,
                  ucell);
}

template <typename T, typename Device>
void ESolver_SDFT_PW<T, Device>::after_all_runners(BaseCell& basecell)
{
    basecell.require_kind(BaseCell::Kind::unit_cell, __FUNCTION__);
    UnitCell& ucell = static_cast<UnitCell&>(basecell);

    // 1) write down etot and eigenvalues (for MDFT) information
    ESolver_FP::after_all_runners(ucell);

    // 2) release memory
    if (this->method_sto == 2)
    {
        stowf.clean_chiallorder(); // release lots of memories
    }

    // 3) write down DOS
    if (this->inp_->out_dos)
    {
        if (!std::is_same<T, std::complex<double>>::value || !std::is_same<Device, base_device::DEVICE_CPU>::value)
        {
            ModuleBase::WARNING_QUIT("ESolver_SDFT_PW", "DOS does not support complex float or GPU yet.");
        }
        Sto_DOS<Real, Device> sto_dos(
            this->pw_wfc,
            &this->kv,
            this->pelec,
            reinterpret_cast<psi::Psi<std::complex<double>>*>(this->stp.psi_cpu),
            reinterpret_cast<hamilt::Hamilt<std::complex<double>>*>(this->p_hamilt),
            this->stoche,
            reinterpret_cast<Stochastic_WF<std::complex<double>, base_device::DEVICE_CPU>*>(&stowf));
        sto_dos.decide_param(this->inp_->dos_nche,
                             this->inp_->emin_sto,
                             this->inp_->emax_sto,
                             PARAM.globalv.dos_setemin,
                             PARAM.globalv.dos_setemax,
                             this->inp_->dos_emin_ev,
                             this->inp_->dos_emax_ev,
                             this->inp_->dos_scale);
        sto_dos.caldos(this->inp_->dos_sigma, this->inp_->dos_edelta_ev, this->inp_->npart_sto);
    }

    // 4) sKG cost memory, and it should be placed at the end of the program
    if (this->inp_->cal_cond)
    {
        Sto_EleCond<Real, Device> sto_elecond(
            &ucell,
            &this->kv,
            this->pelec,
            this->pw_wfc,
            this->stp.template get_psi_t<T, Device>(),
            &this->ppcell,
            static_cast<hamilt::Hamilt<std::complex<double>, Device>*>(this->p_hamilt),
            this->stoche,
            &stowf);
        sto_elecond.decide_nche(this->inp_->cond_dt, 1e-8, this->nche_sto, this->inp_->emin_sto, this->inp_->emax_sto);
        sto_elecond.sKG(this->inp_->cond_smear,
                        this->inp_->cond_fwhm,
                        this->inp_->cond_wcut,
                        this->inp_->cond_dw,
                        this->inp_->cond_dt,
                        this->inp_->cond_nonlocal,
                        this->inp_->npart_sto);
    }
}

// template class ESolver_SDFT_PW<std::complex<float>, base_device::DEVICE_CPU>;
template class ESolver_SDFT_PW<std::complex<double>, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
// template class ESolver_SDFT_PW<std::complex<float>, base_device::DEVICE_GPU>;
template class ESolver_SDFT_PW<std::complex<double>, base_device::DEVICE_GPU>;
#endif
} // namespace ModuleESolver
