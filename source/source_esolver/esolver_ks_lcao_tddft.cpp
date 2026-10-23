#include "esolver_ks_lcao_tddft.h"

#include "source_base/constants.h"
#include "source_lcao/module_rt/boundary_fix.h"

//----------------IO-----------------
#include "source_base/global_variable.h"
#include "source_io/module_ctrl/ctrl_output_td.h"
#include "source_io/module_efield/td_efield_io.h"
#include "source_io/module_efield/td_vector_pot_io.h"
#include "source_io/module_output/output_log.h"
#include "source_io/module_parameter/parameter.h"
#include "source_io/module_wf/read_wfc_nao.h"
//------LCAO HSolver ElecState-------
#include "source_estate/elecstate_tools.h"
#include "source_estate/module_charge/chg_atomic.h"
#include "source_estate/module_charge/chg_symm.h"
#include "source_estate/module_dm/dm_from_psi.h"
#include "source_estate/module_dm/edm_tddft.h"
#include "source_estate/module_pot/potential_new.h"
#include "source_estate/module_pot/td_field_manager.h"
#include "source_estate/write_init.h"
#include "source_hamilt/hamilt_hs_adapter.h"
#include "source_hsolver/hsolver_lcao.h"
#include "source_lcao/module_rt/evolve_elec.h"
#include "source_lcao/rho_tau_lcao.h"
#ifdef __EXX
#include "source_lcao/module_ri/exx_lri_interface.h"
#endif

namespace ModuleESolver
{

template <typename TR, typename Device>
ESolver_KS_LCAO_TDDFT<TR, Device>::ESolver_KS_LCAO_TDDFT()
{
    this->classname = "ESolver_rtTDDFT";
    this->basisname = "LCAO";
}

template <typename TR, typename Device>
ESolver_KS_LCAO_TDDFT<TR, Device>::~ESolver_KS_LCAO_TDDFT()
{
    if (psi_laststep != nullptr)
    {
        delete psi_laststep;
        psi_laststep = nullptr;
    }

    // Destroy operators before the TD context they reference.
    delete this->p_hamilt;
    this->p_hamilt = nullptr;

    if (td_p != nullptr)
    {
        delete td_p;
    }

    if (td_mg_ != nullptr)
    {
        delete td_mg_;
        td_mg_ = nullptr;
    }
}

template <typename TR, typename Device>
void ESolver_KS_LCAO_TDDFT<TR, Device>::before_all_runners(BaseCell& basecell, const Input_para& inp)
{
    basecell.require_kind(BaseCell::Kind::unitcell, __FUNCTION__);
    UnitCell& ucell = static_cast<UnitCell&>(basecell);
    const std::string& output_dir = PARAM.globalv.global_out_dir;

    // If the device is GPU, enable use_tensor; also enable use_lapack unless ks_solver is cuSOLVERMp.
    ct::DeviceType ct_device_type = ct::DeviceTypeToEnum<Device>::value;
    if (ct_device_type == ct::DeviceType::GpuDevice)
    {
        use_tensor = true;
        if (inp.ks_solver != "cusolvermp")
        {
            use_lapack = true;
        }
    }

    // Build the shared field state before the base ESolver registers potential
    // components. Both propagation gauges and the length-gauge potential must
    // observe the same electronic-step counter.
    td_field_manager_ = elecstate::create_td_field_manager(inp);
    if (inp.init_vecpot)
    {
        // An empty directory selects the calculation working directory, not OUT.<suffix>.
        const std::vector<ModuleBase::Vector3<double>> samples = ModuleIO::read_td_vector_pot("");
        td_field_manager_->set_A_samples(samples);
    }
    if (inp.mdp.md_restart)
    {
        td_field_manager_->read_restart(PARAM.globalv.global_readin_dir);
    }
    if (inp.out_efield && GlobalV::MY_RANK == 0)
    {
        ModuleIO::prepare_td_field_output(output_dir, *td_field_manager_, inp.mdp.md_restart);
    }
    if (inp.out_vecpot && GlobalV::MY_RANK == 0)
    {
        if (inp.init_vecpot)
        {
            ModuleIO::check_td_vector_pot_paths("", output_dir);
        }
        ModuleIO::prepare_td_vector_pot_output(output_dir, inp.mdp.md_restart);
    }

    // Run before_all_runners in ESolver_KS_LCAO
    ESolver_KS_LCAO<std::complex<double>, TR>::before_all_runners(ucell, inp);
    this->pelec->pot->set_td_field_manager(td_field_manager_);

    restart_step_ = inp.mdp.md_restart ? td_field_manager_->current_step() : 0;
    td_p = new TD_info(&ucell, this->pv, this->orb_, *td_field_manager_);
    totstep += restart_step_;

    if (this->inp_->init_wfc == "file")
    {
        if (!ModuleIO::read_wfc_nao(PARAM.globalv.global_readin_dir,
                                    this->pv,
                                    *(this->psi),
                                    this->pelec->ekb,
                                    this->pelec->wg,
                                    this->kv.ik2iktot,
                                    this->kv.get_nkstot(),
                                    this->inp_->nspin,
                                    this->inp_->init_wfc_file_format == "binary",
                                    0,
                                    restart_step_))
        {
            ModuleBase::WARNING_QUIT("ESolver_KS_LCAO_TDDFT", "Read electronic wavefunction from file failed!");
        }
    }
}

template <typename TR, typename Device>
void ESolver_KS_LCAO_TDDFT<TR, Device>::runner(BaseCell& basecell, const int istep)
{
    basecell.require_kind(BaseCell::Kind::unitcell, __FUNCTION__);
    UnitCell& ucell = static_cast<UnitCell&>(basecell);

    ModuleBase::TITLE("ESolver_KS_LCAO_TDDFT", "runner");
    ModuleBase::timer::start(this->classname, "runner");

    //----------------------------------------------------------------
    // 1) before_scf (preparation before the electronic-step loop)
    //----------------------------------------------------------------
    // Prepare the endpoint field before before_scf initializes the potential and constructs TD operators.
    const int first_estep = this->totstep + 1;
    if (this->inp_->td_stype == 0 || first_estep == 0)
    {
        td_field_manager_->prepare_sample(first_estep);
    }
    else
    {
        const int left_step = first_estep - 1;
        td_field_manager_->prepare_interval(left_step, first_estep);
    }
    td_p->update_field();
    this->before_scf(ucell, istep, first_estep); // From ESolver_KS_LCAO
    td_p->initialize_phase_hybrid(ucell, dynamic_cast<hamilt::HamiltLCAO<std::complex<double>, TR>*>(this->p_hamilt)->getHR());
    // Refresh existing entries too: phase initialization updates geometry but only creates missing keys.
    td_p->update_field();
    td_p->calculate_grad_overlap(this->pv, ucell, this->gd, this->orb_.cutoffs(), this->two_center_bundle_.overlap_orb.get());
    // Initialize the moving spatial gauge
    if (use_td_moving_gauge && this->td_mg_ == nullptr)
    {
        this->td_mg_ = new module_rt::TD_MovingGauge();
        auto* hamilt_lcao = dynamic_cast<hamilt::HamiltLCAO<std::complex<double>, TR>*>(this->p_hamilt);
        const hamilt::HContainer<TR>* sR_template = hamilt_lcao->getSR();
        this->td_mg_->init_DR(sR_template, &ucell, &this->pv, this->two_center_bundle_.overlap_orb.get());
    }

    if (this->inp_->td_stype == 2)
    {
        this->dmat.dm->cal_dmr_td(td_p->get_phase_hybrid(), td_p->endpoint_A(), -1);
    }
    else
    {
        this->dmat.dm->cal_dmr(-1);
    }
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "INIT SCF");

    // Initialize the velocity operator for electric-current calculation
    if (this->inp_->td_stype != 1 && TD_info::out_current == 1)
    {
        // initialize the velocity operator
        velocity_mat = new Velocity_op<TR>(&ucell, &(this->gd), &this->pv, this->orb_, this->two_center_bundle_.overlap_orb.get());
        // calculate velocity operator
        velocity_mat->calculate_grad_term();
        velocity_mat->calculate_vcomm_r();
    }
    int estep_max = (istep == 0 && !this->inp_->mdp.md_restart) ? 1 : this->inp_->estep_per_md;
    // Include the initial state before fixed-ion propagation; mohan changed md_nstep from 0 to 1, 2026-01-04.
    if (this->inp_->mdp.md_nstep == 1)
    {
        estep_max = this->inp_->estep_per_md + 1;
    }

    // Update the last-step matrices and wave functions if any atom crosses the boundary
    // Apply a phase correction to H, S, and psi to keep consistency when atoms cross periodic boundaries
    const size_t len_hs_ik = use_tensor && use_lapack ? PARAM.globalv.nlocal * PARAM.globalv.nlocal : this->pv.nloc;
    module_rt::reset_matrix_boundary(ucell, this->kv, &(this->pv), this->Hk_laststep, this->Sk_laststep, this->psi_laststep, len_hs_ik);

    for (int estep = 0; estep < estep_max; estep++)
    {
        // Advance the electronic-step index, including the initial state or restart replay.
        this->totstep++;
        this->print_step();
        // update At
        if (this->inp_->td_stype > 0)
        {
            if (this->totstep == 0)
            {
                td_field_manager_->prepare_sample(0);
            }
            else
            {
                const int left_step = this->totstep - 1;
                td_field_manager_->prepare_interval(left_step, this->totstep);
            }
            td_p->update_field();
            TD_info::evolve_once = true;
            // Replaying the saved endpoint rebuilds electronic state without appending another A sample.
            const bool replay_restart = this->inp_->mdp.md_restart && !restart_done;
            if (this->inp_->out_vecpot && GlobalV::MY_RANK == 0 && !replay_restart)
            {
                ModuleIO::write_td_vector_pot(PARAM.globalv.global_out_dir, this->totstep, td_field_manager_->A_right_ha());
            }
            ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "Cartesian vector potential Ax(t)", td_p->endpoint_A()[0]);
            ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "Cartesian vector potential Ay(t)", td_p->endpoint_A()[1]);
            ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "Cartesian vector potential Az(t)", td_p->endpoint_A()[2]);
        }

        if (this->inp_->td_stype == 0)
        {
            td_field_manager_->prepare_sample(this->totstep);
        }
        // Keep E output aligned with A by omitting the same replayed restart endpoint.
        const bool replay_restart = this->inp_->mdp.md_restart && !restart_done;
        if (this->inp_->out_efield && GlobalV::MY_RANK == 0 && !replay_restart)
        {
            ModuleIO::write_td_field_values(*td_field_manager_, PARAM.globalv.global_out_dir);
        }

        if (estep != 0)
        {
            this->CE.update_all_dis(ucell);
            const module_charge::AtomicRhoCfg atomic_rho_cfg_tddft{PARAM.inp.nelec,
                                                                   PARAM.inp.test_charge,
                                                                   PARAM.globalv.domag,
                                                                   PARAM.globalv.domag_z,
                                                                   GlobalV::ofs_warning};
            this->CE.extrapolate_charge(&this->Pgrid,
                                        ucell,
                                        &this->chr,
                                        *this->pw_rhod,
                                        &this->sf,
                                        GlobalV::ofs_running,
                                        GlobalV::ofs_warning,
                                        atomic_rho_cfg_tddft,
                                        PARAM.globalv.has_float_data);
            this->exx_nao.before_scf(ucell, this->kv, this->orb_, this->pv, this->p_chgmix, totstep, *this->inp_, this->exx_info_);
            ModuleIO::write_chg_init(ucell,
                                     this->Pgrid,
                                     this->chr,
                                     this->pelec->eferm,
                                     totstep,
                                     PARAM.globalv.global_out_dir,
                                     *this->inp_,
                                     PARAM.globalv.two_fermi,
                                     GlobalV::ofs_running);
            elecstate::init_scf(ucell,
                                this->Pgrid,
                                this->sf.strucFac,
                                this->locpp.numeric,
                                totstep,
                                PARAM.globalv.global_out_dir,
                                *this->inp_,
                                this->pelec);

            if (totstep <= this->inp_->td_tend + 1)
            {
                TD_info::evolve_once = true;
            }
        }
        //----------------------------------------------------------------
        // 2) SCF iterations
        //----------------------------------------------------------------
        bool conv_esolver = false;
        this->niter = this->maxniter;
        this->diag_ethr = this->inp_->pw_diag_thr;
        for (int iter = 1; iter <= this->maxniter; ++iter)
        {
            ModuleIO::write_head_td(GlobalV::ofs_running, istep, totstep, iter, this->basisname);

            // 2.1) Initialization of the SCF iteration
            this->iter_init(ucell, totstep, iter); // From ESolver_KS_LCAO

            // 2.2) Use Hamiltonian to obtain charge density
            this->hamilt2rho(ucell, totstep, iter, this->diag_ethr); // From ESolver_KS

            // 2.3) Finish the SCF iteration
            this->iter_finish(ucell, totstep, estep, estep_max, iter, conv_esolver);

            // 2.4) Check convergence or density oscillation
            if (conv_esolver || this->oscillate_esolver)
            {
                this->niter = iter;
                if (this->oscillate_esolver)
                {
                    std::cout << " !! Density oscillation is found, STOP HERE !!" << std::endl;
                }
                break;
            }
        } // end SCF iterations

        //----------------------------------------------------------------
        // 3) after_scf
        //----------------------------------------------------------------
        this->after_scf(ucell, totstep, conv_esolver);
        if (this->inp_->out_freq_td > 0 && totstep % this->inp_->out_freq_td == 0 && GlobalV::MY_RANK == 0)
        {
            td_field_manager_->write_restart(PARAM.globalv.global_out_dir);
        }
        if (!restart_done && this->inp_->mdp.md_restart)
        {
            restart_done = true;
            estep += restart_step_ % this->inp_->estep_per_md;
            if (estep == 0)
            {
                break;
            }
            // mohan add 2026-01-04, change md_nstep!=0 to md_nstep!=1
            if (this->inp_->mdp.md_nstep != 1)
            {
                estep -= 1;
            }
        }
    }

    if (this->inp_->td_stype != 1 && TD_info::out_current == 1)
    {
        delete velocity_mat;
    }

    ModuleBase::timer::end(this->classname, "runner");
    return;
}

// Output electronic step information using a 1-based step number
template <typename TR, typename Device>
void ESolver_KS_LCAO_TDDFT<TR, Device>::print_step()
{
    std::cout << " -------------------------------------------" << std::endl;
    std::cout << " STEP OF ELECTRON EVOLVE : " << unsigned(totstep) + 1 << std::endl;
    std::cout << " -------------------------------------------" << std::endl;
}

template <typename TR, typename Device>
void ESolver_KS_LCAO_TDDFT<TR, Device>::hamilt2rho_single(UnitCell& ucell, const int istep, const int iter, const double ethr)
{
    std::vector<module_rt::EijOutputInfo> eij_output(this->kv.get_nks());
    for (int ik = 0; ik < this->kv.get_nks(); ++ik)
    {
        module_rt::EijOutputInfo& info = eij_output[ik];
        info.threshold = this->inp_->td_print_eij;
        info.electronic_step = istep + 1;
        info.scf_iteration = iter;
        info.kpoint = this->kv.ik2iktot[ik] + 1;
        info.nspin = this->inp_->nspin;
        info.spin = this->kv.isk[ik] + 1;
        if (info.nspin == 2)
        {
            info.kpoint -= this->kv.isk[ik] * (this->kv.get_nkstot() / 2);
        }
    }

    // Update the moving spatial gauge
    if (use_td_moving_gauge)
    {
        auto* hamilt_lcao = dynamic_cast<hamilt::HamiltLCAO<std::complex<double>, TR>*>(this->p_hamilt);
        const hamilt::HContainer<TR>* sR_template = hamilt_lcao->getSR();
        this->td_mg_->update_DR(sR_template, &ucell, &this->pv, this->two_center_bundle_.overlap_orb.get());
    }

    if (this->inp_->init_wfc == "file")
    {
        if (istep >= restart_step_ + 1)
        {
            module_rt::Evolve_elec<Device>::solve_psi(eij_output,
                                                      this->inp_->nbands,
                                                      PARAM.globalv.nlocal,
                                                      this->kv.get_nks(),
                                                      static_cast<hamilt::Hamilt<std::complex<double>>*>(this->p_hamilt),
                                                      this->pv,
                                                      this->psi,
                                                      this->psi_laststep,
                                                      this->Hk_laststep,
                                                      this->Sk_laststep,
                                                      this->pelec->ekb,
                                                      GlobalV::ofs_running,
                                                      this->inp_->propagator,
                                                      use_tensor,
                                                      use_lapack,
                                                      this->td_mg_,
                                                      &ucell,
                                                      this->kv.kvec_d,
                                                      use_td_moving_gauge);
        }
        this->weight_dm_rho(ucell);
    }
    else if (istep >= 1)
    {
        module_rt::Evolve_elec<Device>::solve_psi(eij_output,
                                                  this->inp_->nbands,
                                                  PARAM.globalv.nlocal,
                                                  this->kv.get_nks(),
                                                  static_cast<hamilt::Hamilt<std::complex<double>>*>(this->p_hamilt),
                                                  this->pv,
                                                  this->psi,
                                                  this->psi_laststep,
                                                  this->Hk_laststep,
                                                  this->Sk_laststep,
                                                  this->pelec->ekb,
                                                  GlobalV::ofs_running,
                                                  this->inp_->propagator,
                                                  use_tensor,
                                                  use_lapack,
                                                  this->td_mg_,
                                                  &ucell,
                                                  this->kv.kvec_d,
                                                  use_td_moving_gauge);
        this->weight_dm_rho(ucell);
    }
    else
    {
        // For the first step, do normal SCF calculation to get initial state
        this->pelec->f_en.eband = 0.0;
        this->pelec->f_en.demet = 0.0;
        if (this->psi != nullptr)
        {
            const bool skip_charge = false;
            hsolver::HSolverLCAO<std::complex<double>> hsolver_lcao_obj(&this->pv,
                                                                        this->inp_->ks_solver,
                                                                        PARAM.globalv.kpar_lcao,
                                                                        PARAM.globalv.nlocal,
                                                                        this->inp_->nbands,
                                                                        this->inp_->nelec,
                                                                        this->inp_->device == "gpu",
                                                                        GlobalV::NPROC,
                                                                        GlobalV::MY_RANK);
            hamilt::HamiltHSMatrix<std::complex<double>> hs(static_cast<hamilt::Hamilt<std::complex<double>>*>(this->p_hamilt));
            hsolver_lcao_obj.solve(hs, this->psi[0], this->pelec, *this->dmat.dm, this->chr, this->inp_->nspin, ucell.omega, skip_charge);
            if (this->inp_->td_stype == 2)
            {
                // Taoni add 2026-10-11
                // Ground-state diagonalization forms an ordinary DMR. Apply the
                // same hybrid phase convention as propagation before building rho.
                this->dmat.dm->cal_dmr_td(td_p->get_phase_hybrid(), td_p->endpoint_A(), -1);
                LCAO_domain::dm2rho(this->dmat.dm->get_dmr_vec(), this->inp_->nspin, &this->chr, this->inp_->nelec, ucell.omega, false);
            }
        }

        // Symmetrize the charge density only for the ground state
        module_charge::symmetrize_rho(this->inp_->nspin, this->chr, this->pw_rho, ucell.symm);
    }
#ifdef __EXX
    if (this->exx_info_.info_ri.real_number)
        this->exx_nao.exd->exx_hamilt2rho(*this->pelec, this->pv, iter);
    else
        this->exx_nao.exc->exx_hamilt2rho(*this->pelec, this->pv, iter);
#endif

    // Calculate delta energy
    this->pelec->f_en.deband = this->pelec->cal_delta_eband(ucell);
}

template <typename TR, typename Device>
void ESolver_KS_LCAO_TDDFT<TR, Device>::iter_finish(UnitCell& ucell,
                                                    const int istep,
                                                    const int estep,
                                                    const int estep_max,
                                                    int& iter,
                                                    bool& conv_esolver)
{
    // Print the occupation of each band at the first iteration of the first three electronic steps
    if (iter == 1 && istep <= 2)
    {
        GlobalV::ofs_running << " k-point  State   Occupations" << std::endl;
        GlobalV::ofs_running << std::setiosflags(std::ios::showpoint);
        GlobalV::ofs_running << std::left;
        GlobalV::ofs_running << std::setprecision(6);
        for (int ik = 0; ik < this->kv.get_nks(); ik++)
        {
            for (int ib = 0; ib < this->inp_->nbands; ib++)
            {
                GlobalV::ofs_running << " " << std::setw(9) << ik + 1 << std::setw(8) << ib + 1 << std::setw(12) << this->pelec->wg(ik, ib)
                                     << std::endl;
            }
        }
        GlobalV::ofs_running << std::endl;
    }

    ESolver_KS_LCAO<std::complex<double>, TR>::iter_finish(ucell, istep, iter, conv_esolver);

    // Store wave functions, Hamiltonian and overlap matrix, to be used in the next electronic step
    // Store when converged or when the maximum iteration count is reached
    bool force_save = conv_esolver || (iter == this->maxniter);
    this->store_h_s_psi(ucell, istep, iter, force_save);

    // Calculate the RT-TDDFT energy-density matrix at the last converged substep when required by td_edm and the gauge
    if (conv_esolver && estep == estep_max - 1 && istep >= (this->inp_->init_wfc == "file" ? 0 : 1) && this->inp_->td_edm == 0
        && this->inp_->td_stype != 2)
    {
        if (use_tensor && use_lapack)
        {
            module_dm::edm_tddft_lapack<Device>(this->pv,
                                                this->dmat,
                                                this->kv,
                                                static_cast<hamilt::Hamilt<std::complex<double>>*>(this->p_hamilt));
        }
        else
        {
            module_dm::edm_tddft(this->pv, this->dmat, this->kv, static_cast<hamilt::Hamilt<std::complex<double>>*>(this->p_hamilt));
        }
    }
}

template <typename TR, typename Device>
void ESolver_KS_LCAO_TDDFT<TR, Device>::store_h_s_psi(UnitCell& ucell, const int istep, const int iter, const bool conv_esolver)
{
    const int nbands = this->inp_->nbands;
    const int nlocal = PARAM.globalv.nlocal;

    // Store wave functions, Hamiltonian and overlap matrix when the caller requests a save
    if (conv_esolver)
    {
        if (this->psi_laststep == nullptr)
        {
            this->psi_laststep = new psi::Psi<std::complex<double>>(this->kv.get_nks(),
#ifdef __MPI
                                                                    this->pv.ncol_bands,
                                                                    this->pv.nrow,
#else
                                                                    nbands,
                                                                    nlocal,
#endif
                                                                    this->kv.ngk,
                                                                    true);
        }

        // Length of Hk_laststep and Sk_laststep per k-point: nlocal * nlocal for full matrices, nloc for local blocks
        const int len_HS_ik = use_tensor && use_lapack ? nlocal * nlocal : this->pv.nloc;
        const int len_HS_all = this->kv.get_nks() * len_HS_ik;

        // Allocate memory for Hk_laststep; use_tensor && use_lapack requires space for the full matrix
        if (this->Hk_laststep.NumElements() != len_HS_all)
        {
            this->Hk_laststep
                = ct::Tensor(ct::DataType::DT_COMPLEX_DOUBLE, ct_device_type_hs, ct::TensorShape({this->kv.get_nks(), len_HS_ik}));
            this->Hk_laststep.zero();
        }

        // Allocate memory for Sk_laststep; use_tensor && use_lapack requires space for the full matrix
        if (this->Sk_laststep.NumElements() != len_HS_all)
        {
            this->Sk_laststep
                = ct::Tensor(ct::DataType::DT_COMPLEX_DOUBLE, ct_device_type_hs, ct::TensorShape({this->kv.get_nks(), len_HS_ik}));
            this->Sk_laststep.zero();
        }

        // Put information into psi_laststep, Hk_laststep and Sk_laststep
        for (int ik = 0; ik < this->kv.get_nks(); ++ik)
        {
            this->psi->fix_k(ik);
            this->psi_laststep->fix_k(ik);

            // Copy data from psi to psi_laststep at k-point ik
            const int len_psi_ik = this->psi->get_nbands() * this->psi->get_nbasis();
            for (int index = 0; index < len_psi_ik; ++index)
            {
                psi_laststep[0].get_pointer()[index] = this->psi[0].get_pointer()[index];
            }

            // Get H and S matrices at k-point ik
            this->p_hamilt->updateHk(ik);
            hamilt::MatrixBlock<std::complex<double>> h_mat;
            hamilt::MatrixBlock<std::complex<double>> s_mat;
            static_cast<hamilt::Hamilt<std::complex<double>>*>(this->p_hamilt)->matrix(h_mat, s_mat);

            // Store H and S matrices to Hk_laststep and Sk_laststep
            if (use_tensor && use_lapack)
            {
#ifdef __MPI
                int myid = 0;
                int num_procs = 1;
                MPI_Comm_rank(MPI_COMM_WORLD, &myid);
                MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

                std::complex<double>* h_ptr = nullptr;
                std::complex<double>* s_ptr = nullptr;

                // Define containers for gathered data (only needed for multi-process)
                module_rt::Matrix_g<std::complex<double>> h_mat_g;
                module_rt::Matrix_g<std::complex<double>> s_mat_g;

                if (num_procs == 1)
                {
                    // Single process: directly point to local data without gather
                    h_ptr = h_mat.p;
                    s_ptr = s_mat.p;
                }
                else
                {
                    // Multiple processes: gather data to the root process (myid == 0) and point to the gathered data
                    module_rt::gatherMatrix(myid, 0, h_mat, h_mat_g);
                    module_rt::gatherMatrix(myid, 0, s_mat, s_mat_g);
                    if (myid == 0)
                    {
                        h_ptr = h_mat_g.p.get();
                        s_ptr = s_mat_g.p.get();
                    }
                }

                // Only the root process (myid == 0) performs the copy
                if (myid == 0 && h_ptr != nullptr && s_ptr != nullptr)
                {
                    BlasConnector::copy(len_HS_ik, h_ptr, 1, this->Hk_laststep.template data<std::complex<double>>() + ik * len_HS_ik, 1);
                    BlasConnector::copy(len_HS_ik, s_ptr, 1, this->Sk_laststep.template data<std::complex<double>>() + ik * len_HS_ik, 1);
                }
#endif
            }
            else
            {
                BlasConnector::copy(len_HS_ik, h_mat.p, 1, this->Hk_laststep.template data<std::complex<double>>() + ik * len_HS_ik, 1);
                BlasConnector::copy(len_HS_ik, s_mat.p, 1, this->Sk_laststep.template data<std::complex<double>>() + ik * len_HS_ik, 1);
            } // end use_tensor && use_lapack branch
        } // end ik
    } // save requested by the caller
}

template <typename TR, typename Device>
void ESolver_KS_LCAO_TDDFT<TR, Device>::after_scf(UnitCell& ucell, const int istep, const bool conv_esolver)
{
    ModuleBase::TITLE("ESolver_LCAO_TDDFT", "after_scf");
    ModuleBase::timer::start(this->classname, "after_scf");

    ESolver_KS_LCAO<std::complex<double>, TR>::after_scf(ucell, istep, conv_esolver);

    // Output energy for sub-loop (electronic step)
    std::cout << " Potential (Ry): " << std::setprecision(15) << this->pelec->f_en.etot << std::endl;

    // Match field output by skipping the replayed restart endpoint.
    const bool replay_restart = this->inp_->mdp.md_restart && !restart_done;
    if (!replay_restart)
    {
        auto* hamilt_lcao = dynamic_cast<hamilt::HamiltLCAO<std::complex<double>, TR>*>(this->p_hamilt);
        ModuleIO::ctrl_output_td<TR>(ucell,
                                     this->chr.rho_save,
                                     this->chr.rhopw,
                                     istep,
                                     this->psi,
                                     this->pelec,
                                     this->kv,
                                     this->two_center_bundle_.overlap_orb.get(),
                                     &this->pv,
                                     this->orb_,
                                     this->velocity_mat,
                                     this->gd,
                                     hamilt_lcao,
                                     this->RA,
                                     this->td_p,
                                     this->exx_nao,
                                     this->exx_info_);
    }

    ModuleBase::timer::end(this->classname, "after_scf");
}

template <typename TR, typename Device>
void ESolver_KS_LCAO_TDDFT<TR, Device>::weight_dm_rho(const UnitCell& ucell)
{
    if (this->inp_->ocp == 1)
    {
        elecstate::fixed_weights(this->inp_->ocp_kb,
                                 this->inp_->nbands,
                                 this->inp_->nelec,
                                 this->pelec->klist,
                                 this->pelec->wg,
                                 this->pelec->skip_weights);
    }

    // Calculate Eband energy
    elecstate::calEBand(this->pelec->ekb, this->pelec->wg, this->pelec->f_en);

    module_dm::dm_from_psi(&this->pv, this->pelec->wg, this->psi[0], *this->dmat.dm);
    if (this->inp_->td_stype == 2)
    {
        this->dmat.dm->cal_dmr_td(td_p->get_phase_hybrid(), td_p->endpoint_A(), -1);
    }
    else
    {
        this->dmat.dm->cal_dmr(-1);
    }

    // get the real-space charge density, mohan add 2025-10-24
    LCAO_domain::dm2rho(this->dmat.dm->get_dmr_vec(), this->inp_->nspin, &this->chr, this->inp_->nelec, ucell.omega, false);
}

template class ESolver_KS_LCAO_TDDFT<double, base_device::DEVICE_CPU>;
template class ESolver_KS_LCAO_TDDFT<std::complex<double>, base_device::DEVICE_CPU>;
#if ((defined __CUDA) /* || (defined __ROCM) */)
template class ESolver_KS_LCAO_TDDFT<double, base_device::DEVICE_GPU>;
template class ESolver_KS_LCAO_TDDFT<std::complex<double>, base_device::DEVICE_GPU>;
#endif

} // namespace ModuleESolver
