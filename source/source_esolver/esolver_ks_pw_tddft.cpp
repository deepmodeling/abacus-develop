#include "source_esolver/esolver_ks_pw_tddft.h"

#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "source_base/parallel_comm.h"
#include "source_base/parallel_device.h"
#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include "source_estate/elecstate_pw.h"
#include "source_estate/elecstate_tools.h"
#include "source_estate/module_charge/chg_symm.h"
#include "source_hamilt/hamilt_hs_adapter.h"
#include "source_hamilt/module_xc/xc_functional.h"
#include "source_io/module_efield/td_efield_io.h"
#include "source_io/module_efield/td_vector_pot_io.h"
#include "source_io/module_energy/td_eij_io.h"
#include "source_io/module_output/output_log.h"
#include "source_io/module_parameter/parameter.h"
#include "source_pw/module_pwdft/hamilt_pw.h"
#include "source_pw/module_pwdft/td_pw.h"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <limits>

namespace ModuleESolver
{

template <typename T, typename Device>
ESolver_KS_PW_TDDFT<T, Device>::ESolver_KS_PW_TDDFT() : log_(GlobalV::ofs_running)
{
    this->classname = "ESolver_KS_PW_TDDFT";
    this->basisname = "PW";
}

template <typename T, typename Device>
void ESolver_KS_PW_TDDFT<T, Device>::before_all_runners(BaseCell& basecell, const Input_para& inp)
{
    basecell.require_kind(BaseCell::Kind::unitcell, __FUNCTION__);
    UnitCell& ucell = static_cast<UnitCell&>(basecell);
    const std::string& output_dir = PARAM.globalv.global_out_dir;

    this->td_field_manager_ = elecstate::create_td_field_manager(inp);
    if (inp.init_vecpot)
    {
        // An empty directory selects the calculation working directory, not OUT.<suffix>.
        const std::vector<ModuleBase::Vector3<double>> samples = ModuleIO::read_td_vector_pot("");
        this->td_field_manager_->set_A_samples(samples);
    }

    ESolver_KS_PW<T, Device>::before_all_runners(ucell, inp);
    // The parent resolves the XC functional, including pseudopotential defaults.
    pw::check_td_input(inp, ucell, XC_Functional::get_ked_flag(), this->general_exx_info_.cal_exx);
    this->pelec->pot->set_td_field_manager(this->td_field_manager_);
#ifdef __MPI
    const hsolver::diag_comm_info comm(this->pw_wfc->pool_world, this->pw_wfc->poolrank, this->pw_wfc->poolnproc);
#else
    const hsolver::diag_comm_info comm(0, 1);
#endif
    hsolver::PWLinearOptions options;
    options.orthonormal = hsolver::parse_orth_method(inp.td_orthonormal);
    options.out_stat = inp.td_out_stat;
    options.global_k_indices = this->kv.ik2iktot;
    options.linear.method = hsolver::parse_linear_method(inp.lin_solver);
    options.linear.tolerance = inp.lin_thr;
    options.linear.max_iterations = inp.lin_maxiter;
    options.linear.restart = inp.lin_gmres_restart;
    options.linear.reconstruct = inp.lin_reconstruct;
    options.preconditioner = hsolver::parse_pw_precond(inp.lin_precond);
    options.cn_init = inp.td_cn_init;
    options.kinetic_enabled = inp.t_in_h;
    this->td_solver_.reset(new hsolver::HSolverPWTDDFT<T, Device>(*this->pw_wfc, options, comm, log_));
    this->history_.prepare(*this->pelec->pot, XC_Functional::get_ked_flag());
    // Preserve existing field history until input validation and initialization succeed.
    if (GlobalV::MY_RANK == 0)
    {
        if (inp.out_efield)
        {
            ModuleIO::prepare_td_field_output(output_dir, *this->td_field_manager_, false);
        }
        if (inp.out_vecpot)
        {
            if (inp.init_vecpot)
            {
                ModuleIO::check_td_vector_pot_paths("", output_dir);
            }
            ModuleIO::prepare_td_vector_pot_output(output_dir, false);
        }
    }
}

template <typename T, typename Device>
void ESolver_KS_PW_TDDFT<T, Device>::runner(BaseCell& basecell, const int istep)
{
    const int steps_per_md = this->inp_->estep_per_md;
    int first_step = 0;
    int last_step = 0;
    if (this->inp_->mdp.md_nstep == 1)
    {
        // Include the initial state before the fixed-ion propagation sequence.
        last_step = steps_per_md;
    }
    else if (istep > 0)
    {
        first_step = (istep - 1) * steps_per_md + 1;
        last_step = istep * steps_per_md;
    }
    for (int electronic_step = first_step; electronic_step <= last_step; ++electronic_step)
    {
        ESolver_KS::runner(basecell, electronic_step);
    }
}

template <typename T, typename Device>
void ESolver_KS_PW_TDDFT<T, Device>::print_step(const int istep) const
{
    const int evolution_step = istep + 1;
    std::cout << " -------------------------------------------" << std::endl;
    std::cout << " STEP OF ELECTRON EVOLVE : " << evolution_step << std::endl;
    std::cout << " -------------------------------------------" << std::endl;
}

template <typename T, typename Device>
void ESolver_KS_PW_TDDFT<T, Device>::before_scf(UnitCell& ucell, const int istep)
{
    this->print_step(istep);
    this->prepare_td_step(istep);
    this->td_solver_->reset_orth_stats();
    const bool basis_updated = ucell.cell_parameter_updated;
    ESolver_KS_PW<T, Device>::before_scf(ucell, istep);
    if (basis_updated)
    {
        this->td_solver_->invalidate_basis();
    }
    this->history_.prepare(*this->pelec->pot, XC_Functional::get_ked_flag());
}

template <typename T, typename Device>
void ESolver_KS_PW_TDDFT<T, Device>::prepare_init_hamilt(const UnitCell& ucell)
{
    if (this->td_field_manager_->gauge() == 1)
    {
        // The updated basis and endpoint A determine the required nonlocal projector range.
        const double q_unshifted = pw::td_momentum_bound(*this->pw_wfc, ucell.tpiba);
        pw::ensure_td_vnl(ucell, q_unshifted, this->td_field_manager_->A_right_ha(), this->td_field_manager_->A_prop_ha(), &this->ppcell);
        elecstate::Potential& potential = *this->pelec->pot;
        const Real* veff = potential.template get_veff_smooth_data<Real>();
        const Real* vofk = potential.template get_vofk_smooth_data<Real>();
        auto* hamiltonian = static_cast<hamilt::HamiltPW<T, Device>*>(this->p_hamilt);
        // A file can supply nonzero A(0), which must enter H before the initial wavefunctions are prepared.
        hamiltonian->bind_td_state(veff, vofk, this->td_field_manager_->A_right_ha());
    }
}

template <typename T, typename Device>
void ESolver_KS_PW_TDDFT<T, Device>::prepare_td_step(const int istep)
{
    ModuleBase::timer::start("ESolver_KS_PW_TDDFT", "prepare_td_step");
    if (this->td_field_manager_->current_step() == istep)
    {
        ModuleBase::timer::end("ESolver_KS_PW_TDDFT", "prepare_td_step");
        return;
    }
    if (this->td_field_manager_->gauge() == 0 || istep == 0)
    {
        this->td_field_manager_->prepare_sample(istep);
    }
    else
    {
        const int interval_left = istep - 1;
        this->td_field_manager_->prepare_interval(interval_left, istep);
    }
    if (this->td_field_manager_->gauge() == 1)
    {
        const ModuleBase::Vector3<double>& A_right_ha = this->td_field_manager_->A_right_ha();
        ModuleBase::GlobalFunc::OUT(log_, "Cartesian vector potential Ax(t)", A_right_ha[0]);
        ModuleBase::GlobalFunc::OUT(log_, "Cartesian vector potential Ay(t)", A_right_ha[1]);
        ModuleBase::GlobalFunc::OUT(log_, "Cartesian vector potential Az(t)", A_right_ha[2]);
    }
    if (GlobalV::MY_RANK == 0)
    {
        if (this->inp_->out_vecpot)
        {
            ModuleIO::write_td_vector_pot(PARAM.globalv.global_out_dir, istep, this->td_field_manager_->A_right_ha());
        }
        if (this->inp_->out_efield)
        {
            ModuleIO::write_td_field_values(*this->td_field_manager_, PARAM.globalv.global_out_dir);
        }
    }
    ModuleBase::timer::end("ESolver_KS_PW_TDDFT", "prepare_td_step");
}

template <typename T, typename Device>
std::string ESolver_KS_PW_TDDFT<T, Device>::diag_policy(const int istep) const
{
    if (istep == 0)
    {
        return "ksdft";
    }
    return ESolver_KS::diag_policy(istep);
}

template <typename T, typename Device>
void ESolver_KS_PW_TDDFT<T, Device>::iter_init(UnitCell& ucell, const int istep, const int iter)
{
    int ionic_step = 0;
    if (this->inp_->mdp.md_nstep > 1 && istep > 0)
    {
        ionic_step = (istep - 1) / this->inp_->estep_per_md + 1;
    }
    ModuleIO::write_head_td(GlobalV::ofs_running, ionic_step, istep, iter, this->basisname);
    ESolver_KS_PW<T, Device>::iter_init(ucell, istep, iter);
}

template <typename T, typename Device>
ESolver_KS::DensityStage ESolver_KS_PW_TDDFT<T, Device>::density_stage(const int istep, const int iter) const
{
    if (istep > 0 && iter == 1)
    {
        return ESolver_KS::DensityStage::predictor;
    }
    return ESolver_KS::DensityStage::standard;
}

template <typename T, typename Device>
void ESolver_KS_PW_TDDFT<T, Device>::hamilt2rho_single(UnitCell& ucell, const int istep, const int iter, const double ethr)
{
    ModuleBase::timer::start("ESolver_KS_PW_TDDFT", "hamilt2rho_single");
    if (istep == 0)
    {
        ESolver_KS_PW<T, Device>::hamilt2rho_single(ucell, istep, iter, ethr);
        if (this->inp_->td_orthonormal != "none")
        {
            psi::Psi<T, Device>* current = this->stp.template get_psi_t<T, Device>();
            int changed = this->td_solver_->correct_initial(current, iter);
            // Occupations and density construction contain collectives across k-point pools.
            Parallel_Reduce::reduce_all(changed);
            if (changed > 0)
            {
                hamilt::Hamilt<T, Device>* hamiltonian = static_cast<hamilt::Hamilt<T, Device>*>(this->p_hamilt);
                hamilt::HamiltHSOperator<T, Device> op(hamiltonian, this->pw_wfc);
                this->td_solver_->cal_band_energy(op, *current, &this->pelec->ekb);
                elecstate::calculate_weights(this->pelec->ekb,
                                             this->pelec->wg,
                                             this->pelec->klist,
                                             this->pelec->eferm,
                                             this->pelec->f_en,
                                             this->pelec->nelec_spin,
                                             this->inp_->nbands,
                                             this->pelec->skip_weights);
                elecstate::calEBand(this->pelec->ekb, this->pelec->wg, this->pelec->f_en);
                elecstate::ElecStatePW<T, Device>* estate = static_cast<elecstate::ElecStatePW<T, Device>*>(this->pelec);
                estate->psiToRho(*current);
                module_charge::symmetrize_rho(this->inp_->nspin, this->chr, this->pw_rhod, ucell.symm);
            }
        }
        ModuleBase::timer::end("ESolver_KS_PW_TDDFT", "hamilt2rho_single");
        return;
    }

    ModuleBase::TITLE("ESolver_KS_PW_TDDFT", "hamilt2rho_single (RT-TDDFT)");
    psi::Psi<T, Device>* current = this->stp.template get_psi_t<T, Device>();
    elecstate::Potential& potential = *this->pelec->pot;
    hamilt::HamiltPW<T, Device>* hamiltonian = static_cast<hamilt::HamiltPW<T, Device>*>(this->p_hamilt);
    hamilt::HamiltHSOperator<T, Device> op(hamiltonian, this->pw_wfc);

    // Predict with the current potential, then correct with the midpoint.
    const Real* endpoint_veff = potential.template get_veff_smooth_data<Real>();
    const Real* endpoint_vofk = potential.template get_vofk_smooth_data<Real>();
    const ESolver_KS::DensityStage stage = this->density_stage(istep, iter);
    const bool corrector = stage != ESolver_KS::DensityStage::predictor;
    const pw::TDPotentialView<Real> propagation = this->history_.prop_pot(potential, XC_Functional::get_ked_flag(), corrector);
    ModuleBase::Vector3<double> momentum_shift(0.0, 0.0, 0.0);
    if (this->td_field_manager_->gauge() == 1)
    {
        momentum_shift = this->td_field_manager_->A_prop_ha();
    }
    hamiltonian->bind_td_state(propagation.veff, propagation.vofk, momentum_shift);
    // The solver sequence follows electronic steps across MD boundaries.
    this->td_solver_->solve(op,
                            this->history_.previous(),
                            current,
                            this->td_field_manager_->dt_ha(),
                            momentum_shift,
                            istep,
                            iter,
                            this->inp_->out_level == "ie");

    // Restore the endpoint Hamiltonian before evaluating density and energy.
    hamiltonian->bind_td_state(endpoint_veff, endpoint_vofk, this->td_field_manager_->A_right_ha());
    elecstate::ElecStatePW<T, Device>* estate = static_cast<elecstate::ElecStatePW<T, Device>*>(this->pelec);
    // Propagation keeps the occupations of the initial state fixed.
    estate->psiToRho(*current);
    module_charge::symmetrize_rho(this->inp_->nspin, this->chr, this->pw_rhod, ucell.symm);

    // Pair the band sum with KS double-counting terms using this same potential.
    // The later KS potential update belongs to the next SCF iteration.
    this->td_solver_->cal_band_energy(op, *current, &this->pelec->ekb);
    elecstate::calEBand(this->pelec->ekb, this->pelec->wg, this->pelec->f_en);
    ModuleBase::timer::end("ESolver_KS_PW_TDDFT", "hamilt2rho_single");
}

template <typename T, typename Device>
void ESolver_KS_PW_TDDFT<T, Device>::iter_finish(UnitCell& ucell, const int istep, int& iter, bool& conv_esolver)
{
    ESolver_KS_PW<T, Device>::iter_finish(ucell, istep, iter, conv_esolver);
    if (istep > 0 && !conv_esolver && iter >= this->inp_->scf_nmax)
    {
        const int evolution_step = istep + 1;
        const std::string message = "The time-dependent SCF step did not converge. Electronic evolution step: "
                                    + std::to_string(evolution_step) + ", SCF iteration: " + std::to_string(iter) + ".";
        ModuleBase::WARNING_QUIT("ESolver_KS_PW_TDDFT", message);
    }
}

template <typename T, typename Device>
void ESolver_KS_PW_TDDFT<T, Device>::after_scf(UnitCell& ucell, const int istep, const bool conv_esolver)
{
    if (!conv_esolver)
    {
        const int evolution_step = istep + 1;
        const std::string message
            = "Cannot propagate an unconverged electronic state. Electronic evolution step: " + std::to_string(evolution_step) + ".";
        ModuleBase::WARNING_QUIT("ESolver_KS_PW_TDDFT", message);
    }
    if (istep == 0 && this->inp_->td_orthonormal == "none")
    {
        const psi::Psi<T, Device>* current = this->stp.template get_psi_t<T, Device>();
        this->td_solver_->check_initial(*current, this->niter);
    }
    if (this->inp_->td_print_eij >= 0.0)
    {
        print_eij(istep);
    }
    ESolver_KS_PW<T, Device>::after_scf(ucell, istep, conv_esolver);
    const std::streamsize original_precision = std::cout.precision();
    std::cout << " Potential (Ry): " << std::setprecision(15) << this->pelec->f_en.etot << std::endl;
    std::cout.precision(original_precision);
    if (istep >= 0)
    {
        psi::Psi<T, Device>* current = this->stp.template get_psi_t<T, Device>();
        this->history_.save(*current, *this->pelec->pot, XC_Functional::get_ked_flag());

        if (this->inp_->out_current == 1)
        {
            const ModuleBase::Vector3<double>& A_right_ha = this->td_field_manager_->A_right_ha();
            this->current_output_.write(istep,
                                        ucell,
                                        this->pw_wfc,
                                        current,
                                        this->pelec,
                                        this->kv,
                                        &this->ppcell,
                                        this->td_field_manager_->gauge(),
                                        A_right_ha,
                                        this->inp_->out_current_k,
                                        PARAM.globalv.global_out_dir,
                                        GlobalV::MY_RANK);
        }
    }
    if (this->inp_->td_out_stat)
    {
        report_orth(ucell, istep);
    }
}

template <typename T, typename Device>
void ESolver_KS_PW_TDDFT<T, Device>::print_eij(const int istep)
{
    ModuleBase::timer::start("ESolver_KS_PW_TDDFT", "print_eij");
    psi::Psi<T, Device>* current = this->stp.template get_psi_t<T, Device>();
    const int original_k = current->get_current_k();
    const int nband = current->get_nbands();
    const size_t elements = static_cast<size_t>(nband) * nband;
    const int nspin = this->inp_->nspin;
    const int nk = this->kv.get_nkstot() / (nspin == 2 ? 2 : 1);
    const Parallel_Kpoints& pools = this->kv.para_k;
    const bool main_rank = pools.my_pool == 0 && pools.rank_in_pool == 0;
    auto* hamiltonian = static_cast<hamilt::HamiltPW<T, Device>*>(this->p_hamilt);
    const Real* veff = this->pelec->pot->template get_veff_smooth_data<Real>();
    const Real* vofk = this->pelec->pot->template get_vofk_smooth_data<Real>();
    // cal_converged() updates the potential from the final density, after the last band-energy evaluation.
    hamiltonian->bind_td_state(veff, vofk, this->td_field_manager_->A_right_ha());
    hamilt::HamiltHSOperator<T, Device> op(hamiltonian, this->pw_wfc);
    for (int global_k = 0; global_k < this->kv.get_nkstot(); ++global_k)
    {
        const int physical_k = global_k % nk;
        const int spin = global_k / nk;
        int owner = 0;
        int local_k = global_k;
#ifdef __MPI
        owner = pools.whichpool[physical_k];
        local_k = physical_k - pools.startk_pool[owner] + spin * pools.nks_pool[owner];
#endif
        std::vector<std::complex<double>> eij;
        if (pools.my_pool == owner)
        {
            eij = this->td_solver_->cal_eij(op, *current, local_k);
        }
#ifdef __MPI
        // Each pool already has the full matrix; transfer only one copy to the main log.
        if (owner != 0 && pools.rank_in_pool == 0 && (main_rank || pools.my_pool == owner))
        {
            if (main_rank)
            {
                eij.resize(elements);
            }
            MPI_Comm world = MPI_COMM_WORLD;
            const int source = pools.get_startpro_pool(owner);
            const size_t max_count = std::numeric_limits<int>::max();
            for (size_t offset = 0; offset < elements; offset += max_count)
            {
                const int count = static_cast<int>(std::min(max_count, elements - offset));
                if (main_rank)
                {
                    Parallel_Common::recv_data(eij.data() + offset, count, source, 0, world, MPI_STATUS_IGNORE);
                }
                else
                {
                    Parallel_Common::send_data(eij.data() + offset, count, 0, 0, world);
                }
            }
        }
#endif
        if (main_rank)
        {
            const ModuleIO::EijOutputInfo info{this->inp_->td_print_eij, istep + 1, this->niter, physical_k + 1, nspin, spin + 1};
            ModuleIO::print_eij(eij.data(), nband, nband, info, log_);
        }
    }
    current->fix_k(original_k);
    op.update_k(original_k);
    ModuleBase::timer::end("ESolver_KS_PW_TDDFT", "print_eij");
}

template <typename T, typename Device>
void ESolver_KS_PW_TDDFT<T, Device>::report_orth(const UnitCell& ucell, int istep)
{
    ModuleBase::timer::start("ESolver_KS_PW_TDDFT", "report_orth");
    double electrons = this->td_solver_->wave_electrons(this->pelec->wg);
    // Norms are replicated within a pool; only its root contributes to the total.
    if (this->pw_wfc->poolrank != 0)
    {
        electrons = 0.0;
    }
    Parallel_Reduce::reduce_all(electrons);
    double rho_electrons = 0.0;
    const int density_spins = this->chr.nspin == 2 ? 2 : 1;
    for (int spin = 0; spin < density_spins; ++spin)
    {
        for (int ir = 0; ir < this->chr.nrxx; ++ir)
        {
            rho_electrons += this->chr.rho[spin][ir];
        }
    }
    rho_electrons *= ucell.omega / this->chr.rhopw->nxyz;
    Parallel_Reduce::reduce_pool(rho_electrons);
    if (istep == 0)
    {
        initial_wave_electrons_ = electrons;
        initial_rho_electrons_ = rho_electrons;
    }
    const hsolver::TDOrthStats& stats = this->td_solver_->orth_stats();
    const bool full_gram = this->inp_->td_orthonormal != "none";
    double maxima[] = {stats.before, stats.after};
    if (full_gram)
    {
        Parallel_Reduce::reduce_max(maxima, 2);
    }
    if (!log_.good())
    {
        ModuleBase::timer::end("ESolver_KS_PW_TDDFT", "report_orth");
        return;
    }
    const char* stage = istep == 0 ? "initial state" : "propagation";
    const int evolution_step = istep + 1;
    const double delta_wave = electrons - initial_wave_electrons_;
    const double delta_rho = rho_electrons - initial_rho_electrons_;
    const std::ios::fmtflags saved_flags = log_.flags();
    const std::streamsize saved_precision = log_.precision();
    log_ << " PW RT-TDDFT conservation: evolution step " << evolution_step << ", " << stage << ", " << this->inp_->td_orthonormal << '\n';
    log_ << std::right << std::fixed << std::setprecision(10);
    ModuleBase::GlobalFunc::OUT(log_, "Npsi", electrons);
    ModuleBase::GlobalFunc::OUT(log_, "Nrho", rho_electrons);
    log_ << std::scientific << std::setprecision(6);
    ModuleBase::GlobalFunc::OUT(log_, "dNpsi", delta_wave);
    ModuleBase::GlobalFunc::OUT(log_, "dNrho", delta_rho);
    if (full_gram)
    {
        if (istep == 0)
        {
            ModuleBase::GlobalFunc::OUT(log_, "initial_orth_error", maxima[1]);
        }
        else
        {
            ModuleBase::GlobalFunc::OUT(log_, "orth_before_step_max", maxima[0]);
            ModuleBase::GlobalFunc::OUT(log_, "orth_after_step_max", maxima[1]);
        }
    }
    log_.flags(saved_flags);
    log_.precision(saved_precision);
    ModuleBase::timer::end("ESolver_KS_PW_TDDFT", "report_orth");
}

template class ESolver_KS_PW_TDDFT<std::complex<float>, base_device::DEVICE_CPU>;
template class ESolver_KS_PW_TDDFT<std::complex<double>, base_device::DEVICE_CPU>;
#if defined(__CUDA) || defined(__ROCM)
template class ESolver_KS_PW_TDDFT<std::complex<float>, base_device::DEVICE_GPU>;
template class ESolver_KS_PW_TDDFT<std::complex<double>, base_device::DEVICE_GPU>;
#endif

} // namespace ModuleESolver
