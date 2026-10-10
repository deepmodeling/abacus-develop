// LR execution flow and final spectrum/density output.
#include "esolver_lr_lcao_tddft.h"
#include "source_lcao/module_lr/utils/lr_io.h"
#include "source_lcao/module_lr/utils/lr_util.h"
#include "source_lcao/module_lr/hamilt_ulr.hpp"
#include "source_lcao/module_lr/hsolver_lrtd.hpp"
#include "source_lcao/module_lr/lr_spectrum.h"
#include "source_lcao/module_lr/lr_density.hpp"
#include "source_lcao/module_lr/ri_benchmark/ri_benchmark.h"
#include "source_lcao/module_lr/operator_casida/operator_lr_diag.h"
#include "source_lcao/module_lr/zeq_solver.h"
#include "source_hamilt/module_xc/xc_functional.h"
#include <algorithm>
#include <iomanip>

using namespace LR;

template <typename T, typename TR>
ModuleESolver::ESolver_LR<T, TR>::ESolver_LR(const Input_para& inp,
                                            const std::string& in_dir,
                                            const std::string& out_dir)
    : in_dir(in_dir), out_dir(out_dir)
{
#ifdef __EXX
    init_exx_info(this->exx_info, inp);
#endif
}

template <typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::runner(BaseCell& basecell, const int istep)
{
    basecell.require_kind(BaseCell::Kind::unitcell, __FUNCTION__);
    UnitCell& ucell = static_cast<UnitCell&>(basecell);

    ModuleBase::TITLE("ESolver_LR", "runner");
    ModuleBase::timer::start("ESolver_LR", "runner");

    if (this->ks_)
    {
        // `init_pot_groundstate` leaves the global XC type on the LR kernel, so put it back before
        // the ground-state SCF. `ESolver_KS::runner` is re-entrant: its `before_scf` rebuilds the
        // neighbour lists, the grid tables and the Hamiltonian for the geometry of this step.
        XC_Functional::set_xc_type(ucell.atoms[0].ncpp.xc_func);
        this->ks_->runner(ucell, istep);
        this->etot_gs_ = this->ks_->cal_energy();
        if (!this->ks_initialized_)
        {
            this->initialize_from_ks_(ucell, *this->inp_);
            this->setup_relax_target_();   // needs the final dimensions, i.e. `openshell`
        }
        else { this->refresh_from_ks_(ucell); }
    }

    //allocate 2-particle state and setup 2d division
    this->setup_eigenvectors_X();
    this->pelec->ekb.create(nspin, this->nstates);

    // set once here rather than inside the solver branches: the open-shell `spectrum` branch used
    // to leave it empty, and both `after_all_runners` and the gradient index it
    this->spin_types = this->openshell ? std::vector<std::string>({ "updown" })
                                       : std::vector<std::string>({ "singlet", "triplet" });

    // a relaxation writes these once per ionic step; without the suffix every step would overwrite
    // the last, and the trajectory would be impossible to inspect afterwards
    const std::string step_suffix = this->excited_relax_ ? "_step" + std::to_string(istep) : "";
    auto efile_out = [&](const std::string& label)->std::string {return this->out_dir + "Excitation_Energy_" + label + step_suffix + ".dat";};
    auto vfile_out = [&](const std::string& label)->std::string {return this->out_dir + "Excitation_Amplitude_" + label + step_suffix + "_" + std::to_string(this->my_rank_+1) + ".dat";};
    if (this->inp_->lr_solver == "elpa")
    {
        ModuleBase::WARNING_QUIT("ESolver_LR", "ESolver_LR doesn't support elpa now.");

    }
    else if (this->inp_->lr_solver != "spectrum")
    {
        auto write_states = [&](const std::string& label, const Real<T>* e, const T* v, const int& dim, const int& nst, const int& prec = 8)->void
            {
                if (this->my_rank_ == 0) { assert(nst == LR_Util::write_value(efile_out(label), prec, e, nst)); }
                assert(nst * dim == LR_Util::write_value(vfile_out(label), prec, v, nst, dim));
            };
        std::vector<double> precondition(this->inp_->lr_solver == "lapack" ? 0 : nloc_per_state, 1.0);
        // allocate and initialize A matrix and density matrix
        if (openshell)
        {
            for (int is : {0, 1})
            {
                if (this->inp_->lr_solver != "lapack") {
                    const int offset_is = is * this->paraX_[0].get_local_size();
                    OperatorLRDiag<double> pre_op(this->eig_ks.c + is * nk * (nocc[0] + nvirt[0]), this->paraX_[is], this->nk, this->nocc[is], this->nvirt[is]);
                    pre_op.act(1, offset_is, 1, precondition.data() + offset_is, precondition.data() + offset_is);
                }
            }
            std::cout << "Solving spin-conserving excitation for open-shell system." << std::endl;
            HamiltULR<T> hulr(xc_kernel,
                              nspin,
                              this->nbasis,
                              this->nocc,
                              this->nvirt,
                              *this->ucell_,
                              orb_cutoff_,
                              this->gd(),
                              *this->psi_ks,
                              this->eig_ks,
#ifdef __EXX
                              this->exx_lri,
                              this->exx_info.info_global.hybrid_alpha,
#endif
                              this->pot,
                              this->kv,
                              this->paraX_,
                              this->paraC_,
                              this->paraMat_);
            LR::HSolver::solve(hulr, this->X[0].template data<T>(), nloc_per_state, nstates,
                               this->nk, this->nocc, this->nvirt, this->paraX_,
                               this->pelec->ekb.c, this->inp_->lr_solver,
                               this->inp_->lr_thr, precondition);
            if (this->inp_->out_wfc_lr) { write_states("openshell", this->pelec->ekb.c, this->X[0].template data<T>(), nloc_per_state, nstates); }
        }
        else
        {
            if (this->inp_->lr_solver != "lapack") {
                OperatorLRDiag<double> pre_op(this->eig_ks.c, this->paraX_[0], this->nk, this->nocc[0], this->nvirt[0]);
                pre_op.act(1, nloc_per_state, 1, precondition.data(), precondition.data()); 
            }
            const std::vector<std::string>& spin_types = this->spin_types;
            for (int is = 0;is < nspin;++is)
            {
                std::cout << " Calculating " << spin_types[is] << " excitations" << std::endl;
                HamiltLR<T> hlr(xc_kernel,
                                nspin,
                                this->nbasis,
                                this->nocc,
                                this->nvirt,
                                *this->ucell_,
                                orb_cutoff_,
                                this->gd(),
                                *this->psi_ks,
                                this->eig_ks,
#ifdef __EXX
                                this->exx_lri,
                                this->exx_info.info_global.hybrid_alpha,
#endif
                                this->pot[is],
                                this->kv,
                                this->paraX_,
                                this->paraC_,
                                this->paraMat_,
                                spin_types[is],
                                this->in_dir,
                                this->out_dir,
                                this->inp_->ri_hartree_benchmark,
                                (this->inp_->ri_hartree_benchmark == "aims" ? this->inp_->aims_nbasis : std::vector<int>({})));
                LR::HSolver::solve(hlr, this->X[is].template data<T>(), nloc_per_state, nstates,
                                   this->nk, this->nocc, this->nvirt, this->paraX_,
                                   this->pelec->ekb.c + is * nstates,
                                   this->inp_->lr_solver,
                                   this->inp_->lr_thr,
                                   precondition);
                if (this->inp_->out_wfc_lr) { write_states(spin_types[is], this->pelec->ekb.c + is * nstates, this->X[is].template data<T>(), nloc_per_state, nstates); }
            }
        }
    }
    else    // lr_solver == "spectrum", read the eigenvalues
    {
        auto efile_in = [&](const std::string& label)->std::string {return this->in_dir + "Excitation_Energy_" + label + ".dat";};
        auto vfile_in = [&](const std::string& label)->std::string {return this->in_dir + "Excitation_Amplitude_" + label + "_" + std::to_string(this->my_rank_+1) + ".dat";};
    
        auto read_states = [&](const std::string& label, Real<T>* e, T* v, const int& dim, const int& nst)->void
            {
                if (this->my_rank_ == 0) {
                    assert(nst == LR_Util::read_value(efile_in(label), e, nst));
                    std::cout <<"Rank "<< this->my_rank_ << ": finish reading " << efile_in(label) << std::endl;
                }
#ifdef __MPI
// in velocity gauge, the eigenvalues may be used to calculate the transition dipole, so we'd better broadcast them
                MPI_Bcast(e, nst, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
                assert(nst * dim == LR_Util::read_value(vfile_in(label), v, nst, dim));
                std::cout <<"Rank "<< this->my_rank_ << ": finish reading " << vfile_in(label) << std::endl;
            };
        std::cout << "reading the excitation states from file: \n";
        if (openshell)
        {
            read_states("openshell", this->pelec->ekb.c, this->X[0].template data<T>(), nloc_per_state, nstates);
        }
        else
        {
            const std::vector<std::string>& spin_types = this->spin_types;
            for (int is = 0;is < nspin;++is) { read_states(spin_types[is], this->pelec->ekb.c + is * nstates, this->X[is].template data<T>(), nloc_per_state, nstates); }
        }
    }
    if (this->excited_relax_)
    {
        // Re-select which root to follow BEFORE taking the gradient, so the force belongs to
        // the same diabatic state as the previous step's.
        this->follow_target_state_(this->ofs_running_);

        // The LR terms only carry the Omega part of the force; the ground-state part is separate
        // and comes straight from the KS solver.
        this->ks_->cal_force(ucell, this->force_gs_);
        this->lr_force_ = this->cal_lr_force_relax_(this->ofs_running_);

        // One line per ionic step with the two halves of the energy and of the gradient. Without
        // it the relaxation only reports a force, and whether E_gs + Omega actually goes down --
        // the thing being minimised -- cannot be read off the log at all.
        const double omega = this->target_omega_();
        auto max_abs = [](const ModuleBase::matrix& m) -> double
            { double v = 0.0; for (int i = 0;i < m.nr * m.nc;++i) { v = std::max(v, std::abs(m.c[i])); } return v; };
        this->ofs_running_ << std::setprecision(8) << std::fixed
            << " EXCITED-STATE RELAX step " << istep
            << ": E_gs = " << this->etot_gs_ * ModuleBase::Ry_to_eV
            << " eV, Omega = " << omega * ModuleBase::Ry_to_eV
            << " eV, E_exc = " << (this->etot_gs_ + omega) * ModuleBase::Ry_to_eV << " eV"
            << std::setprecision(6)
            << " | |F_gs|max = " << max_abs(this->force_gs_) * ModuleBase::Ry_to_eV / ModuleBase::BOHR_TO_A
            << ", |F_Omega|max = " << max_abs(this->lr_force_) * ModuleBase::Ry_to_eV / ModuleBase::BOHR_TO_A
            << " eV/Angstrom" << std::defaultfloat << std::endl;
    }

    ModuleBase::timer::end("ESolver_LR", "runner");
    return;
}

template <typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::after_all_runners(BaseCell& basecell)
{
    basecell.require_kind(BaseCell::Kind::unitcell, __FUNCTION__);
    UnitCell& ucell = static_cast<UnitCell&>(basecell);

    ModuleBase::TITLE("ESolver_LR", "after_all_runners");
    if (this->inp_->ri_hartree_benchmark != "none") { return; } //no need to calculate the spectrum in the benchmark routine

    // cal electron-hole density
    if (this->inp_->out_chg[0])
    {
        LR_Density<T> lr_density(*this->ucell_, kv, gd(), *psi_ks, orb_cutoff_, pgrid(),
            nspin, nocc, nvirt, nbasis,
            paraX_, paraC_, paraMat_, openshell);

        if (openshell)
        {
            for (int is = 0;is < this->nspin;++is)
            {
                const T* const amplitudes = this->X[0].template data<T>();
                lr_density.output_eh_density_all_states(amplitudes, is, nstates, this->ofs_running_);
            }
        }
        else
        {
            for (int is = 0;is < this->X.size();++is)
            {
                const T* const amplitudes = this->X[is].template data<T>();
                lr_density.output_eh_density_all_states(amplitudes, is, nstates, this->ofs_running_);
            }
        }
    }

    //cal spectrum
    if (LR_Util::tolower(this->inp_->abs_gauge) == "velocity" )
    {
        const int nspin_tmp = this->inp_->nspin == 2 ? 2 : 1;
        this->velocity_mo = LR_Util::cal_velocity_mo(*this->ucell_, this->gd(), this->tcb(), 
                                                    this->paraMat_, this->paraC_, this->kv, *this->psi_ks, 
                                                    this->nk, nspin_tmp, this->nbasis, this->nocc, this->nvirt);
    }
    std::vector<double> freq(100);
    std::vector<double> abs_wavelen_range({ 20, 200 });//default range
    if (this->inp_->abs_wavelen_range.size() >= 2 && std::abs(this->inp_->abs_wavelen_range[1] - this->inp_->abs_wavelen_range[0]) > 0.02)
    {
        abs_wavelen_range = this->inp_->abs_wavelen_range;
    }
    double lambda_diff = std::abs(abs_wavelen_range[1] - abs_wavelen_range[0]);
    double lambda_min = std::min(abs_wavelen_range[1], abs_wavelen_range[0]);
    for (int i = 0;i < freq.size();++i) { freq[i] = 91.126664 / (lambda_min + 0.01 * static_cast<double>(i + 1) * lambda_diff); }
    // auto spin_types = (nspin == 2 && !openshell) ? std::vector<std::string>({ "singlet", "triplet" }) : std::vector<std::string>({ "updown" });
    // for (int is = 0;is < this->X.size() - 1;++is)
    for (int is = 0;is < this->X.size();++is)
    {
        LR_Spectrum<T> spectrum(nspin, this->nbasis, this->nocc, this->nvirt, *this->pw_rho, *this->psi_ks,
            *this->ucell_, this->kv, this->gd(), this->orb_cutoff_, this->tcb(),
            this->paraX_, this->paraC_, this->paraMat_,
            &this->pelec->ekb.c[is * nstates], this->eig_ks.c, this->X[is].template data<T>(), nstates, openshell,
            LR_Util::tolower(this->inp_->abs_gauge), this->my_rank_, this->out_dir);
        if (LR_Util::tolower(this->inp_->abs_gauge) == "velocity" ) {spectrum.set_vmo(this->velocity_mo.data());}
        spectrum.cal_spectrum();
        spectrum.transition_analysis(spin_types[is]+"_tda");
        if (spin_types[is] != "triplet")        // triplets has no transition dipole and no contribution to the spectrum
        {
            spectrum.optical_absorption_method1(freq, this->inp_->abs_broadening);
            spectrum.write_transition_dipole(this->out_dir +
                "trans_dipole_" + spin_types[is] + "_tda.dat");
            // =============================================== for test ====================================================
            // spectrum.optical_absorption_method2(freq, this->inp_->abs_broadening);
            // if (LR_Util::tolower(this->inp_->abs_gauge) == "velocity")
            // {   // TEST the formula v/omega rather than v/(e_a-e_i)
            //     spectrum.test_transition_dipoles_velocity_omega();
            //     spectrum.write_transition_dipole(this->out_dir +
            //         "trans_dipole_" + spin_types[is] + "_vomega_tda.dat");
            // }
            // =============================================== for test ====================================================
        }
        if (this->inp_->cal_force && !this->excited_relax_) { this->cal_force_and_grad_matrix_(is, this->ofs_running_); }
    }
}
template class ModuleESolver::ESolver_LR<double, double>;
template class ModuleESolver::ESolver_LR<std::complex<double>, double>;
