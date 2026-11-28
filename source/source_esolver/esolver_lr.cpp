#include "esolver_lr.h"
#include "source_hamilt/module_xc/xc_functional.h"
#include "source_lcao/module_lr/hamilt_ulr.h"
#include "source_lcao/module_lr/hsolver_lrtd.h"
#include "source_lcao/module_lr/lr_density.h"
#include "source_lcao/module_lr/lr_spectrum.h"
#include "source_lcao/module_lr/operator_casida/operator_lr_diag.h"
#include "source_lcao/module_lr/ri_benchmark/ri_benchmark.h"
#include "source_lcao/module_lr/utils/lr_io.h"
#include "source_lcao/module_lr/utils/lr_util.h"
#include "source_lcao/module_lr/utils/lr_util_print.h"
#include "source_lcao/module_lr/utils/spectrum_weights.h"
#include "source_lcao/module_lr/zeq_solver.h"
#include <algorithm>
#include <iomanip>
#include "source_lcao/module_lr/kernel_pots.h"
#include "source_io/module_parameter/parameter.h"
#include "source_lcao/module_lr/lr_dim.h"
#include <type_traits>
#include <vector>
#include "source_base/tool_quit.h"
#include "source_base/tool_title.h"
#include "source_lcao/module_lr/lr_gs.h"
#include <cassert>
#include <complex>
#include <memory>
#include "source_estate/module_dm/dm_from_psi.h"
#include "source_io/module_output/output_log.h"
#include "source_lcao/module_lr/grad_degen.h"
#include "source_lcao/module_lr/gradient_inputs.h"
#include "source_lcao/module_lr/gradient_output.h"
#include "source_lcao/module_lr/lr_amp.h"
#include "source_lcao/module_lr/lr_dm.h"
#include "source_lcao/module_lr/lr_force.h"
#include "source_base/timer.h"
#include "source_lcao/module_lr/init_helpers.h"
#include "source_base/module_external/scalapack_connector.h"
#include "source_basis/module_pw/pw_basis_big.h"
#include "source_cell/module_neighbor/sltk_atom_arrange.h"
#include "source_hamilt/module_gint/gint.h"
#include "source_io/module_output/print_info.h"
#ifdef __JSON
#include "source_io/module_json/output_info.h"
#endif
#ifdef __EXX
#include "source_hamilt/module_xc/exx_info.h"
#include "source_lcao/module_ri/exx_lri_interface.h"
#endif

using namespace LR;

// LR execution flow and final spectrum/density output.


template <typename T, typename TR>
ModuleESolver::ESolver_LR<T, TR>::ESolver_LR(
    const Input_para &inp,
    const std::string &in_dir,
    const std::string &out_dir,
    const LRRuntime &runtime)
    : in_dir(in_dir), out_dir(out_dir), ofs_running_(runtime.running_log), ofs_warning_(runtime.warning_log),
      my_rank_(runtime.rank)
{
#ifdef __EXX
    init_exx_info(this->exx_info, inp);
#endif
}

template <typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::runner(BaseCell &basecell, const int istep)
{
    basecell.require_kind(BaseCell::Kind::unitcell, __FUNCTION__);
    UnitCell &ucell = static_cast<UnitCell &>(basecell);

    ModuleBase::TITLE("ESolver_LR", "runner");
    ModuleBase::timer::start("ESolver_LR", "runner");

    if (this->ks_)
    {
        // `init_pot_groundstate` leaves the global XC type on the LR kernel, so put it back before
        // the ground-state SCF. `ESolver_KS::runner` is re-entrant: its `before_scf` rebuilds the
        // neighbour lists, the grid tables and the Hamiltonian for the geometry of this step.
        XC_Functional::set_xc_type(ucell.atoms[0].ncpp.xc_func);
        this->ks_->runner(ucell, istep);
        this->conv_esolver = this->ks_->conv_esolver;
        if (this->inp_->calculation == "relax" && !this->conv_esolver)
        {
            ModuleBase::WARNING_QUIT(
                "ESolver_LR", "Excited-state relaxation requires a converged ground-state SCF; ionic step rejected");
        }
        this->lr_.relax.etot_gs = this->ks_->cal_energy();
        if (!this->ks_initialized_)
        {
            this->initialize_from_ks_(ucell, *this->inp_);
            this->lr_.setup_relax(); // needs the final dimensions, i.e. `openshell`
        }
        else
        {
            this->refresh_from_ks_(ucell);
        }
    }

    // allocate 2-particle state and setup 2d division
    this->setup_eigenvectors_X();
    this->pelec->ekb.create(this->lr_.dim.nspin, this->lr_.dim.nstates);

    // set once here rather than inside the solver branches: the open-shell `spectrum` branch used
    // to leave it empty, and both `after_all_runners` and the gradient index it
    this->lr_.amp.spin_types =
        this->lr_.dim.openshell ? std::vector<std::string>({"updown"}) : std::vector<std::string>({"singlet", "triplet"});

    // a relaxation writes these once per ionic step; without the suffix every step would overwrite
    // the last, and the trajectory would be impossible to inspect afterwards
    const std::string step_suffix = this->lr_.relax.target.active() ? "_step" + std::to_string(istep) : "";
    auto efile_out = [this, &step_suffix](const std::string &label) -> std::string
    { return this->out_dir + "Excitation_Energy_" + label + step_suffix + ".dat"; };
    auto vfile_out = [this, &step_suffix](const std::string &label) -> std::string
    {
        return this->out_dir + "Excitation_Amplitude_" + label + step_suffix + "_" +
               std::to_string(this->my_rank_ + 1) + ".dat";
    };
    if (this->inp_->lr_solver == "elpa")
    {
        ModuleBase::WARNING_QUIT("ESolver_LR", "ESolver_LR doesn't support elpa now.");
    }
    else if (this->inp_->lr_solver != "spectrum")
    {
        auto write_states = [this, &efile_out, &vfile_out](
                                const std::string &label, const Real<T> *e, const T *v, const int &dim, const int &nst,
                                const int &prec = 8) -> void
        {
            if (this->my_rank_ == 0)
            {
                assert(nst == LR_Util::write_value(efile_out(label), prec, e, nst));
            }
            assert(nst * dim == LR_Util::write_value(vfile_out(label), prec, v, nst, dim));
        };
        std::vector<double> precondition(this->inp_->lr_solver == "lapack" ? 0 : this->lr_.amp.nloc, 1.0);
        // allocate and initialize A matrix and density matrix
        if (this->lr_.dim.openshell)
        {
            for (int is : {0, 1})
            {
                if (this->inp_->lr_solver != "lapack")
                {
                    const int offset_is = is * this->lr_.amp.para_x[0].size();
                    OperatorLRDiag<double> pre_op(
                        this->lr_.gs.eig.c + is * this->lr_.dim.nk * (this->lr_.dim.nocc[0] + this->lr_.dim.nvirt[0]), this->lr_.amp.para_x[is], this->lr_.dim.nk, this->lr_.dim.nocc[is],
                        this->lr_.dim.nvirt[is], false);
                    pre_op.act(1, offset_is, 1, precondition.data() + offset_is, precondition.data() + offset_is);
                }
            }
            std::cout << "Solving spin-conserving excitation for open-shell system." << std::endl;
            HamiltULR<T> hulr(
                xc_kernel, this->lr_.dim.nspin, this->lr_.dim.nbasis, this->lr_.dim.nocc, this->lr_.dim.nvirt, *this->ucell_, orb_cutoff_, this->gd(),
                *this->lr_.gs.psi, this->lr_.gs.eig,
#ifdef __EXX
                this->exx_lri, this->exx_info.info_global.hybrid_alpha,
#endif
                this->lr_.pot.pot, this->kv, this->lr_.amp.para_x, this->lr_.amp.para_c, this->lr_.gs.para_mat);
            LR::HSolver::solve(
                hulr, this->lr_.amp.x[0].template data<T>(), this->lr_.amp.nloc, this->lr_.dim.nstates, this->lr_.dim.nk, this->lr_.dim.nocc, this->lr_.dim.nvirt,
                this->lr_.amp.para_x, this->pelec->ekb.c, this->inp_->lr_solver, this->inp_->lr_thr, precondition, true);
            if (this->inp_->out_wfc_lr)
            {
                write_states("openshell", this->pelec->ekb.c, this->lr_.amp.x[0].template data<T>(), this->lr_.amp.nloc, this->lr_.dim.nstates);
            }
        }
        else
        {
            if (this->inp_->lr_solver != "lapack")
            {
                OperatorLRDiag<double> pre_op(
                    this->lr_.gs.eig.c, this->lr_.amp.para_x[0], this->lr_.dim.nk, this->lr_.dim.nocc[0], this->lr_.dim.nvirt[0], false);
                pre_op.act(1, this->lr_.amp.nloc, 1, precondition.data(), precondition.data());
            }
            const std::vector<std::string> &spin_types = this->lr_.amp.spin_types;
            for (int is = 0; is < this->lr_.dim.nspin; ++is)
            {
                std::cout << " Calculating " << spin_types[is] << " excitations" << std::endl;
                HamiltLR<T> hlr(
                    xc_kernel, this->lr_.dim.nspin, this->lr_.dim.nbasis, this->lr_.dim.nocc, this->lr_.dim.nvirt, *this->ucell_, orb_cutoff_, this->gd(),
                    *this->lr_.gs.psi, this->lr_.gs.eig,
#ifdef __EXX
                    this->exx_lri, this->exx_info.info_global.hybrid_alpha,
#endif
                    this->lr_.pot.pot[is], this->kv, this->lr_.amp.para_x, this->lr_.amp.para_c, this->lr_.gs.para_mat, spin_types[is], this->in_dir,
                    this->out_dir, this->inp_->ri_hartree_benchmark,
                    (this->inp_->ri_hartree_benchmark == "aims" ? this->inp_->aims_nbasis : std::vector<int>({})));
                LR::HSolver::solve(
                    hlr, this->lr_.amp.x[is].template data<T>(), this->lr_.amp.nloc, this->lr_.dim.nstates, this->lr_.dim.nk, this->lr_.dim.nocc, this->lr_.dim.nvirt,
                    this->lr_.amp.para_x, this->pelec->ekb.c + is * this->lr_.dim.nstates, this->inp_->lr_solver, this->inp_->lr_thr,
                    precondition, true);
                if (this->inp_->out_wfc_lr)
                {
                    write_states(
                        spin_types[is], this->pelec->ekb.c + is * this->lr_.dim.nstates, this->lr_.amp.x[is].template data<T>(),
                        this->lr_.amp.nloc, this->lr_.dim.nstates);
                }
            }
        }
    }
    else // lr_solver == "spectrum", read the eigenvalues
    {
        auto efile_in = [&](const std::string &label) -> std::string
        { return this->in_dir + "Excitation_Energy_" + label + ".dat"; };
        auto vfile_in = [this](const std::string &label) -> std::string
        { return this->in_dir + "Excitation_Amplitude_" + label + "_" + std::to_string(this->my_rank_ + 1) + ".dat"; };

        auto read_states = [&](const std::string &label, Real<T> *e, T *v, const int &dim, const int &nst) -> void
        {
            if (this->my_rank_ == 0)
            {
                assert(nst == LR_Util::read_value(efile_in(label), e, nst));
                std::cout << "Rank " << this->my_rank_ << ": finish reading " << efile_in(label) << std::endl;
            }
#ifdef __MPI
            // in velocity gauge, the eigenvalues may be used to calculate the transition dipole, so we'd better
            // broadcast them
            MPI_Bcast(e, nst, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
            assert(nst * dim == LR_Util::read_value(vfile_in(label), v, nst, dim));
            std::cout << "Rank " << this->my_rank_ << ": finish reading " << vfile_in(label) << std::endl;
        };
        std::cout << "reading the excitation states from file: \n";
        if (this->lr_.dim.openshell)
        {
            read_states("openshell", this->pelec->ekb.c, this->lr_.amp.x[0].template data<T>(), this->lr_.amp.nloc, this->lr_.dim.nstates);
        }
        else
        {
            const std::vector<std::string> &spin_types = this->lr_.amp.spin_types;
            for (int is = 0; is < this->lr_.dim.nspin; ++is)
            {
                read_states(
                    spin_types[is], this->pelec->ekb.c + is * this->lr_.dim.nstates, this->lr_.amp.x[is].template data<T>(), this->lr_.amp.nloc,
                    this->lr_.dim.nstates);
            }
        }
    }
    if (this->lr_.relax.target.active())
    {
        // Re-select which root to follow BEFORE taking the gradient, so the force belongs to
        // the same diabatic state as the previous step's.
        this->lr_.follow_target(this->pelec->ekb, this->ofs_running_);

        // The LR terms only carry the Omega part of the force; the ground-state part is separate
        // and comes straight from the KS solver.
        this->ks_->cal_force(ucell, this->lr_.relax.force_gs);
        this->lr_.relax.lr_force = this->lr_.cal_lr_force(this->pelec->ekb, this->ofs_running_);

        // One line per ionic step with the two halves of the energy and of the gradient. Without
        // it the relaxation only reports a force, and whether E_gs + Omega actually goes down --
        // the thing being minimised -- cannot be read off the log at all.
        const double omega = this->lr_.relax.target.omega(this->pelec->ekb.c, this->lr_.dim.nstates);
        auto max_abs = [](const ModuleBase::matrix &m) -> double
        {
            double v = 0.0;
            for (int i = 0; i < m.nr * m.nc; ++i)
            {
                v = std::max(v, std::abs(m.c[i]));
            }
            return v;
        };
        this->ofs_running_ << std::setprecision(8) << std::fixed << " EXCITED-STATE RELAX step " << istep
                           << ": E_gs = " << this->lr_.relax.etot_gs * ModuleBase::Ry_to_eV
                           << " eV, Omega = " << omega * ModuleBase::Ry_to_eV
                           << " eV, E_exc = " << (this->lr_.relax.etot_gs + omega) * ModuleBase::Ry_to_eV << " eV"
                           << std::setprecision(6) << " | |F_gs|max = "
                           << max_abs(this->lr_.relax.force_gs) * ModuleBase::Ry_to_eV / ModuleBase::BOHR_TO_A
                           << ", |F_Omega|max = "
                           << max_abs(this->lr_.relax.lr_force) * ModuleBase::Ry_to_eV / ModuleBase::BOHR_TO_A << " eV/Angstrom"
                           << std::defaultfloat << std::endl;
    }

    ModuleBase::timer::end("ESolver_LR", "runner");
    return;
}

template <typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::after_all_runners(BaseCell &basecell)
{
    basecell.require_kind(BaseCell::Kind::unitcell, __FUNCTION__);
    UnitCell &ucell = static_cast<UnitCell &>(basecell);

    ModuleBase::TITLE("ESolver_LR", "after_all_runners");
    if (this->inp_->ri_hartree_benchmark != "none")
    {
        return;
    } // no need to calculate the spectrum in the benchmark routine

    // cal electron-hole density
    if (this->inp_->out_chg[0])
    {
        LR_Density<T> lr_density(
            *this->ucell_, kv, gd(), *this->lr_.gs.psi, orb_cutoff_, pgrid(), this->lr_.dim.nspin, this->lr_.dim.nocc, this->lr_.dim.nvirt, this->lr_.dim.nbasis, this->lr_.amp.para_x, this->lr_.amp.para_c,
            paraMat_, this->lr_.dim.openshell, this->inp_->out_chg[1], this->out_dir);

        if (this->lr_.dim.openshell)
        {
            for (int is = 0; is < this->lr_.dim.nspin; ++is)
            {
                const T *const amplitudes = this->lr_.amp.x[0].template data<T>();
                lr_density.output_eh_density_all_states(amplitudes, is, this->lr_.dim.nstates, this->ofs_running_);
            }
        }
        else
        {
            for (int is = 0; is < this->lr_.amp.x.size(); ++is)
            {
                const T *const amplitudes = this->lr_.amp.x[is].template data<T>();
                lr_density.output_eh_density_all_states(amplitudes, is, this->lr_.dim.nstates, this->ofs_running_);
            }
        }
    }

    const std::string gauge = LR_Util::tolower(this->inp_->abs_gauge);
    const bool calculate_spectrum = gauge != "none";
    const std::vector<double> spatial_weights = LR_Util::spatial_k_weights(this->kv.wk, this->lr_.dim.nspin);
    // cal spectrum
    if (LR_Util::tolower(this->inp_->abs_gauge) == "velocity")
    {
        const int nspin_tmp = this->inp_->nspin == 2 ? 2 : 1;
        this->velocity_mo = LR_Util::cal_velocity_mo(
            *this->ucell_, this->gd(), this->tcb(), this->lr_.gs.para_mat, this->lr_.amp.para_c, this->kv, *this->lr_.gs.psi, this->lr_.dim.nk,
            nspin_tmp, this->lr_.dim.nbasis, this->lr_.dim.nocc, this->lr_.dim.nvirt);
    }
    std::vector<double> freq;
    if (calculate_spectrum)
    {
        freq.resize(100);
        std::vector<double> abs_wavelen_range({20, 200}); // default range
        if (this->inp_->abs_wavelen_range.size() >= 2 &&
            std::abs(this->inp_->abs_wavelen_range[1] - this->inp_->abs_wavelen_range[0]) > 0.02)
        {
            abs_wavelen_range = this->inp_->abs_wavelen_range;
        }
        double lambda_diff = std::abs(abs_wavelen_range[1] - abs_wavelen_range[0]);
        double lambda_min = std::min(abs_wavelen_range[1], abs_wavelen_range[0]);
        for (int i = 0; i < freq.size(); ++i)
        {
            freq[i] = 91.126664 / (lambda_min + 0.01 * static_cast<double>(i + 1) * lambda_diff);
        }
    }
    // auto spin_types = (nspin == 2 && !openshell) ? std::vector<std::string>({ "singlet", "triplet" }) :
    // std::vector<std::string>({ "updown" }); for (int is = 0;is < this->lr_.amp.x.size() - 1;++is)
    for (int is = 0; is < this->lr_.amp.x.size(); ++is)
    {
        if (!calculate_spectrum)
        {
            const double *const energies = &this->pelec->ekb.c[is * this->lr_.dim.nstates];
            LR_Util::print_excitation_energies(this->ofs_running_, energies, this->lr_.dim.nstates, this->lr_.amp.spin_types[is]);
        }
        else
        {
            LR_Spectrum<T> spectrum(
                this->lr_.dim.nspin, this->lr_.dim.nbasis, this->lr_.dim.nocc, this->lr_.dim.nvirt, *this->pw_rho, *this->lr_.gs.psi, *this->ucell_, this->kv,
                spatial_weights, this->gd(), this->orb_cutoff_, this->tcb(), this->lr_.amp.para_x, this->lr_.amp.para_c, this->lr_.gs.para_mat,
                &this->pelec->ekb.c[is * this->lr_.dim.nstates], this->lr_.gs.eig.c, this->lr_.amp.x[is].template data<T>(), this->lr_.dim.nstates, this->lr_.dim.openshell,
                LR_Util::tolower(this->inp_->abs_gauge), this->my_rank_, this->out_dir);
            if (LR_Util::tolower(this->inp_->abs_gauge) == "velocity")
            {
                spectrum.set_vmo(this->velocity_mo.data());
            }
            spectrum.cal_spectrum();
            spectrum.transition_analysis(this->lr_.amp.spin_types[is] + "_tda");
            if (this->lr_.amp.spin_types[is] != "triplet") // triplets has no transition dipole and no contribution to the spectrum
            {
                spectrum.optical_absorption_method1(freq, this->inp_->abs_broadening);
                spectrum.write_transition_dipole(this->out_dir + "trans_dipole_" + this->lr_.amp.spin_types[is] + "_tda.dat");
                // =============================================== for test
                // ==================================================== spectrum.optical_absorption_method2(freq,
                // this->inp_->abs_broadening); if (LR_Util::tolower(this->inp_->abs_gauge) == "velocity") {   // TEST
                // the formula v/omega rather than v/(e_a-e_i)
                //     spectrum.test_transition_dipoles_velocity_omega();
                //     spectrum.write_transition_dipole(this->out_dir +
                //         "trans_dipole_" + spin_types[is] + "_vomega_tda.dat");
                // }
                // =============================================== for test
                // ====================================================
            }
        }
        if (this->inp_->cal_force && !this->lr_.relax.target.active())
        {
            this->lr_.cal_force_grad(is, this->pelec->ekb, this->ofs_running_);
        }
    }
}

// Resolve LR dimensions and initialize kernel potentials.

template <typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::set_dim()
{
    LR::DimSrc src;
    src.inp = this->inp_;
    src.nlocal = PARAM.globalv.nlocal;
    src.nks = this->kv.get_nks();
    src.nkstot = this->kv.get_nkstot();
    src.rank = this->my_rank_;
    src.gamma_only = std::is_same<T, double>::value;
    src.wg_ks = this->ks_ ? &this->ks_->pelec->wg : nullptr;
    src.kpar = this->kv.para_k.kpar;
    src.rank_in_pool = this->kv.para_k.rank_in_pool;
    src.in_dir = this->in_dir;
    src.ucell = const_cast<UnitCell *>(this->ucell_);
    src.log = &this->ofs_running_;

    this->lr_.dim = LR::resolve_dim(src);
}

template <typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::reset_dim_spin()
{
    LR::apply_spin_pol(this->lr_.dim, *this->inp_);
}

template <typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::init_pot(const Charge &chg_gs)
{
    this->lr_.pot.pot.assign(this->lr_.dim.nspin, nullptr);
    if (this->inp_->ri_hartree_benchmark != "none")
    {
        return;
    }
    const auto &grid = pgrid();
    auto prepared =
        LR::prepare_lr_potentials(*this->inp_, xc_kernel, *this->pw_rho, *this->ucell_, chg_gs, grid, this->lr_.dim.nspin, this->lr_.dim.openshell);
    this->lr_.pot.pot = std::move(prepared.pot);
    if (this->inp_->cal_force)
    {
        this->init_pot_groundstate(chg_gs);
        this->lr_.pot.pot_hxc_gs = LR::prepare_gs_potential(
            *this->inp_, xc_kernel, *this->pw_rho, *this->ucell_, chg_gs, grid, this->lr_.dim.nspin, this->lr_.dim.openshell, prepared);
    }
}

// LR orbital distributions, excitation amplitudes, and ground-state file input.


template <typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::setup_eigenvectors_X()
{
    ModuleBase::TITLE("ESolver_LR", "setup_eigenvectors_X");
    LR::XAmp<T> amp = LR::setup_x<T>(this->lr_.amp.para_c, this->lr_.dim.nocc, this->lr_.dim.nvirt, this->lr_.dim.npairs, this->lr_.gs.eig, this->lr_.dim.nstates,
                                     this->lr_.dim.nk, this->lr_.dim.nspin, this->lr_.dim.openshell, this->inp_->lr_solver == "spectrum",
                                     this->ofs_running_);
    this->lr_.amp.x = std::move(amp.x);
    this->lr_.amp.para_x = std::move(amp.para_x);
    this->lr_.amp.nloc = amp.nloc;
}

template <typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::read_ks_wfc()
{
    this->lr_.read_wfc(this->in_dir, *this->inp_, this->kv, this->ofs_running_);
}

template <typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::fill_z_window_(const int *desc_src)
{
    this->lr_.fill_z(desc_src, this->kv, this->inp_->cal_force, !this->ks_initialized_, this->ofs_running_);
}

template <typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::read_ks_chg(Charge &chg_gs)
{
    LR::read_gs_chg(this->in_dir, *this->inp_, this->lr_.dim.nspin, *this->ucell_, pgrid(), this->pw_rho, chg_gs,
                    this->my_rank_, this->ofs_running_);
}


///========================= excited-state geometry relaxation =========================

template <typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::init_pot_groundstate(const Charge &chg_gs)
{
    ModuleBase::TITLE("ESolver_LR", "init_pot_gs");
    this->lr_.pot.gs_pot.refresh(*this->inp_, this->xc_kernel, *this->ucell_, this->pw_rhod, this->pw_rho, pgrid(), chg_gs,
                          this->sfac(), this->vloc(), this->solvent, !this->ks_, PARAM.globalv.has_float_data);
}


template <typename T, typename TR>
double ModuleESolver::ESolver_LR<T, TR>::cal_energy()
{
    // Outside a relaxation nothing consumes this, and returning a non-zero value would change
    // what the existing single-point outputs report.
    if (!this->lr_.relax.target.active())
    {
        return 0.0;
    }
    // `omega` is the multiplet average under `lr_degen_mode = average` and the single state
    // otherwise, matching whatever `cal_lr_force_relax_` produced the gradient of. The
    // energy-based optimisers line-search on this, so the two must not describe different
    // surfaces.
    return this->lr_.relax.etot_gs + this->lr_.relax.target.omega(this->pelec->ekb.c, this->lr_.dim.nstates);
}

template <typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::cal_force(BaseCell &basecell, ModuleBase::matrix &force)
{
    basecell.require_kind(BaseCell::Kind::unitcell, __FUNCTION__);
    const UnitCell &ucell = static_cast<const UnitCell &>(basecell);
    if (!this->lr_.relax.target.active())
    { // single-point runs print the gradients of every state from `after_all_runners` instead
        return;
    }
    if (this->lr_.relax.lr_force.nr != ucell.nat)
    {
        ModuleBase::WARNING_QUIT(
            "ESolver_LR::cal_force", "the excited-state gradient has not been computed for this geometry.");
    }
    // Both halves already follow the ABACUS force convention F = -dE/dR (Ry/Bohr), so they add.
    // The "Gradients of each excited state" heading that `cal_force(int)` prints under is a
    // misnomer: the finite-difference reference it was validated against (`abacus-fd lr-custom`)
    // computes (E(-h) - E(+h))/h, which is -d(Omega)/dR, and the two agree in sign.
    force.create(ucell.nat, 3);
    force = this->lr_.relax.force_gs + this->lr_.relax.lr_force;
    ModuleIO::print_force(this->ofs_running_, ucell, "EXCITED-STATE TOTAL-FORCE (eV/Angstrom)", force, false);
}

template <typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::cal_stress(BaseCell &basecell, ModuleBase::matrix &stress)
{
    static_cast<void>(stress);
    basecell.require_kind(BaseCell::Kind::unitcell, __FUNCTION__);
    ModuleBase::WARNING_QUIT(
        "ESolver_LR::cal_stress", "the excited-state stress is not implemented (every LR gradient stress term is a "
                                  "dummy passed with isstress=false).");
}

// Initialize LR from an embedded KS solver or standalone ground-state files.

#ifdef __EXX
namespace
{

/// One `Exx_LRI` carries ONE Coulomb operator, and it may be needed for two different
/// reasons: the LR kernel (when `xc_kernel` is a hybrid) and the ground-state force (when
/// `dft_functional` is a hybrid). That operator is NOT chosen here -- `Exx_LRI` reads
/// `info_ri.coulomb_param`, which `input_conv` builds from `dft_functional` alone. So when
/// the two disagree, the kernel silently gets the ground state's screening, not its own.
///
/// This used to be decided by a local `exx_ccp_type()` returning Erfc for "hse" and bare
/// Hf otherwise, written onto `info_global.ccp_type`. With general range-separated hybrids
/// that two-way split is not even expressible ($\alpha/r+\beta\,\mathrm{erfc}(\mu r)/r$
/// is both at once), and the write was in any case dead for the RI path: only `Exx_LRI`
/// runs here, and it never looks at `ccp_type`.
void warn_if_kernel_differs_from_gs(
    const std::string &xc_kernel, const std::string &dft_functional, std::ofstream &ofs_running)
{
    const bool k = LR::exx_kernel_list().count(xc_kernel) > 0;
    const bool g = LR::exx_kernel_list().count(dft_functional) > 0;
    if (k && xc_kernel != dft_functional)
    {
        ofs_running << " WARNING: xc_kernel (" << xc_kernel << ") and dft_functional (" << dft_functional
                    << ") are not the same functional. The Coulomb operator of"
                       " Exx_LRI follows dft_functional"
                    << (g ? ""
                          : ", which is not a hybrid at all"
                            " (no coulomb_param, so the LR exchange kernel vanishes)")
                    << "; set both to the"
                       " same hybrid to get the kernel you asked for."
                    << std::endl;
    }
}
} // namespace
#endif

template <typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::before_all_runners(BaseCell &basecell, const Input_para &inp)
{
    basecell.require_kind(BaseCell::Kind::unitcell, __FUNCTION__);
    UnitCell &ucell = static_cast<UnitCell &>(basecell);
    this->ucell_ = &ucell;
    this->inp_ = &inp;
    if (inp.esolver_type == "ks-lr")
    {
#ifdef __JSON
        // The embedded KS run happens before Relax_Driver starts its first step.
        Json::init_output_array_obj();
#endif
        // the ground-state solver is a member, not a temporary: `runner` re-runs its SCF on every
        // ionic step, and the objects aliased from it have to stay alive as long as it does.
        // Its SCF is deliberately NOT run here -- `before_all_runners` is called once, outside the
        // relaxation loop, so the ground state has to be recomputed from `runner(istep)` instead.
        this->ks_ = LR_Util::make_unique<ModuleESolver::ESolver_KS_LCAO<T, TR>>();
        this->ks_->before_all_runners(basecell, inp);
        LR_Util::check_force_pp(ucell, inp.cal_force, inp.calculation);
    }
    else
    {
        this->initialize_from_unitcell_(ucell, inp);
    }
}

template <typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::bind_ground_state_aliases_()
{
    if (this->ks_)
    {
        // `ESolver_FP::before_all_runners` was never called on this object, so its own `pw_rhod`,
        // `Pgrid`, `sf` and `locpp` are empty -- the excited-state force needs all four. Point at
        // the ground-state solver's, which `before_scf` refreshes for the current geometry every
        // ionic step. `pw_rho_flag` stays false so the destructor does not free what it borrowed.
        this->pw_rho = this->ks_->pw_rho;
        this->pw_rhod = this->ks_->pw_rhod;
        this->pw_big = this->ks_->pw_big;
        this->pw_rho_flag = false;
        this->pgrid_ptr_ = &this->ks_->Pgrid;
        this->sf_ptr_ = &this->ks_->sf;
        this->locpp_ptr_ = &this->ks_->locpp;
        this->gd_ptr_ = &this->ks_->gd;
        // the two-center tables depend on the orbitals only, so the ground-state solver's single
        // build serves every geometry. This used to be moved over only for `abs_gauge velocity`,
        // leaving `LR_Force` and `cal_hs_grad` with an empty bundle in the length gauge.
        this->tcb_ptr_ = &this->ks_->two_center_bundle_;
    }
    else
    {
        this->pgrid_ptr_ = &this->Pgrid;
        this->sf_ptr_ = &this->sf;
        this->locpp_ptr_ = &this->locpp;
        this->gd_ptr_ = &this->gd_own_;
        this->tcb_ptr_ = &this->two_center_bundle_own_;
    }
}

template <typename T, typename TR>
LR::Refs<T> ModuleESolver::ESolver_LR<T, TR>::make_refs() const
{
    LR::Refs<T> refs;
    refs.ucell = this->ucell_;
    refs.kv = &this->kv;
    refs.gd = &this->gd();
    refs.orb_cutoff = &this->orb_cutoff_;
    refs.rho = this->pw_rho;
    refs.rhod = this->pw_rhod;
    refs.vloc = &this->vloc();
    refs.sfac = &this->sfac();
    refs.tcb = &this->tcb();
    refs.inp = this->inp_;
    refs.log = &this->ofs_running_;
    refs.pgrid = &this->pgrid();
    refs.in_dir = &this->in_dir;
    refs.out_dir = &this->out_dir;
    refs.xc_kernel = &this->xc_kernel;
    refs.my_rank = this->my_rank_;
    refs.pv_all = this->ks_ ? &this->ks_->pv : &this->lr_.gs.para_mat_all;
#ifdef __EXX
    refs.exx_lri = &this->exx_lri;
    refs.hybrid_alpha = this->exx_info.info_global.hybrid_alpha;
#endif
    return refs;
}

template <typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::initialize_from_ks_(UnitCell &ucell, const Input_para &inp)
{
    ModuleBase::TITLE("ESolver_LR", "ESolver_LR(KS)");
    ModuleESolver::ESolver_KS_LCAO<T, TR> &ks_sol = *this->ks_;
    this->bind_ground_state_aliases_();

    if (this->inp_->lr_solver == "spectrum")
    {
        throw std::invalid_argument("when lr_solver==spectrum, esolver_type must be `lr` to skip KS calculation.");
    }

    // xc kernel
    this->xc_kernel = LR_Util::tolower(inp.xc_kernel);
    // kv
    this->kv = ks_sol.kv; // copy: cheap, and the KS solver keeps using its own

    LR::check_lr_parameters(*this->inp_, this->xc_kernel, this->lr_.dim.nspin);

    this->set_dim();

    // KS and LR must use the same AO distribution, including its BLACS grid.
    const int ks_block_size = ks_sol.pv.para_2d().nb();
    LR_Util::setup_2d_division(
        this->lr_.gs.para_mat.para_2d(), ks_block_size, this->lr_.dim.nbasis, this->lr_.dim.nbasis
#ifdef __MPI
        ,
        ks_sol.pv.para_2d().ictxt()
#endif
    );
    LR::set_band_layout(this->lr_.gs.para_mat, this->lr_.dim.nbasis, this->lr_.dim.nbands, *this->ucell_,
                        this->inp_->ri_hartree_benchmark == "aims");
    if (this->inp_->cal_force)
    {
        LR_Util::setup_2d_division(
            this->lr_.gs.para_mat_all.para_2d(), ks_block_size, this->lr_.dim.nbasis, this->lr_.dim.nbasis
#ifdef __MPI
            ,
            ks_sol.pv.para_2d().ictxt()
#endif
        );
        LR::set_band_layout(this->lr_.gs.para_mat_all, this->lr_.dim.nbasis, this->inp_->nbands, *this->ucell_,
                            this->inp_->ri_hartree_benchmark == "aims");
    }

    this->lr_.gs.para_mat.atom_begin_row = ks_sol.pv.atom_begin_row;
    this->lr_.gs.para_mat.atom_begin_col = ks_sol.pv.atom_begin_col;
    this->lr_.gs.para_mat.iat2iwt_ = ucell.get_iat2iwt();

    LR_Util::setup_2d_division(
        this->lr_.amp.para_c, ks_block_size, this->lr_.dim.nbasis, this->lr_.dim.nbands
#ifdef __MPI
        ,
        this->lr_.gs.para_mat.para_2d().ictxt()
#endif
    );

    // allocate psi_ks and eig_ks in the [nocc, nvirt] window
#ifdef __MPI
    this->lr_.gs.psi.reset(
        new psi::Psi<T>(this->kv.get_nks(), this->lr_.amp.para_c.ncloc(), this->lr_.amp.para_c.nrloc(), this->kv.ngk, true));
#else
    this->lr_.gs.psi.reset(new psi::Psi<T>(this->kv.get_nks(), this->lr_.dim.nbands, this->lr_.dim.nbasis, this->kv.ngk, true));
#endif
    this->lr_.gs.eig.create(this->kv.get_nks(), this->lr_.dim.nbands);
    this->pelec = new elecstate::ElecStateLCAO<T>();
    orb_cutoff_ = ks_sol.orb_.cutoffs();

#ifdef __EXX
    // Two independent reasons to need an Exx_LRI: the LR exchange kernel, and the ground-state
    // EXX terms of the gradient (the H_gs[T+Z] and W multipliers, gated on gs_is_hybrid).
    // `initialize_from_unitcell_` has always covered both; this path used to test only the first,
    // so a local kernel on top of a hybrid ground state had no exx_lri at all when asked for forces.
    if (exx_kernel_list().count(xc_kernel) || (this->inp_->cal_force && gs_is_hybrid(this->inp_->dft_functional)))
    {
        std::string dft_functional = LR_Util::tolower(this->inp_->dft_functional);
        // Either object would be built from the same `info_ri.coulomb_param`, which `input_conv`
        // derives from dft_functional alone -- so whenever the ground-state solver has one of the
        // right type its geometry tensors are already up to date. Share those tensors in a
        // separate LR electronic workspace, skipping a `cal_exx_ions` per ionic step.
        const bool share = (ks_sol.exx_nao.exd && std::is_same<T, double>::value) ||
                           (ks_sol.exx_nao.exc && std::is_same<T, std::complex<double>>::value);
        warn_if_kernel_differs_from_gs(xc_kernel, dft_functional, this->ofs_running_);
        if (share)
        {
            this->exx_owned_ = false;
        } // `refresh_from_ks_` refreshes the geometry snapshot
        else // construct C, V from scratch
        {
            // `input_conv` already filled `info_ri.coulomb_param` from INPUT.
            exx_info.sync_from_global();
            // populate ABFs/JLE file lists from UnitCell; keep in sync with Exx_NAO::init
            exx_info.info_ri.files_abfs = ucell.abfs_orbital_files;
            exx_info.info_opt_abfs.files_abfs = ucell.abfs_orbital_files;
            exx_info.info_opt_abfs.files_jles = ucell.jle_orbital_files;
            this->exx_lri = std::make_shared<Exx_LRI<T>>(exx_info.info_ri);
            this->exx_lri->init(MPI_COMM_WORLD, ucell, this->kv, ks_sol.orb_);
            this->exx_owned_ = true; // the position-dependent `cal_exx_ions` is left to `refresh_from_ks_`
        }
    }
#endif

    refresh_from_ks_(ucell);
    this->lr_.bind(this->make_refs());
    this->ks_initialized_ = true;
}

template <typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::refresh_from_ks_(UnitCell &ucell)
{
    ModuleBase::TITLE("ESolver_LR", "refresh_from_ks_");
    ModuleESolver::ESolver_KS_LCAO<T, TR> &ks_sol = *this->ks_;
    // The ground-state solver owns `psi` and reuses it across ionic steps, so it cannot be stolen.
    // `eig_ks_all` / `wg_ks_all` are nspin x nbands matrices -- a few kB, copied rather than aliased
    // so that they survive the KS solver overwriting `pelec` on the next step.
    this->lr_.gs.psi_all = ks_sol.psi;
    this->lr_.gs.eig_all = ks_sol.pelec->ekb;
    this->lr_.gs.wg_all = ks_sol.pelec->wg;
    const int max_occupied_bands = *std::max_element(this->lr_.dim.nocc.begin(), this->lr_.dim.nocc.end());
    const int start_band = this->lr_.dim.nocc_max - max_occupied_bands;

    for (int ik = 0; ik < this->kv.get_nks(); ++ik)
    {
        // copy the KS orbitals in the [nocc, nvirt] window
#ifdef __MPI
        this->lr_.gs.psi_all->fix_k(ik);
        this->lr_.gs.psi->fix_k(ik);
        const int first_band = start_band + 1;
        Cpxgemr2d(
            this->lr_.dim.nbasis, this->lr_.dim.nbands, this->lr_.gs.psi_all->get_pointer(), 1, first_band, ks_sol.pv.desc_wfc,
            this->lr_.gs.psi->get_pointer(), 1, 1, const_cast<int *>(this->lr_.amp.para_c.desc()), this->lr_.amp.para_c.ictxt());
#else
        // serial: each band is `nbasis` contiguous coefficients, so the window is a plain
        // band-by-band copy (this loop used to compute the two pointers and copy nothing,
        // leaving `psi_ks` uninitialized in every non-MPI build)
        for (int ib = 0; ib < this->lr_.dim.nbands; ++ib)
        {
            const auto *start = &(*this->lr_.gs.psi_all)(ik, start_band + ib, 0);
            auto *to = &(*this->lr_.gs.psi)(ik, ib, 0);
            std::copy(start, start + this->lr_.dim.nbasis, to);
        }
#endif
        // copy the KS bands in the [nocc, nvirt] window
        for (int ib = 0; ib < this->lr_.dim.nbands; ++ib)
        {
            this->lr_.gs.eig(ik, ib) = this->lr_.gs.eig_all(ik, start_band + ib);
        }
    }

    if (this->lr_.dim.nspin == 2)
    {
        const int nupdown_now =
            LR_Util::spin_population_difference(ks_sol.pelec->wg, this->kv.para_k.kpar, this->kv.para_k.rank_in_pool);
        if (!this->ks_initialized_)
        {
            this->lr_.dim.nupdown = nupdown_now;
            reset_dim_spin(); // shifts nocc/nvirt between the spin channels: must run exactly once
        }
        else if (nupdown_now != this->lr_.dim.nupdown)
        {
            ModuleBase::WARNING_QUIT(
                "ESolver_LR::refresh_from_ks_",
                "the ground-state spin population changed between ionic steps, but nocc/nvirt and"
                " the distributions built from them were fixed at the first step.");
        }
    }
    // `gint_info_` stays owned by the KS solver: its `before_scf` rebuilds and re-publishes it
    // init potential and calculate kernels using ground state charge
    init_pot(*ks_sol.pelec->charge);

#ifdef __EXX
    if (exx_kernel_list().count(xc_kernel) || (this->inp_->cal_force && gs_is_hybrid(this->inp_->dft_functional)))
    {
        if (this->exx_owned_)
        { // Cs/Vs follow the atoms, so they are rebuilt for every geometry
            this->exx_lri->cal_exx_ions(ucell, this->inp_->out_ri_cv);
        }
        else if (ks_sol.exx_nao.exd)
        {
            this->exx_lri = LR::make_lr_exx<T>(ks_sol.exx_nao.exd->exx_ptr, ucell, this->kv);
        }
        else if (ks_sol.exx_nao.exc)
        {
            this->exx_lri = LR::make_lr_exx<T>(ks_sol.exx_nao.exc->exx_ptr, ucell, this->kv);
        }
    }
#endif
    // the grid-integration tables hang off a static pointer that the ground-state solver
    // re-publishes in its `before_scf`; make sure it names the object we integrate on
    ModuleGint::Gint::set_gint_info(this->ks_->gint_info_.get());

    // the Z-vector window: after `reset_dim_spin`, so nocc/nvirt/openshell are final
#ifdef __MPI
    this->fill_z_window_(ks_sol.pv.desc_wfc);
#else
    this->fill_z_window_(nullptr);
#endif
}

template <typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::initialize_from_unitcell_(UnitCell &ucell, const Input_para &inp)
{
    ModuleBase::TITLE("ESolver_LR", "ESolver_LR(from scratch)");
    this->bind_ground_state_aliases_();
    // xc kernel
    this->xc_kernel = LR_Util::tolower(inp.xc_kernel);

    // necessary steps in ESolver_FP
    ESolver_FP::before_all_runners(ucell, inp);
    LR_Util::check_force_pp(ucell, inp.cal_force, inp.calculation);
    this->pelec = new elecstate::ElecStateLCAO<T>();

    // necessary steps in ESolver_KS::before_all_runners : symmetry and k-points
    if (ModuleSymmetry::Symmetry::symm_flag == 1)
    {
        const int cal_symm_repr[2] = {this->inp_->cal_symm_repr[0], this->inp_->cal_symm_repr[1]};
        ucell.symm.analy_sys(
            ucell.lat, ucell.atoms, ucell.nat, ucell.ntype, ucell.iat2it, ucell.iat2ia, ucell.itia2iat,
            this->ofs_running_, this->inp_->symmetry_prec, this->inp_->nspin, this->inp_->calculation, cal_symm_repr);
        ModuleBase::GlobalFunc::DONE(this->ofs_running_, "SYMMETRY");
    }
    const bool use_ibz = false;
    const bool gamma_only_local = PARAM.globalv.gamma_only_local;
    const double kspacing[3] = {this->inp_->kspacing[0], this->inp_->kspacing[1], this->inp_->kspacing[2]};
    const double koffset[3] = {this->inp_->koffset[0], this->inp_->koffset[1], this->inp_->koffset[2]};
    this->kv.set(
        ucell, ucell.symm, this->inp_->kpoint_file, this->inp_->nspin, ucell.G, ucell.latvec, this->ofs_running_,
        this->ofs_warning_, use_ibz, this->out_dir, gamma_only_local, kspacing, this->inp_->kmesh_type, koffset);
    ModuleBase::GlobalFunc::DONE(this->ofs_running_, "INIT K-POINTS");
    ModuleIO::print_parameters(ucell, this->kv, inp);

    LR::check_lr_parameters(*this->inp_, this->xc_kernel, this->lr_.dim.nspin);

    /// read orbitals and build the interpolation table
    two_center_bundle_own_.build_orb(ucell.ntype, ucell.orbital_fn.data(), inp.orbital_dir);

    LCAO_Orbitals orb;
    two_center_bundle_own_.to_LCAO_Orbitals(
        orb, inp.lcao_ecut, inp.lcao_dk, inp.lcao_dr, inp.lcao_rmax, inp.out_element_info, inp.cal_force);
    orb_cutoff_ = orb.cutoffs();
    // Standalone LR neighbor search needs projector metadata for every optical gauge.
    LR::setup_nonlocal_info(orb, ucell, *this->inp_, this->my_rank_, this->ofs_running_);
    const bool need_force_tables = this->inp_->cal_force || this->inp_->calculation == "relax";
    const bool need_velocity_tables = LR_Util::tolower(this->inp_->abs_gauge) == "velocity";
    if (need_force_tables || need_velocity_tables)
    {
        LR::setup_2center_table(this->two_center_bundle_own_, ucell, *this->inp_);
    }

    this->set_dim();
    //  setup 2d-block distribution for AO-matrix and KS wfc
    LR_Util::setup_2d_division(this->lr_.gs.para_mat.para_2d(), 1, this->lr_.dim.nbasis, this->lr_.dim.nbasis);
    LR::set_band_layout(this->lr_.gs.para_mat, this->lr_.dim.nbasis, this->lr_.dim.nbands, *this->ucell_,
                        this->inp_->ri_hartree_benchmark == "aims");
    if (this->inp_->cal_force)
    {
        LR_Util::setup_2d_division(this->lr_.gs.para_mat_all.para_2d(), 1, this->lr_.dim.nbasis, this->lr_.dim.nbasis);
        LR::set_band_layout(this->lr_.gs.para_mat_all, this->lr_.dim.nbasis, this->inp_->nbands, *this->ucell_,
                            this->inp_->ri_hartree_benchmark == "aims");
    }

    // read the ground state info
    // now ModuleIO::read_wfc_nao needs `Parallel_Orbitals` and can only read all the bands
    // it need improvement to read only the bands needed
    this->lr_.gs.psi.reset(new psi::Psi<T>(
        this->kv.get_nks(), this->lr_.gs.para_mat.ncol_bands, this->lr_.gs.para_mat.para_2d().nrloc(), this->kv.ngk, true));
    this->read_ks_wfc();

    if (this->lr_.dim.nspin == 2)
    {
        // Complete file populations were read before selecting the occupied window.
        reset_dim_spin();
    }
    // the Z-vector window: after `reset_dim_spin`, so nocc/nvirt/openshell are final
#ifdef __MPI
    this->fill_z_window_(paraMat_all_.desc_wfc);
#else
    this->fill_z_window_(nullptr);
#endif

    LR_Util::setup_2d_division(
        this->lr_.amp.para_c, 1, this->lr_.dim.nbasis, this->lr_.dim.nbands
#ifdef __MPI
        ,
        paraMat_.para_2d().ictxt()
#endif
    );

    // clear ks info, new elecstate for excition
    delete this->pelec; // the ElecStateLCAO allocated above, only needed while reading the KS data
    this->pelec = new elecstate::ElecState();

    // read the ground state charge density and calculate xc kernel
    pgrid().init(
        this->pw_rho->nx, this->pw_rho->ny, this->pw_rho->nz, this->pw_rho->nplane, this->pw_rho->nrxx, pw_big->nbz,
        pw_big->bz, GlobalV::NPROC);
    Charge chg_gs;
    if (this->inp_->ri_hartree_benchmark == "none")
    {
        this->read_ks_chg(chg_gs);
    }
    this->init_pot(chg_gs);

    // search adjacent atoms and init Gint
    double search_radius = -1.0;
    search_radius = atom_arrange::set_sr_NL(
        this->ofs_running_, this->inp_->out_level, orb.get_rcutmax_Phi(), ucell.infoNL->get_rcutmax_Beta(),
        PARAM.globalv.gamma_only_local);
    atom_arrange::search(
        PARAM.globalv.search_pbc, this->ofs_running_, this->gd(), *this->ucell_, search_radius,
        this->inp_->test_atom_input);
    gint_info_.reset(new ModuleGint::GintInfo(
        this->pw_big->nbx, this->pw_big->nby, this->pw_big->nbz, this->pw_rho->nx, this->pw_rho->ny, this->pw_rho->nz,
        0, 0, this->pw_big->nbzp_start, this->pw_big->nbx, this->pw_big->nby, this->pw_big->nbzp, orb.Phi, ucell,
        this->gd(), this->inp_->nspin, PARAM.globalv.gamma_only_local, PARAM.globalv.domag, this->inp_->device == "gpu",
        this->inp_->nstream));
    ModuleGint::Gint::set_gint_info(gint_info_.get());
    // if EXX from scratch, init 2-center integral and calculate Cs, Vs
    // when:
    // 1. EXX xc_kernel
    // 2. cal_force with ground state with EXX functional
#ifdef __EXX
    if (((exx_kernel_list().count(xc_kernel)) && this->inp_->lr_solver != "spectrum") ||
        (this->inp_->cal_force && gs_is_hybrid(this->inp_->dft_functional)))
    {
        warn_if_kernel_differs_from_gs(xc_kernel, LR_Util::tolower(this->inp_->dft_functional), this->ofs_running_);
        // `input_conv` already filled `info_ri.coulomb_param` from INPUT.
        exx_info.sync_from_global();
        // populate ABFs/JLE file lists from UnitCell; keep in sync with Exx_NAO::init
        exx_info.info_ri.files_abfs = ucell.abfs_orbital_files;
        exx_info.info_opt_abfs.files_abfs = ucell.abfs_orbital_files;
        exx_info.info_opt_abfs.files_jles = ucell.jle_orbital_files;
        this->exx_lri = std::make_shared<Exx_LRI<T>>(exx_info.info_ri);
        this->exx_lri->init(MPI_COMM_WORLD, ucell, this->kv, orb);
        this->exx_lri->cal_exx_ions(ucell, this->inp_->out_ri_cv);
    }
    // else
#endif
    // ModuleBase::Ylm::set_coefficients() is deprecated
    this->lr_.bind(this->make_refs());
}

template class ModuleESolver::ESolver_LR<double, double>;
template class ModuleESolver::ESolver_LR<std::complex<double>, double>;
