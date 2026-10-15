#include "relax_driver.h"
#include "relax_history.h"
#include "relax_stru_io.h"
#include "socket_driver.h"
#include "source_base/formatter.h"
#include "source_base/global_file.h"
#include "source_io/module_json/output_info.h"
#include "source_io/module_output/output_log.h"
#include "source_io/module_output/print_info.h"
#include "source_base/module_out/read_exit_file.h"
#include "source_io/module_parameter/parameter.h"

void Relax_Driver::relax_driver(
        ModuleESolver::ESolver* p_esolver,
        UnitCell& ucell,
        const Input_para& inp,
        std::ofstream& ofs_running)
{
    ModuleBase::TITLE("Relax_Driver", "relax_driver");
    ModuleBase::timer::start("Relax_Driver", "relax_driver");

    // Cache global configuration once so downstream helpers do not each reach
    // into PARAM / GlobalV.
    out_dir_ = PARAM.globalv.global_out_dir;
    deepks_setorb_ = PARAM.globalv.deepks_setorb;
    my_rank_ = GlobalV::MY_RANK;

    if (inp.socket_driver)
    {
        Socket_Driver socket_driver;
        socket_driver.socket_driver(p_esolver, ucell, inp, ofs_running);
        ModuleBase::timer::end("Relax_Driver", "relax_driver");
        return;
    }

    this->init_relax(ucell.nat, inp);

    // steps[0]: istep (main iteration step)
    // steps[1]: force_step
    // steps[2]: stress_step
    std::vector<int> steps = {0, 1, 1};

    // Main iteration loop for relaxation calculations
    // For scf/nscf calculations, relax_step returns true immediately,
    // so the loop exits after one iteration
    double etot = 0.0;
    ModuleBase::matrix stress(3, 3);
    ModuleBase::matrix force(ucell.nat, 3);

    // Track whether the current geometry has been evaluated by esolve().
    // After relax_step() proposes a new geometry, it is not evaluated until
    // the next esolve() call; if we exit the loop early, force/stress are stale.
    bool geometry_evaluated = false;

    while (steps[0] < inp.relax_nmax)
    {
        this->iter_info(steps, inp);
        this->esolve(steps[0], p_esolver, ucell, inp, force, stress, etot);
        geometry_evaluated = true;
        this->stru_out(steps[0], ucell, inp, etot, stress, force);
        bool converged = this->relax_step(steps, p_esolver, ucell, inp, force, stress, etot, ofs_running);
        if (!converged)
        {
            geometry_evaluated = false;
        }
        this->json_out(p_esolver, ucell, inp, force, stress);

        // Check stop conditions
        if (converged)
        {
            // Relaxation converged, exit loop immediately
            break;
        }
        else if (ModuleIO::read_exit_file(my_rank_, "EXIT", ofs_running))
        {
            // EXIT file detected, exit loop
            break;
        }

        ++steps[0];
    }

    this->final_out(steps[0], ucell, inp, etot, stress, force, geometry_evaluated, ofs_running);

    ModuleBase::timer::end("Relax_Driver", "relax_driver");
    return;
}

void Relax_Driver::init_relax(const int nat, const Input_para& inp)
{
    if (inp.calculation == "relax" || inp.calculation == "cell-relax")
    {
        if (!inp.uses_simultaneous_relaxation())
        {
            this->rl_old.init_relax(nat, inp);
        }
        else
        {
            this->rl.init_relax(nat, inp);
        }
    }
}

void Relax_Driver::iter_info(const std::vector<int>& steps, const Input_para& inp)
{
    if (inp.out_level == "ie"
            && (inp.calculation == "relax"
                || inp.calculation == "cell-relax"
                || inp.calculation == "scf"
                || inp.calculation == "nscf")
            && (inp.esolver_type != "lr"))
    {
        ModuleIO::print_screen(steps[2], steps[1], steps[0]+1);
    }

#ifdef __JSON
    // ks-lr runs an embedded KS calculation in before_all_runners(), which
    // already starts the first output record.
    if (inp.esolver_type != "ks-lr" || steps[0] != 0)
    {
        Json::init_output_array_obj();
    }
#endif
}

void Relax_Driver::esolve(const int istep,
        ModuleESolver::ESolver* p_esolver,
        UnitCell& ucell,
        const Input_para& inp,
        ModuleBase::matrix& force,
        ModuleBase::matrix& stress,
        double& etot)
{
    p_esolver->runner(ucell, istep);

    etot = p_esolver->cal_energy();

    if (inp.cal_force)
    {
        p_esolver->cal_force(ucell, force);
    }

    if (inp.cal_stress)
    {
        p_esolver->cal_stress(ucell, stress);
    }
}

bool Relax_Driver::relax_step(std::vector<int>& steps,
        ModuleESolver::ESolver* p_esolver,
        UnitCell& ucell,
        const Input_para& inp,
        const ModuleBase::matrix& force,
        const ModuleBase::matrix& stress,
        const double etot,
        std::ofstream& ofs_running)
{
    // Guard: For non-relaxation calculations (scf, nscf, etc.), return true immediately
    // to ensure the main loop exits after one iteration. This provides robustness
    // even if relax_nmax is set to a large value.
    if (inp.calculation != "relax" && inp.calculation != "cell-relax")
    {
        return true;
    }

    bool converged = false;

    if (inp.uses_simultaneous_relaxation())
    {
        converged = this->rl.relax_step(ucell, force, stress, etot, ofs_running);
        // stress step +1
        steps[2]++;
        // fix force step to 1
        steps[1] = 1;
    }
    else
    {
        converged = this->rl_old.relax_step(steps[0]+1, etot, ucell, force,
            stress, steps[1], steps[2], ofs_running);
    }

    ModuleIO::output_after_relax(converged, p_esolver->conv_esolver, ofs_running);

    return converged;
}

void Relax_Driver::stru_out(const int istep, UnitCell& ucell, const Input_para& inp, const double etot, const ModuleBase::matrix& stress, const ModuleBase::matrix& force)
{
    // out_stru is effective for scf/nscf/relax/cell-relax (md writes STRU_MD_* via md_restartfreq)
    if (inp.calculation != "relax" && inp.calculation != "cell-relax"
        && inp.calculation != "scf" && inp.calculation != "nscf")
    {
        return;
    }

    const bool is_relax = (inp.calculation == "relax" || inp.calculation == "cell-relax");

    // out_stru: -1 no output, 0 final only, 1 STRU format, 2 CIF format
    // For -1 and 0, no per-step structure output
    if (inp.out_stru <= 0)
    {
        return;
    }

    // stru_out is called right after esolve(), so the geometry has been
    // evaluated and forces/stress are consistent with the written structure.
    const bool geometry_evaluated = true;
    const std::string header = relax_stru_io::build_stru_header(istep, etot, stress, inp, false, geometry_evaluated);
    const bool need_orb = relax_stru_io::need_orbital(inp);
    const bool freq_ok = (inp.out_freq_ion > 0 && istep % inp.out_freq_ion == 0);

    // STRU_NOW: overwrite each step (for out_stru 1 and 2)
    // For scf/nscf the structure is identical to STRU_FINAL; only STRU_FINAL
    // is written in final_out() to avoid a duplicate file.
    if (is_relax)
    {
        const std::string now_file = out_dir_ + (inp.out_stru == 1 ? "STRU_NOW" : "STRU_NOW.cif");
        relax_stru_io::write_stru(ucell, inp, now_file, header, force, need_orb,
                                  deepks_setorb_, my_rank_, inp.cal_force);
    }

    // Numbered files per out_freq_ion: only meaningful for relaxation calculations
    if (is_relax && freq_ok)
    {
        const std::string step_file = out_dir_ + "STRU" + std::to_string(istep + 1)
                                      + (inp.out_stru == 1 ? "" : ".cif");
        relax_stru_io::write_stru(ucell, inp, step_file, header, force, need_orb,
                                  deepks_setorb_, my_rank_, inp.cal_force);
    }
}

void Relax_Driver::json_out(ModuleESolver::ESolver* p_esolver, UnitCell& ucell, const Input_para& inp, const ModuleBase::matrix& force, const ModuleBase::matrix& stress)
{
#ifdef __JSON
    Json::add_output_energy(p_esolver->cal_energy() * ModuleBase::Ry_to_eV);

    double unit_transform = ModuleBase::RYDBERG_SI / pow(ModuleBase::BOHR_RADIUS_SI, 3) * 1.0e-8;
    double fac = ModuleBase::Ry_to_eV / 0.529177;
    Json::add_output_cell_coo_stress_force(ucell,
                                           force,
                                           fac,
                                           stress,
                                           unit_transform,
                                           inp.cal_force,
                                           inp.cal_stress);
#endif
}

void Relax_Driver::final_out(const int istep,
                             UnitCell& ucell,
                             const Input_para& inp,
                             const double etot,
                             const ModuleBase::matrix& stress,
                             const ModuleBase::matrix& force,
                             const bool geometry_evaluated,
                             std::ofstream& ofs_running)
{
    // Structure final output is effective for scf/nscf/relax/cell-relax;
    // relax-specific screen messages remain guarded below.
    const bool is_relax = (inp.calculation == "relax" || inp.calculation == "cell-relax");
    const bool stru_effective = is_relax || inp.calculation == "scf" || inp.calculation == "nscf";

    // out_stru: 0 no output, 1 STRU format, 2 CIF format
    // 1: write STRU_FINAL; 2: write STRU_FINAL.cif
    if (stru_effective && (inp.out_stru == 1 || inp.out_stru == 2))
    {
        const std::string header = relax_stru_io::build_stru_header(istep, etot, stress, inp, true, geometry_evaluated);
        const bool need_orb = relax_stru_io::need_orbital(inp);
        const std::string final_file = out_dir_ + (inp.out_stru == 1 ? "STRU_FINAL" : "STRU_FINAL.cif");
        // Only write forces when they were actually computed and belong
        // to the geometry being written.
        const bool write_force = inp.cal_force && geometry_evaluated;
        relax_stru_io::write_stru(ucell, inp, final_file, header, force, need_orb,
                                  deepks_setorb_, my_rank_, write_force);
    }

    // relax_nmax == 0 is a valid dry-run mode: no relaxation step was ever
    // taken, so it must not be reported as a failed relaxation.
    if (istep == inp.relax_nmax && inp.relax_nmax > 0)
    {
        if (is_relax)
        {
            print_not_converged_summary(inp, ofs_running);

            std::cout << "\n ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
            std::cout << " Geometry relaxation stops here due to reaching the maximum      " << std::endl;
            std::cout << " relaxation steps. More steps are needed to converge the results " << std::endl;
            std::cout << " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
        }
    }
    else if (inp.relax_nmax > 0)
    {
        if (is_relax)
        {
            std::cout << "\n ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
            std::cout << " Geometry relaxation thresholds are reached within " << istep << " steps." << std::endl;
            std::cout << " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
        }
    }

    if (is_relax && inp.relax_nmax == 0)
    {
        print_dry_run_message();
    }
}

void Relax_Driver::print_not_converged_summary(const Input_para& inp, std::ofstream& ofs_running) const
{
    // Unified not-converged summary for both relaxation paths so ASE and
    // users can read the method, the actual number of steps taken, and
    // the per-step force / stress history from the running log.
    const std::vector<double>& force_hist = inp.uses_simultaneous_relaxation()
        ? rl.get_max_force_history()
        : rl_old.get_max_force_history();
    const std::vector<double>& stress_hist = inp.uses_simultaneous_relaxation()
        ? rl.get_max_stress_history()
        : rl_old.get_max_stress_history();
    const int ionic_steps = static_cast<int>(force_hist.size());
    const std::string method = inp.relax_method.empty() ? "unknown" : inp.relax_method[0];
    ofs_running << " Relaxation method: " << method << std::endl;
    ofs_running << " Relaxation stopped after " << ionic_steps << " ionic step(s) (relax_nmax = "
                << inp.relax_nmax << " reached)." << std::endl;
    if (!force_hist.empty())
    {
        ofs_running << " Largest force per step (eV/Angstrom):" << format_relax_history(force_hist);
    }
    if (!stress_hist.empty())
    {
        ofs_running << " Largest stress per step (kbar):" << format_relax_history(stress_hist);
    }
    ofs_running << " Relaxation is not converged after reaching relax_nmax!" << std::endl;
}

void Relax_Driver::print_dry_run_message() const
{
    std::cout << "-----------------------------------------------" << std::endl;
    std::cout << " relax_nmax = 0, DRY RUN TEST SUCCEEDS :)" << std::endl;
    std::cout << "-----------------------------------------------" << std::endl;
}
