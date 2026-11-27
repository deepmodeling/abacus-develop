#include "run_md.h"

#include "fire.h"
#include "langevin.h"
#include "md_func.h"
#include "source_base/fs_compat.h"
#include "source_base/global_file.h"
#include "source_base/global_variable.h"
#include "source_base/timer.h"
#include "source_cell/dom_decomp.h"
#include "source_cell/mdcell.h"
#include "source_cell/mdcell_cif_io.h"
#include "source_io/module_parameter/parameter.h"
#ifdef __JSON
#include "source_io/module_json/output_info.h"
#endif
#include "msst.h"
#include "nhchain.h"
#include "source_cell/print_cell.h"
#include "source_cell/update_cell.h"
#include "source_io/module_output/frame_append.h"
#include "source_io/module_output/print_info.h"
#include "source_md/module_plumed/plumed_interface.h"
#include "verlet.h"

#include <sys/stat.h>
#include <unistd.h>
#include <vector>
#ifdef __MPI
#include <mpi.h>
#endif

namespace Run_MD
{

void md_line(
    MDCell &mdcell,
    ModuleESolver::ESolver *p_esolver,
    const Parameter &param_in,
    DomainDecomposition &decomp)
{
    ModuleBase::TITLE("Run_MD", "md_line");
    ModuleBase::timer::start("Run_MD", "md_line");
    MD_base *mdrun = nullptr;
    const MD_para &mdp = param_in.mdp;
    const bool cal_stress = param_in.inp.cal_stress;
    const bool init_vel = param_in.inp.init_vel;
    const int my_rank = param_in.globalv.myrank;
    if (param_in.mdp.md_type == "fire")
    {
        mdrun = new FIRE(mdp, cal_stress, init_vel, my_rank, param_in.inp.force_thr, mdcell);
    }
    else if ((param_in.mdp.md_type == "nvt" && param_in.mdp.md_thermostat == "nhc") || param_in.mdp.md_type == "npt")
    {
        mdrun = new Nose_Hoover(mdp, cal_stress, init_vel, my_rank, mdcell);
    }
    else if (param_in.mdp.md_type == "nve" || param_in.mdp.md_type == "nvt")
    {
        mdrun = new Verlet(mdp, cal_stress, init_vel, my_rank, mdcell);
    }
    else if (param_in.mdp.md_type == "langevin")
    {
        mdrun = new Langevin(mdp, cal_stress, init_vel, my_rank, mdcell);
    }
    else if (param_in.mdp.md_type == "msst")
    {
        mdrun = new MSST(mdp, cal_stress, init_vel, my_rank, mdcell);
    }
    else
    {
        ModuleBase::WARNING_QUIT("md_line", "no such md_type!");
    }

    ModulePlumed::PlumedInterface plumed;
    std::vector<ModuleBase::Vector3<double>> plumed_positions;
    std::vector<ModuleBase::Vector3<double>> plumed_forces;
    double plumed_cell[9];

    auto run_plumed_step = [&](const int istep)
    {
        std::vector<LocalAtom> &atoms = mdcell.owned_atoms();
        const double lat0 = mdcell.lat0();
        plumed_positions.resize(atoms.size());
        plumed_forces.resize(atoms.size());
        for (std::size_t i = 0; i < atoms.size(); ++i)
        {
            plumed_positions[i] = atoms[i].cart * lat0; // Bohr
            plumed_forces[i] = atoms[i].force;          // Hartree/Bohr
        }

        const double *latvec = &mdcell.latvec().e11;
        for (int i = 0; i < 9; ++i)
        {
            plumed_cell[i] = latvec[i] * lat0; // row-major, Bohr
        }

        const double *virial = param_in.inp.cal_stress ? mdrun->virial.c : nullptr;
        plumed.move(istep, plumed_positions.data(), plumed_forces.data(), plumed_cell, mdrun->potential, virial);
        for (std::size_t i = 0; i < atoms.size(); ++i)
        {
            atoms[i].force = plumed_forces[i];
        }
    };

    const int output_freq = param_in.inp.out_freq_ion;
    const int restart_freq = param_in.mdp.md_restartfreq;
    const std::string stru_out_dir = param_in.globalv.global_stru_dir;
    const std::string restart_stru = param_in.globalv.global_out_dir + "STRU_restart";

    auto log_written_file = [&](const std::string &filename)
    {
        if (my_rank != 0)
        {
            return;
        }
        const std::string full_path = ModuleBase::absolute_path_for_log(filename);
        GlobalV::ofs_running << "Write file: " << full_path << '\n';
    };

    auto write_stru_frame = [&](const int current_step, const std::string &filename, const bool append,
                                const bool include_force)
    {
        if (mdcell.has_backing_unitcell())
        {
            const bool positions_updated = false;
            mdcell.sync_backing_unitcell(positions_updated);
        }
        const int stru_step = current_step + 1;
        const std::string frame_comment = "# MD step " + std::to_string(stru_step) + "\n\n";
        mdcell::print_stru_file(mdcell, mdcell.stru_meta(), filename, frame_comment,
                                append, include_force);
        log_written_file(filename);
    };

    /// md cycle, mohan update 2026-01-04, change '<=' to '<'
    int last_processed_step = -1;
    while ((mdrun->step_ + mdrun->step_rst_) < param_in.mdp.md_nstep && !mdrun->stop)
    {
        const int md_step = mdrun->step_ + mdrun->step_rst_;
        /// Mark the final ionic step on the cell passed to the solver. The
        /// output module then writes the final snapshot even when this step is
        /// off the out_freq_ion cadence or its SCF did not converge.
        const bool is_final_step = is_last_step(md_step, param_in.mdp.md_nstep);
        if (mdcell.has_backing_unitcell())
        {
            mdcell.backing_unitcell().last_ion_step = is_final_step;
        }
        bool print_thermo = param_in.mdp.out_thermo_freq > 0 && md_step % param_in.mdp.out_thermo_freq == 0;
#ifdef __JSON
        // Start the UnitCell-backed JSON record before the solver writes SCF data.
        if (mdcell.has_backing_unitcell())
        {
            Json::init_output_array_obj();
        }
#endif

        if (mdrun->step_ == 0)
        {
            mdrun->setup(p_esolver, param_in.globalv.global_readin_dir, decomp);

            // setup() may change step_rst_ when resuming a run.
            const int restart_step = mdrun->step_ + mdrun->step_rst_;
            print_thermo = param_in.mdp.out_thermo_freq > 0 && restart_step % param_in.mdp.out_thermo_freq == 0;

            if (param_in.mdp.plumed)
            {
                int nproc = 1;
#ifdef __MPI
                MPI_Comm_size(MPI_COMM_WORLD, &nproc);
#endif
                if (nproc > 1)
                {
                    ModuleBase::WARNING_QUIT(
                        "md_line", "the PLUMED interface currently supports a single MPI rank only; "
                                   "run the MD with one process (OpenMP threads are fine)");
                }
                std::vector<double> plumed_masses(mdcell.owned_atoms().size());
                for (std::size_t i = 0; i < mdcell.owned_atoms().size(); ++i)
                {
                    plumed_masses[i] = mdcell.owned_atoms()[i].mass;
                }
                plumed.init(
                    param_in.mdp.plumed_file, static_cast<int>(mdcell.owned_atoms().size()), mdrun->get_md_dt(),
                    plumed_masses.data());
                // The initial state already has positions and forces.
                run_plumed_step(mdrun->step_ + mdrun->step_rst_);
            }
        }
        else
        {
            const int stress_step = 0;
            const int force_step = 0;
            const int istep_print = md_step + 1;
            if (print_thermo)
            {
                ModuleIO::print_screen(stress_step, force_step, istep_print);
            }
            mdrun->first_half(GlobalV::ofs_running);

            // The solver needs the local step for initialization; output uses the global step.
            // first_half has moved the atoms, so charge extrapolation may advance.
            const bool positions_updated = true;
            MD_func::force_virial(
                p_esolver, mdrun->step_, positions_updated, mdcell, decomp, mdrun->potential, param_in.inp.cal_stress,
                mdrun->virial);

            if (plumed.active())
            {
                // PLUMED adds forces before the integrator's second half.
                run_plumed_step(mdrun->step_ + mdrun->step_rst_);
            }

            mdrun->second_half();

            MD_func::compute_stress(mdcell, param_in.inp.cal_stress, mdrun->virial, mdrun->stress);
            mdrun->t_current = MD_func::current_temp(mdrun->kinetic, mdcell, mdrun->frozen_freedom_);
        }

        if (print_thermo)
        {
            mdrun->print_md(GlobalV::ofs_running, param_in.inp.cal_stress);
        }
        const int current_step = mdrun->step_ + mdrun->step_rst_;
        const bool ordinary_output_enabled = param_in.inp.out_stru > 0;
        const bool do_output_step = is_structure_output_step(current_step, output_freq, ordinary_output_enabled);
        const bool do_restart_step = is_cadence_step(current_step, restart_freq);

        if (do_output_step)
        {
            ModuleBase::ensure_directory(stru_out_dir);
            std::string periodic_dir = stru_out_dir;
            if (param_in.inp.out_stru == 1 && param_in.inp.out_app_flag != 2)
            {
                // Older MD runs used OUT.<suffix>/STRU/ as a directory.
                // Keep that directory and write the fixed-name frame inside it.
                struct stat stru_path;
                const std::string legacy_path = stru_out_dir + "STRU";
                if (stat(legacy_path.c_str(), &stru_path) == 0 && S_ISDIR(stru_path.st_mode))
                {
                    periodic_dir = legacy_path + "/";
                }
            }
            const std::string filename = ModuleIO::stru_frame_filename(
                periodic_dir, param_in.inp.out_app_flag, current_step);
            if (param_in.inp.out_stru == 2)
            {
                const int cif_step = current_step + 1;
                const std::string frame_comment = "# MD step " + std::to_string(cif_step);
                const std::string cif_filename = filename + ".cif";
                mdcell::print_cif_file(mdcell, cif_filename, frame_comment);
                log_written_file(cif_filename);
            }
            else
            {
                const bool append = ModuleIO::frame_appends(param_in.inp.out_app_flag, current_step);
                const bool include_force = true;
                write_stru_frame(current_step, filename, append, include_force);
            }
        }
        last_processed_step = current_step;

        if (do_restart_step)
        {
            const bool append = false;
            const bool include_force = false;
            write_stru_frame(current_step, restart_stru, append, include_force);
            mdrun->write_restart(param_in.globalv.global_out_dir);
            const std::string restart_metadata = param_in.globalv.global_out_dir + "Restart_md.txt";
            log_written_file(restart_metadata);
        }

        mdrun->step_++;
    }

    // Final geometry is a separate output and is written even when the last
    // step was also selected for periodic output.
    if (last_processed_step >= 0 && param_in.inp.out_stru > 0)
    {
        ModuleBase::ensure_directory(stru_out_dir);
        const std::string filename = stru_out_dir + "STRU_FINAL";
        if (param_in.inp.out_stru == 2)
        {
            const int cif_step = last_processed_step + 1;
            const std::string frame_comment = "# MD step " + std::to_string(cif_step);
            const std::string cif_filename = filename + ".cif";
            mdcell::print_cif_file(mdcell, cif_filename, frame_comment);
            log_written_file(cif_filename);
        }
        else
        {
            const bool append = false;
            const bool include_force = true;
            write_stru_frame(last_processed_step, filename, append, include_force);
        }
    }

    /// flush the COLVAR and other PLUMED output files
    plumed.finalize();

    delete mdrun;
    ModuleBase::timer::end("Run_MD", "md_line");
    return;
}

} // namespace Run_MD
