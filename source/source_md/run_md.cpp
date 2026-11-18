#include "run_md.h"

#include "source_base/fs_compat.h"
#include "source_base/global_variable.h"
#include "source_cell/mdcell.h"
#include "source_cell/dom_decomp.h"
#include "source_io/module_parameter/parameter.h"
#include "fire.h"
#include "langevin.h"
#include "md_func.h"
#include "source_base/global_file.h"
#include "source_base/timer.h"
#ifdef __JSON
#include "source_io/module_json/output_info.h"
#endif
#include "source_io/module_output/print_info.h"
#include "msst.h"
#include "nhchain.h"
#include "verlet.h"
#include "source_cell/update_cell.h"
#include "source_cell/print_cell.h"
#include "source_md/module_plumed/plumed_interface.h"

#include <vector>
#include <sys/stat.h>
#ifdef __MPI
#include <mpi.h>
#endif

namespace Run_MD
{

void md_line(MDCell& mdcell,
             ModuleESolver::ESolver* p_esolver,
             const Parameter& param_in,
             DomainDecomposition& decomp)
{
    ModuleBase::TITLE("Run_MD", "md_line");
    ModuleBase::timer::start("Run_MD", "md_line");
    /// determine the md_type
    MD_base* mdrun = nullptr;
    /// the integrators take the values they use, not the whole Parameter
    const MD_para& mdp = param_in.mdp;
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

    /// interface to PLUMED: no-op unless `plumed` is set in INPUT and ABACUS
    /// was compiled with PLUMED support (-DPLUMED_ROOT=...)
    ModulePlumed::PlumedInterface plumed;
    std::vector<ModuleBase::Vector3<double>> plumed_positions;
    std::vector<ModuleBase::Vector3<double>> plumed_forces;
    double plumed_cell[9];

    auto run_plumed_step = [&](const int istep) {
        std::vector<LocalAtom>& atoms = mdcell.owned_atoms();
        const double lat0 = mdcell.lat0();
        plumed_positions.resize(atoms.size());
        plumed_forces.resize(atoms.size());
        for (std::size_t i = 0; i < atoms.size(); ++i)
        {
            plumed_positions[i] = atoms[i].cart * lat0; // Bohr
            plumed_forces[i] = atoms[i].force; // Hartree/Bohr
        }

        const double* latvec = &mdcell.latvec().e11;
        for (int i = 0; i < 9; ++i)
        {
            plumed_cell[i] = latvec[i] * lat0; // row-major, Bohr
        }

        const double* virial = param_in.inp.cal_stress ? mdrun->virial.c : nullptr;
        plumed.move(istep,
                    plumed_positions.data(),
                    plumed_forces.data(),
                    plumed_cell,
                    mdrun->potential,
                    virial);
        for (std::size_t i = 0; i < atoms.size(); ++i)
        {
            atoms[i].force = plumed_forces[i];
        }
    };

    /// Effective ionic-step output frequency for STRU and other per-step
    /// output files. Prefer out_freq_ion; fall back to md_restartfreq for
    /// backward compatibility when out_freq_ion is not set (0).
    const int eff_output_freq = (param_in.inp.out_freq_ion > 0)
                                    ? param_in.inp.out_freq_ion
                                    : param_in.mdp.md_restartfreq;

    /// Restart-file cadence is independent of the output cadence so that
    /// out_freq_ion and md_restartfreq can be set to different values.
    const int restart_freq = param_in.mdp.md_restartfreq;

    /// out_app_flag == 1 writes one fixed-name "STRU" file. A previous
    /// release kept per-step STRU_MD_{step} files inside a STRU/
    /// subdirectory of the output directory (the mode-2 directory has
    /// since been renamed to STRUC/), so an MD restarted from those
    /// files still sees that directory. Writing the fixed name on top of a
    /// directory aborts the run, therefore the legacy numbered layout is
    /// kept whenever such a directory is present.
    const bool append_stru_allowed = (param_in.inp.out_app_flag == 1);
    bool has_legacy_stru_dir = false;
    std::string stru_out_dir = PARAM.globalv.global_stru_dir;
    if (append_stru_allowed)
    {
        struct stat stru_path;
        if (stat((stru_out_dir + "STRU").c_str(), &stru_path) == 0 && S_ISDIR(stru_path.st_mode))
        {
            has_legacy_stru_dir = true;
            stru_out_dir += "STRU/";
        }
    }

    /// md cycle, mohan update 2026-01-04, change '<=' to '<'
    while ((mdrun->step_ + mdrun->step_rst_) < param_in.mdp.md_nstep && !mdrun->stop)
    {
        const int md_step = mdrun->step_ + mdrun->step_rst_;
        bool print_thermo = param_in.mdp.md_thermo_freq > 0
                            && md_step % param_in.mdp.md_thermo_freq == 0;
#ifdef __JSON
        // JSON output currently follows the UnitCell-backed electronic-structure path.
        // Start one output record before the solver appends SCF information for this MD step.
        if (mdcell.has_backing_unitcell())
        {
            Json::init_output_array_obj();
        }
#endif

        if (mdrun->step_ == 0)
        {
            mdrun->setup(p_esolver, param_in.globalv.global_readin_dir, decomp);

            // setup() may load step_rst_ from Restart_md.txt.
            const int restart_step = mdrun->step_ + mdrun->step_rst_;
            print_thermo = param_in.mdp.md_thermo_freq > 0
                           && restart_step % param_in.mdp.md_thermo_freq == 0;

            if (param_in.mdp.plumed)
            {
                int nproc = 1;
#ifdef __MPI
                MPI_Comm_size(MPI_COMM_WORLD, &nproc);
#endif
                if (nproc > 1)
                {
                    ModuleBase::WARNING_QUIT("md_line",
                                             "the PLUMED interface currently supports a single MPI rank only; "
                                             "run the MD with one process (OpenMP threads are fine)");
                }
                std::vector<double> plumed_masses(mdcell.owned_atoms().size());
                for (std::size_t i = 0; i < mdcell.owned_atoms().size(); ++i)
                {
                    plumed_masses[i] = mdcell.owned_atoms()[i].mass;
                }
                plumed.init(param_in.mdp.plumed_file,
                            static_cast<int>(mdcell.owned_atoms().size()),
                            mdrun->get_md_dt(),
                            plumed_masses.data());
                /// the initial configuration already carries positions and
                /// forces, so PLUMED can be called right away
                run_plumed_step(mdrun->step_ + mdrun->step_rst_);
            }
        }
        else
        {
            // mohan add 2026-01-04
            const int stress_step = 0;
            const int force_step = 0;
            const int istep_print = md_step + 1;
            if (print_thermo)
            {
                ModuleIO::print_screen(stress_step, force_step, istep_print);
            }
            mdrun->first_half(GlobalV::ofs_running);

            /// update force and virial due to the update of atom positions.
            /// The electronic solver interprets the step index as the local
            /// step of this run: istep == 0/1 enables its initialization
            /// branches (init_chg dm/hr, init_wfc file, PW need_subspace),
            /// so pass the local step; the global step is used only by the
            /// MD-side output writers below (print_md, dump_info, STRU
            /// naming, write_restart, PLUMED).
            /// first_half has already drifted the positions in this run, so
            /// the charge-extrapolation history may advance.
            const bool positions_updated = true;
            MD_func::force_virial(p_esolver,
                                  mdrun->step_,
                                  positions_updated,
                                  mdcell,
                                  decomp,
                                  mdrun->potential,
                                  param_in.inp.cal_stress,
                                  mdrun->virial);

            if (plumed.active())
            {
                /// let PLUMED compute the collective variables and add the
                /// biasing forces; they are used by the second half of the
                /// current step and by the first half of the next one
                run_plumed_step(mdrun->step_ + mdrun->step_rst_);
            }

            mdrun->second_half();

            MD_func::compute_stress(mdcell,
                                    param_in.inp.cal_stress,
                                    mdrun->virial,
                                    mdrun->stress);
            mdrun->t_current = MD_func::current_temp(mdrun->kinetic,
                                                     mdcell,
                                                     mdrun->frozen_freedom_);
        }

        if (print_thermo)
        {
            mdrun->print_md(GlobalV::ofs_running, param_in.inp.cal_stress);
        }
        if (param_in.mdp.md_dumpfreq > 0
            && (mdrun->step_ + mdrun->step_rst_) % param_in.mdp.md_dumpfreq == 0)
        {
            MD_func::dump_info(mdrun->step_ + mdrun->step_rst_,
                               PARAM.globalv.global_out_dir,
                               mdcell,
                               mdp,
                               cal_stress,
                               mdrun->virial);
        }

        const int current_step = mdrun->step_ + mdrun->step_rst_;
        const bool do_output_step = eff_output_freq > 0 && current_step % eff_output_freq == 0;
        const bool do_restart_step = restart_freq > 0 && current_step % restart_freq == 0;
        const bool write_any_frame = (do_output_step && param_in.inp.out_stru > 0) || do_restart_step;

        /// A checkpoint STRU frame must exist on every restart step so that the
        /// state written to Restart_md.txt (coordinates + thermostat state +
        /// step number) can actually be resumed from. This frame is paired
        /// with Restart_md.txt on the md_restartfreq cadence and is independent
        /// of both out_freq_ion and out_stru: out_stru=0 silences the regular
        /// per-step output only, not the restart-ready frame. When the same
        /// step fires both the output cadence and the restart cadence we write
        /// one frame shared by both paths, so the cumulative STRU is not
        /// duplicated in out_app_flag=1 mode.
        if (write_any_frame)
        {
            if (mdcell.has_backing_unitcell())
            {
                const bool positions_updated = false;
                mdcell.sync_backing_unitcell(positions_updated);
            }
            // Mode-2 STRUC/ subdirectory is created lazily at write time.
            ModuleBase::ensure_directory(stru_out_dir);
            const int stru_step = current_step + 1;
            // File-name selection and append/truncate are two independent
            // decisions. The file name is fixed to "STRU" whenever mode 1
            // accumulation is allowed (and no legacy STRU/ directory is in
            // the way); modes 0/2 keep the numbered STRU{stru_step} files.
            // Whether to append depends on whether this is the very first
            // step of a fresh (non-restart) run: a fresh first frame
            // truncates the fixed-name file so it does not keep frames
            // from a previous run, while later steps and restarted runs
            // append one frame each, introduced by a one-based
            // "# MD step N" comment line.
            const bool use_fixed_name = append_stru_allowed && !has_legacy_stru_dir;
            const bool is_first_step_of_fresh_run = (mdrun->step_ == 0) && (mdrun->step_rst_ == 0);
            const bool append_stru = use_fixed_name && !is_first_step_of_fresh_run;
            const std::string frame_comment = use_fixed_name
                ? "# MD step " + std::to_string(stru_step) + "\n\n"
                : "";
            std::stringstream file;
            file << stru_out_dir << "STRU";
            if (!use_fixed_name)
            {
                file << stru_step;
            }
            mdcell::print_stru_file(mdcell, mdcell.stru_meta(), file.str(), frame_comment, append_stru);
        }

        /// Restart files follow md_restartfreq independently of the output
        /// cadence, so setting out_freq_ion does not suppress or re-aim
        /// md_restartfreq. A matching checkpoint STRU frame was written in
        /// the block above whenever do_restart_step was true.
        if (do_restart_step)
        {
            mdrun->write_restart(PARAM.globalv.global_out_dir);
        }

        mdrun->step_++;
    }

    /// flush the COLVAR and other PLUMED output files
    plumed.finalize();

    delete mdrun;
    ModuleBase::timer::end("Run_MD", "md_line");
    return;
}

} // namespace Run_MD
