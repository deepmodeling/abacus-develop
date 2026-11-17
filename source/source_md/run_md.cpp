#include "run_md.h"

#include "source_base/constants.h"
#include "source_base/fs_compat.h"
#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "source_base/parallel_cell.h"
#include "source_cell/mdcell_reader.h"
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

void prepare_mdcell(MDCell& mdcell, const Parameter& param_in, DomainDecomposition& decomp)
{
    const Input_para& input = param_in.inp;
    std::vector<int> effective_replicate = input.cell_replica;
    if (input.mdp.md_restart)
    {
        effective_replicate = {1, 1, 1};
    }

    const ModuleBase::CommunicationDomain comm_domain = ModuleBase::world_comm_domain();
    mdcell = MDCellReader::read_stru(param_in.globalv.global_in_stru,
                                     effective_replicate,
                                     input.mdp.md_neighbor_skin / ModuleBase::BOHR_TO_A,
                                     comm_domain,
                                     decomp);
    GlobalV::ofs_running << std::endl;
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "TOTAL ATOM NUMBER", mdcell.nat());
    GlobalV::ofs_running << std::endl;
}

void prepare_mdcell(MDCell& mdcell, UnitCell& ucell, DomainDecomposition& decomp)
{
    const ModuleBase::CommunicationDomain comm_domain = ModuleBase::world_comm_domain();
    decomp.init(comm_domain, ucell.latvec, ucell.lat0, 0.0, 0.0);
    const std::vector<LocalAtom> owned_atoms = decomp.split_owned_atoms_from_ucell(ucell);
    std::vector<std::string> type_labels;
    std::vector<double> type_masses;
    std::vector<std::int64_t> type_atom_counts;
    for (int it = 0; it < ucell.ntype; ++it)
    {
        type_labels.push_back(ucell.atoms[it].label);
        type_masses.push_back(ucell.atoms[it].mass);
        type_atom_counts.push_back(ucell.atoms[it].na);
    }
    mdcell.initialize_from_owned_atoms(ucell.latvec,
                                       ucell.GT,
                                       ucell.lat0,
                                       ucell.omega,
                                       ucell.nat,
                                       owned_atoms,
                                       type_labels,
                                       type_masses,
                                       type_atom_counts,
                                       0.0,
                                       comm_domain);
    mdcell.set_backing_unitcell(ucell);
    mdcell.mutable_stru_meta() = unitcell::make_stru_meta(ucell);
}

namespace
{
/// @brief collect the absolute Cartesian positions (in Bohr) of the atoms
/// owned by this rank; the cell stores cart in units of the lattice constant
void collect_plumed_positions(const MDCell& mdcell, std::vector<ModuleBase::Vector3<double>>& positions)
{
    const double lat0 = mdcell.lat0();
    const std::vector<LocalAtom>& atoms = mdcell.owned_atoms();
    positions.resize(atoms.size());
    for (std::size_t i = 0; i < atoms.size(); ++i)
    {
        positions[i] = atoms[i].cart * lat0;
    }
}

/// @brief gather the atomic forces (in Hartree/Bohr)
void collect_plumed_forces(const MDCell& mdcell, std::vector<ModuleBase::Vector3<double>>& forces)
{
    const std::vector<LocalAtom>& atoms = mdcell.owned_atoms();
    forces.resize(atoms.size());
    for (std::size_t i = 0; i < atoms.size(); ++i)
    {
        forces[i] = atoms[i].force;
    }
}

/// @brief write the (PLUMED-biased) forces back into the cell
void apply_plumed_forces(MDCell& mdcell, const std::vector<ModuleBase::Vector3<double>>& forces)
{
    std::vector<LocalAtom>& atoms = mdcell.owned_atoms();
    for (std::size_t i = 0; i < atoms.size(); ++i)
    {
        atoms[i].force = forces[i];
    }
}

/// @brief 3x3 cell matrix (in Bohr, row-major) expected by PLUMED
void collect_plumed_cell(const MDCell& mdcell, double* cell)
{
    const double lat0 = mdcell.lat0();
    const double* latvec = &mdcell.latvec().e11; // unitless lattice vectors
    for (int i = 0; i < 9; ++i)
    {
        cell[i] = latvec[i] * lat0;
    }
}
} // namespace

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
        collect_plumed_positions(mdcell, plumed_positions);
        collect_plumed_forces(mdcell, plumed_forces);
        collect_plumed_cell(mdcell, plumed_cell);
        plumed.move(istep,
                    plumed_positions.data(),
                    plumed_forces.data(),
                    plumed_cell,
                    mdrun->potential,
                    param_in.inp.cal_stress ? mdrun->virial.c : nullptr);
        apply_plumed_forces(mdcell, plumed_forces);
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

    /// Write one STRU frame for `current_step`. Cadence frames and the final
    /// frame share the same naming/append rules, so both call sites use this
    /// lambda. `current_step == 0` marks the first frame of a fresh run, which
    /// truncates the fixed-name file instead of appending to it.
    auto write_stru_frame = [&](const int current_step) {
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
        const bool is_first_step_of_fresh_run = (current_step == 0);
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
        const bool do_output_step = is_cadence_step(current_step, eff_output_freq);
        const bool do_restart_step = is_cadence_step(current_step, restart_freq);
        const bool write_any_frame = is_frame_write_step(current_step,
                                                         param_in.inp.out_freq_ion,
                                                         restart_freq,
                                                         param_in.inp.out_stru > 0);

        /// A checkpoint STRU frame must exist on every restart step so that the
        /// state written to Restart_md.txt (coordinates + thermostat state +
        /// step number) can actually be resumed from. This frame is paired
        /// with Restart_md.txt on the md_restartfreq cadence and is independent
        /// of both out_freq_ion and out_stru: out_stru=0 silences the regular
        /// per-step output only, not the restart-ready frame. When the same
        /// step fires both the output cadence and the restart cadence we write
        /// one frame shared by both paths, so the cumulative STRU is not
        /// duplicated in out_app_flag=1 mode. See is_frame_write_step for the
        /// formal predicate behind write_any_frame.
        if (write_any_frame)
        {
            write_stru_frame(current_step);
        }
        last_processed_step = current_step;

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

    /// Supplementary final frame. When the last step is off both the output
    /// and the restart cadences the loop above never writes its STRU, so write
    /// it here to guarantee the final geometry is persisted. Steps already
    /// covered by a cadence are not written twice. out_stru=0 keeps the
    /// existing contract of silencing regular per-step geometry output, so the
    /// final frame is skipped as well.
    if (last_processed_step >= 0 && param_in.inp.out_stru > 0
        && !is_frame_write_step(last_processed_step,
                                param_in.inp.out_freq_ion,
                                restart_freq,
                                true))
    {
        write_stru_frame(last_processed_step);
    }

    /// flush the COLVAR and other PLUMED output files
    plumed.finalize();

    delete mdrun;
    ModuleBase::timer::end("Run_MD", "md_line");
    return;
}

} // namespace Run_MD
