#include "run_md.h"

#include "source_base/constants.h"
#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "source_base/parallel_cell.h"
#include "source_cell/mdcell_reader.h"
#include "source_cell/mdcell.h"
#include "source_cell/module_neighlist/domain_decomposition.h"
#include "source_io/module_parameter/parameter.h"
#include "fire.h"
#include "langevin.h"
#include "md_func.h"
#include "source_base/global_file.h"
#include "source_base/timer.h"
#include "source_io/module_output/print_info.h"
#include "msst.h"
#include "nhchain.h"
#include "verlet.h"
#include "source_cell/update_cell.h"
#include "source_cell/print_cell.h"

#include <vector>

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

void md_line(MDCell& mdcell,
             ModuleESolver::ESolver* p_esolver,
             const Parameter& param_in,
             DomainDecomposition& decomp)
{
    ModuleBase::TITLE("Run_MD", "md_line");
    ModuleBase::timer::start("Run_MD", "md_line");
    /// determine the md_type
    MD_base* mdrun = nullptr;
    if (param_in.mdp.md_type == "fire")
    {
        mdrun = new FIRE(param_in, mdcell);
    }
    else if ((param_in.mdp.md_type == "nvt" && param_in.mdp.md_thermostat == "nhc") || param_in.mdp.md_type == "npt")
    {
        mdrun = new Nose_Hoover(param_in, mdcell);
    }
    else if (param_in.mdp.md_type == "nve" || param_in.mdp.md_type == "nvt")
    {
        mdrun = new Verlet(param_in, mdcell);
    }
    else if (param_in.mdp.md_type == "langevin")
    {
        mdrun = new Langevin(param_in, mdcell);
    }
    else if (param_in.mdp.md_type == "msst")
    {
        mdrun = new MSST(param_in, mdcell);
    }
    else
    {
        ModuleBase::WARNING_QUIT("md_line", "no such md_type!");
    }

    /// md cycle, mohan update 2026-01-04, change '<=' to '<'
    while ((mdrun->step_ + mdrun->step_rst_) < param_in.mdp.md_nstep && !mdrun->stop)
    {
        if (mdrun->step_ == 0)
        {
            mdrun->setup(p_esolver, param_in.globalv.global_readin_dir, decomp);
        }
        else
        {
            // mohan add 2026-01-04
            const int stress_step = 0;
            const int force_step = 0;
            const int istep_print = mdrun->step_ + mdrun->step_rst_ + 1;
            ModuleIO::print_screen(stress_step, force_step, istep_print);
            mdrun->first_half(GlobalV::ofs_running);

            /// update force and virial due to the update of atom positions
            MD_func::force_virial(p_esolver,
                                  mdrun->step_,
                                  mdcell,
                                  decomp,
                                  mdrun->potential,
                                  param_in.inp.cal_stress,
                                  mdrun->virial,
                                  param_in.mdp.md_out_force);

            mdrun->second_half();

            MD_func::compute_stress(mdcell,
                                    param_in.inp.cal_stress,
                                    mdrun->virial,
                                    mdrun->stress);
            mdrun->t_current = MD_func::current_temp(mdrun->kinetic,
                                                     mdcell,
                                                     mdrun->frozen_freedom_);
        }

        mdrun->print_md(GlobalV::ofs_running, PARAM.inp.cal_stress);
        if (param_in.mdp.md_dumpfreq > 0
            && (mdrun->step_ + mdrun->step_rst_) % param_in.mdp.md_dumpfreq == 0)
        {
            MD_func::dump_info(mdrun->step_ + mdrun->step_rst_,
                               PARAM.globalv.global_out_dir,
                               mdcell,
                               param_in,
                               mdrun->virial);
        }

        if (param_in.mdp.md_restartfreq > 0
            && (mdrun->step_ + mdrun->step_rst_) % param_in.mdp.md_restartfreq == 0)
        {
            if (mdcell.has_backing_unitcell())
            {
                mdcell.sync_backing_unitcell();
            }
            std::stringstream file;
            file << PARAM.globalv.global_stru_dir << "STRU_MD_" << mdrun->step_ + mdrun->step_rst_;
            mdcell::print_stru_file(mdcell, mdcell.stru_meta(), file.str());
            mdrun->write_restart(PARAM.globalv.global_out_dir);
        }

        mdrun->step_++;
    }

    delete mdrun;
    ModuleBase::timer::end("Run_MD", "md_line");
    return;
}

} // namespace Run_MD
