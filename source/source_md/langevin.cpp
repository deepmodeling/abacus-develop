#include "langevin.h"

#include "md_func.h"
#include "source_cell/unitcell.h"
#include "source_base/parallel_common.h"
#include "source_base/timer.h"

#ifdef __MPI
#include <mpi.h>
#endif

Langevin::Langevin(const Parameter& param_in, MDCell& mdcell_in) : MD_base(param_in, mdcell_in)
{
    /// convert to a.u. unit
    assert(ModuleBase::AU_to_FS!=0.0);

    md_damp = mdp.md_damp / ModuleBase::AU_to_FS;

    total_force.resize(static_cast<std::size_t>(mdcell.nlocal()));
}


void Langevin::setup(ModuleESolver::ESolver* p_esolver, const std::string& global_readin_dir)
{
    ModuleBase::TITLE("Langevin", "setup");
    ModuleBase::timer::start("Langevin", "setup");

    MD_base::setup(p_esolver, global_readin_dir);

    post_force();

    ModuleBase::timer::end("Langevin", "setup");
    return;
}


void Langevin::first_half(std::ofstream& ofs)
{
    ModuleBase::TITLE("Langevin", "first_half");
    ModuleBase::timer::start("Langevin", "first_half");

    for (int i = 0; i < mdcell.nlocal(); ++i)
    {
        LocalAtom& atom = mdcell.mutable_owned_atoms()[static_cast<std::size_t>(i)];
        for (int k = 0; k < 3; ++k)
        {
            if (atom.mbl[k]) atom.vel[k] += 0.5 * total_force[i][k] * md_dt / atom.mass;
        }
    }
    MD_base::update_pos();

    ModuleBase::timer::end("Langevin", "first_half");
    return;
}


void Langevin::second_half()
{
    ModuleBase::TITLE("Langevin", "second_half");
    ModuleBase::timer::start("Langevin", "second_half");

    post_force();
    for (int i = 0; i < mdcell.nlocal(); ++i)
    {
        LocalAtom& atom = mdcell.mutable_owned_atoms()[static_cast<std::size_t>(i)];
        for (int k = 0; k < 3; ++k)
        {
            if (atom.mbl[k]) atom.vel[k] += 0.5 * total_force[i][k] * md_dt / atom.mass;
        }
    }

    ModuleBase::timer::end("Langevin", "second_half");
    return;
}


void Langevin::print_md(std::ofstream& ofs, const bool& cal_stress)
{
    MD_base::print_md(ofs, cal_stress);
    return;
}


void Langevin::write_restart(const std::string& global_out_dir)
{
    MD_base::write_restart(global_out_dir);
    return;
}


void Langevin::restart(const std::string& global_readin_dir)
{
    MD_base::restart(global_readin_dir);
    return;
}


void Langevin::post_force()
{
    double t_target = MD_func::target_temp(step_ + step_rst_, mdp.md_nstep, md_tfirst, md_tlast);
    total_force.resize(static_cast<std::size_t>(mdcell.nlocal()));

    std::vector<ModuleBase::Vector3<double> > random_values(static_cast<std::size_t>(mdcell.nlocal()));
#ifdef __MPI
    if (mdcell.mpi_size() > 1)
    {
        const std::vector<int>& type_atom_counts = mdcell.type_atom_counts();
        for (std::size_t it = 0; it < type_atom_counts.size(); ++it)
        {
            for (int ia = 0; ia < type_atom_counts[it]; ++ia)
            {
                ModuleBase::Vector3<double> random_value;
                if (mdcell.mpi_rank() == 0)
                {
                    for (int k = 0; k < 3; ++k)
                    {
                        random_value[k] = static_cast<double>(std::rand()) / RAND_MAX - 0.5;
                    }
                }
                MPI_Bcast(&random_value.x, 3, MPI_DOUBLE, 0, mdcell.communicator());
                for (int i = 0; i < mdcell.nlocal(); ++i)
                {
                    const LocalAtom& atom = mdcell.owned_atoms()[static_cast<std::size_t>(i)];
                    if (atom.type == it && atom.type_index == ia)
                    {
                        random_values[static_cast<std::size_t>(i)] = random_value;
                        break;
                    }
                }
            }
        }
    }
    else
#endif
    for (int i = 0; i < mdcell.nlocal(); ++i)
    {
        for (int k = 0; k < 3; ++k)
        {
            random_values[static_cast<std::size_t>(i)][k] = static_cast<double>(std::rand()) / RAND_MAX - 0.5;
        }
    }

    for (int i = 0; i < mdcell.nlocal(); ++i)
    {
        const LocalAtom& atom = mdcell.owned_atoms()[static_cast<std::size_t>(i)];
        ModuleBase::Vector3<double> fictitious_force = -atom.mass * atom.vel / md_damp;
        for (int j = 0; j < 3; ++j)
        {
            fictitious_force[j] += sqrt(24.0 * t_target * atom.mass / md_damp / md_dt)
                                   * random_values[static_cast<std::size_t>(i)][j];
        }
        total_force[i] = atom.force + fictitious_force;
    }
}
