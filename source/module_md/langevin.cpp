#include "langevin.h"

#include "md_func.h"
#include "module_base/parallel_common.h"
#include "module_base/timer.h"

Langevin::Langevin(MD_parameters& MD_para_in, UnitCell& unit_in) : MD_base(MD_para_in, unit_in)
{
    // convert to a.u. unit
    mdp.md_damp /= ModuleBase::AU_to_FS;

    total_force = new ModuleBase::Vector3<double>[ucell.nat];
}

Langevin::~Langevin()
{
    delete[] total_force;
}

void Langevin::setup(ModuleESolver::ESolver* p_esolver, const int& my_rank, const std::string& global_readin_dir)
{
    ModuleBase::TITLE("Langevin", "setup");
    ModuleBase::timer::tick("Langevin", "setup");

    MD_base::setup(p_esolver, my_rank, global_readin_dir);

    post_force();

    ModuleBase::timer::tick("Langevin", "setup");
}

void Langevin::first_half(const int& my_rank, std::ofstream& ofs)
{
    ModuleBase::TITLE("Langevin", "first_half");
    ModuleBase::timer::tick("Langevin", "first_half");

    MD_base::update_vel(total_force, my_rank);
    MD_base::update_pos(my_rank);

    ModuleBase::timer::tick("Langevin", "first_half");
}

void Langevin::second_half(const int& my_rank)
{
    ModuleBase::TITLE("Langevin", "second_half");
    ModuleBase::timer::tick("Langevin", "second_half");

    post_force();
    MD_base::update_vel(total_force, my_rank);

    ModuleBase::timer::tick("Langevin", "second_half");
}

void Langevin::outputMD(std::ofstream& ofs, const bool& cal_stress, const int& my_rank)
{
    MD_base::outputMD(ofs, cal_stress, my_rank);
}

void Langevin::write_restart(const int& my_rank, const std::string& global_out_dir)
{
    MD_base::write_restart(my_rank, global_out_dir);
}

void Langevin::restart(const int& my_rank, const std::string& global_readin_dir)
{
    MD_base::restart(my_rank, global_readin_dir);
}

void Langevin::post_force()
{
    if (GlobalV::MY_RANK == 0)
    {
        double t_target = MD_func::target_temp(step_ + step_rst_, mdp.md_nstep, mdp.md_tfirst, mdp.md_tlast);
        ModuleBase::Vector3<double> fictitious_force;
        for (int i = 0; i < ucell.nat; ++i)
        {
            fictitious_force = -allmass[i] * vel[i] / mdp.md_damp;
            for (int j = 0; j < 3; ++j)
            {
                fictitious_force[j] += sqrt(24.0 * t_target * allmass[i] / mdp.md_damp / mdp.md_dt)
                                       * (static_cast<double>(std::rand()) / RAND_MAX - 0.5);
            }
            total_force[i] = force[i] + fictitious_force;
        }
    }

#ifdef __MPI
    MPI_Bcast(total_force, ucell.nat * 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
}