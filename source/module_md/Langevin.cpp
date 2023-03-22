#include "Langevin.h"
#include "MD_func.h"
#include "module_base/parallel_common.h"
#include "module_base/timer.h"

Langevin::Langevin(MD_parameters& MD_para_in, UnitCell &unit_in) : MDrun(MD_para_in, unit_in)
{
    // convert to a.u. unit
    mdp.md_damp /= ModuleBase::AU_to_FS;

    total_force = new ModuleBase::Vector3<double> [ucell.nat];
}

Langevin::~Langevin()
{
    delete []total_force;
}

void Langevin::setup(ModuleESolver::ESolver *p_ensolve)
{
    ModuleBase::TITLE("Langevin", "setup");
    ModuleBase::timer::tick("Langevin", "setup");
    
    MDrun::setup(p_ensolve);

    post_force();

    ModuleBase::timer::tick("Langevin", "setup");
}

void Langevin::first_half()
{
    ModuleBase::TITLE("Langevin", "first_half");
    ModuleBase::timer::tick("Langevin", "first_half");

    MDrun::update_vel(total_force);
    MDrun::update_pos();

    ModuleBase::timer::tick("Langevin", "first_half");
}

void Langevin::second_half()
{
    ModuleBase::TITLE("Langevin", "second_half");
    ModuleBase::timer::tick("Langevin", "second_half");

    post_force();
    MDrun::update_vel(total_force);

    ModuleBase::timer::tick("Langevin", "second_half");
}

void Langevin::outputMD(std::ofstream &ofs, bool cal_stress)
{
    MDrun::outputMD(ofs, cal_stress);
}

void Langevin::write_restart()
{
    MDrun::write_restart();
}

void Langevin::restart()
{
    MDrun::restart();
}

void Langevin::post_force()
{
    if(GlobalV::MY_RANK==0)
    {
        double t_target = MD_func::target_temp(step_ + step_rst_, mdp.md_tfirst, mdp.md_tlast);
        ModuleBase::Vector3<double> fictitious_force;
        for(int i=0; i<ucell.nat; ++i)
        {
            fictitious_force = - allmass[i] * vel[i] / mdp.md_damp;
            for(int j=0; j<3; ++j)
            {
                fictitious_force[j] += sqrt(24.0 * t_target * allmass[i] / mdp.md_damp / mdp.md_dt) * (rand()/double(RAND_MAX) - 0.5);
            }
            total_force[i] = force[i] + fictitious_force;
        }
    }

#ifdef __MPI
    MPI_Bcast(total_force, ucell.nat*3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
}