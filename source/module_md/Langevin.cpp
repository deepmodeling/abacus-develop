#include "Langevin.h"
#include "MD_func.h"
#include "../src_parallel/parallel_common.h"
#include "../module_base/timer.h"

Langevin::Langevin(MD_parameters& MD_para_in, UnitCell &unit_in) : MDrun(MD_para_in, unit_in)
{
    // convert to a.u. unit
    mdp.md_damp /= ModuleBase::AU_to_FS;

    fictitious_force = new ModuleBase::Vector3<double> [ucell.nat];
}

Langevin::~Langevin()
{
    delete []fictitious_force;
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

    if(GlobalV::MY_RANK==0)
    {
        for(int i=0; i<ucell.nat; ++i)
        {
            for(int k=0; k<3; ++k)
            {
                if(ionmbl[i][k])
                {
                    vel[i][k] += 0.5 * (force[i][k] + fictitious_force[i][k]) * mdp.md_dt / allmass[i];
                    pos[i][k] += vel[i][k] * mdp.md_dt;
                }
            }
        }
    }

#ifdef __MPI
    MPI_Bcast(pos , ucell.nat*3,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(vel , ucell.nat*3,MPI_DOUBLE,0,MPI_COMM_WORLD);
#endif

    ucell.update_pos_tau(pos);
    ucell.periodic_boundary_adjustment();
    MD_func::InitPos(ucell, pos);

    ModuleBase::timer::tick("Langevin", "first_half");
}

void Langevin::second_half()
{
    ModuleBase::TITLE("Langevin", "second_half");
    ModuleBase::timer::tick("Langevin", "second_half");

    post_force();

    for(int i=0; i<ucell.nat; ++i)
    {
        for(int k=0; k<3; ++k)
        {
            if(ionmbl[i][k])
            {
                vel[i][k] += 0.5 * (force[i][k] + fictitious_force[i][k]) * mdp.md_dt / allmass[i];
            }
        }
    }

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
    double t_target = MD_func::target_temp(step_ + step_rst_, mdp.md_tfirst, mdp.md_tlast);

    if(GlobalV::MY_RANK==0)
    {
        for(int i=0; i<ucell.nat; ++i)
        {
            fictitious_force[i] = - allmass[i] * vel[i] / mdp.md_damp;
            for(int j=0; j<3; ++j)
            {
                fictitious_force[i][j] += sqrt(24.0 * t_target * allmass[i] / mdp.md_damp / mdp.md_dt) * (rand()/double(RAND_MAX) - 0.5);
            }
        }
    }

#ifdef __MPI
    MPI_Bcast(fictitious_force, ucell.nat*3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
}