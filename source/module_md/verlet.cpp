#include "verlet.h"
#include "MD_func.h"
#include "../module_base/timer.h"

Verlet::Verlet(MD_parameters& MD_para_in, UnitCell &unit_in) : MDrun(MD_para_in, unit_in){}

Verlet::~Verlet(){}

void Verlet::setup(ModuleESolver::ESolver *p_ensolve)
{
    ModuleBase::TITLE("Verlet", "setup");
    ModuleBase::timer::tick("Verlet", "setup");
    
    MDrun::setup(p_ensolve);

    ModuleBase::timer::tick("Verlet", "setup");
}

void Verlet::first_half()
{
    ModuleBase::TITLE("Verlet", "first_half");
    ModuleBase::timer::tick("Verlet", "first_half");

    MDrun::first_half();

    ModuleBase::timer::tick("Verlet", "first_half");
}

void Verlet::second_half()
{
    ModuleBase::TITLE("Verlet", "second_half");
    ModuleBase::timer::tick("Verlet", "second_half");

    MDrun::second_half();
    apply_thermostat();

    ModuleBase::timer::tick("Verlet", "second_half");
}

void Verlet::apply_thermostat()
{
    double t_target = 0;
    t_current = MD_func::current_temp(kinetic, ucell.nat, frozen_freedom_, allmass, vel);

    if(mdp.md_thermostat == "NVE"){}
    else if(mdp.md_thermostat == "Rescaling")
    {
        t_target = MD_func::target_temp(step_ + step_rst_, mdp.md_tfirst, mdp.md_tlast) * ModuleBase::Hartree_to_K;
        if(abs(t_target - t_current) > mdp.md_tolerance)
        {
            thermalize(0, t_current, t_target);
        }
    }
    else if(mdp.md_thermostat == "Rescale_v")
    {
        if((step_+step_rst_) % mdp.md_nraise == 0)
        {
            t_target = MD_func::target_temp(step_ + step_rst_, mdp.md_tfirst, mdp.md_tlast) * ModuleBase::Hartree_to_K;
            thermalize(0, t_current, t_target);
        }
    }
    else if(mdp.md_thermostat == "Anderson")
    {
        if(GlobalV::MY_RANK==0)
        {
            double deviation;
            for(int i=0; i<ucell.nat; ++i)
            {  
                if(rand()/double(RAND_MAX) <= 1.0 / mdp.md_nraise)
                {
                    deviation = sqrt(mdp.md_tlast / allmass[i]);
                    for(int k=0; k<3; ++k)
                    {
                        if(ionmbl[i][k]) 
                        {
                            vel[i][k] = deviation * MD_func::gaussrand();
                        }
                    }
                }
            }
        }
#ifdef __MPI
        MPI_Bcast(vel, ucell.nat*3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
    }
    else if(mdp.md_thermostat == "Berendsen")
    {
        t_target = MD_func::target_temp(step_ + step_rst_, mdp.md_tfirst, mdp.md_tlast) * ModuleBase::Hartree_to_K;
        thermalize(mdp.md_nraise, t_current, t_target);
    }
    else
    {
        ModuleBase::WARNING_QUIT("Verlet", "No such thermostat !");
    }
}

void Verlet::thermalize(const int &nraise, const double &current_temp, const double &target_temp)
{
    double fac = 0;
    if(nraise > 0 && current_temp > 0 && target_temp >0)
    {
        fac = sqrt(1 + (target_temp / current_temp - 1) / nraise);
    }
    else if(nraise == 0 && current_temp > 0 && target_temp >0)
    {
        fac = sqrt(target_temp / current_temp);
    }

    for(int i=0; i < ucell.nat; ++i)
    {
        vel[i] *= fac;
    }
}

void Verlet::outputMD(std::ofstream &ofs, bool cal_stress)
{
    MDrun::outputMD(ofs, cal_stress);
}

void Verlet::write_restart()
{
    MDrun::write_restart();
}

void Verlet::restart()
{
    MDrun::restart();
}