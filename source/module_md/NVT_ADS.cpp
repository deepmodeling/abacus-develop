#include "NVT_ADS.h"
#include "MD_func.h"
#ifdef __MPI
#include "mpi.h"
#endif
#include "../module_base/timer.h"

NVT_ADS::NVT_ADS(MD_parameters& MD_para_in, UnitCell_pseudo &unit_in) : Verlet(MD_para_in, unit_in)
{
    nraise = mdp.md_tfreq * mdp.md_dt;
}

NVT_ADS::~NVT_ADS(){}

void NVT_ADS::setup(ModuleESolver::ESolver *p_ensolve)
{
    ModuleBase::TITLE("NVT_ADS", "setup");
    ModuleBase::timer::tick("NVT_ADS", "setup");

    Verlet::setup(p_ensolve);

    ModuleBase::timer::tick("NVT_ADS", "setup");
}

void NVT_ADS::first_half()
{
    ModuleBase::TITLE("NVT_ADS", "first_half");
    ModuleBase::timer::tick("NVT_ADS", "first_half");

    Verlet::first_half();

    ModuleBase::timer::tick("NVT_ADS", "first_half");
}

void NVT_ADS::second_half()
{
    ModuleBase::TITLE("NVT_ADS", "second_half");
    ModuleBase::timer::tick("NVT_ADS", "second_half");

    Verlet::second_half();

    if(GlobalV::MY_RANK==0)
    {
        double deviation;
        for(int i=0; i<ucell.nat; ++i)
        {
            if(rand()/double(RAND_MAX) <= nraise)
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

    ModuleBase::timer::tick("NVT_ADS", "second_half");
}

void NVT_ADS::outputMD(std::ofstream &ofs)
{
    Verlet::outputMD(ofs);
}

void NVT_ADS::write_restart()
{
    Verlet::write_restart();
}

void NVT_ADS::restart()
{
    Verlet::restart();
}