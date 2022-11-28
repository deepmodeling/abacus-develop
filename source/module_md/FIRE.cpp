#include "FIRE.h"
#include "MD_func.h"
#ifdef __MPI
#include "mpi.h"
#endif
#include "../module_base/timer.h"

FIRE::FIRE(MD_parameters& MD_para_in, UnitCell &unit_in) : MDrun(MD_para_in, unit_in)
{
    dt_max = -1.0;
    alpha_start = 0.10;
    alpha = alpha_start;

    finc = 1.1;
    fdec = 0.5;
    f_alpha = 0.99;
    N_min = 4;
    negative_count = 0;
}

FIRE::~FIRE(){}

void FIRE::setup(ModuleESolver::ESolver *p_ensolve)
{
    ModuleBase::TITLE("FIRE", "setup");
    ModuleBase::timer::tick("FIRE", "setup");
    
    MDrun::setup(p_ensolve);

    check_force();

    ModuleBase::timer::tick("FIRE", "setup");
}

void FIRE::first_half()
{
    ModuleBase::TITLE("FIRE", "first_half");
    ModuleBase::timer::tick("FIRE", "first_half");
    if(GlobalV::MY_RANK == 0)
    {
    for(int i=0; i<ucell.nat; ++i)
    {
        for(int k=0; k<3; ++k)
        {
            if(ionmbl[i][k])
            {
                vel[i][k] += 0.5*force[i][k]*mdp.md_dt/allmass[i];
            }
        }
    }

    check_FIRE();

    for(int i=0; i<ucell.nat; ++i)
    {
        for(int k=0; k<3; ++k)
        {
            if(ionmbl[i][k])
            {
                pos[i][k] += vel[i][k]*mdp.md_dt;
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

    ModuleBase::timer::tick("FIRE", "first_half");
}

void FIRE::second_half()
{
    ModuleBase::TITLE("FIRE", "second_half");
    ModuleBase::timer::tick("FIRE", "second_half");

    MDrun::second_half();

    check_force();

    ModuleBase::timer::tick("FIRE", "second_half");
}

void FIRE::outputMD(std::ofstream &ofs, bool cal_stress)
{
    MDrun::outputMD(ofs, cal_stress);

    ofs << " LARGEST GRAD (eV/A)  : " 
        << max * ModuleBase::Hartree_to_eV * ModuleBase::ANGSTROM_AU << std::endl;
    std::cout << " LARGEST GRAD (eV/A)  : " 
        << max * ModuleBase::Hartree_to_eV * ModuleBase::ANGSTROM_AU << std::endl;
}

void FIRE::write_restart()
{
    if(!GlobalV::MY_RANK)
    {
		std::stringstream ssc;
		ssc << GlobalV::global_out_dir << "Restart_md.dat";
		std::ofstream file(ssc.str().c_str());

        file << step_ + step_rst_ << std::endl;
        file << alpha << std::endl;
        file << negative_count << std::endl;
        file << dt_max << std::endl;
        file << mdp.md_dt << std::endl;
		file.close();
	}
#ifdef __MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif
}

void FIRE::restart()
{
    bool ok = true;

    if(!GlobalV::MY_RANK)
    {
        std::stringstream ssc;
        ssc << GlobalV::global_readin_dir << "Restart_md.dat";
        std::ifstream file(ssc.str().c_str());

        if(!file)
        {
            ok = false;
        }

        if(ok)
        {
            file >> step_rst_ >> alpha >> negative_count >> dt_max >> mdp.md_dt;
            file.close();
        }
    }

#ifdef __MPI
    MPI_Bcast(&ok, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

    if(!ok)
    {
        ModuleBase::WARNING_QUIT("mdrun", "no Restart_md.dat !");
    }

#ifdef __MPI
    MPI_Bcast(&step_rst_, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&alpha, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&negative_count, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&dt_max, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&mdp.md_dt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
}

void FIRE::check_force()
{
    max = 0;

    for(int i=0; i<ucell.nat; ++i)
    {
        for(int j=0; j<3; ++j)
        {
            if( max < abs(force[i][j]) )
            {
                max = abs(force[i][j]);
            }
        }
    }

    if(2.0*max < GlobalV::FORCE_THR)
    {
        stop = true;
    }
}

void FIRE::check_FIRE()
{
	double P = 0;
	double sumforce = 0;
	double normvel = 0;

    if(dt_max < 0) dt_max = 2.5 * mdp.md_dt;   //initial dt_max
	
	for(int i =0; i<ucell.nat; ++i)
    {
		P += vel[i].x*force[i].x+vel[i].y*force[i].y+vel[i].z*force[i].z;
		sumforce += force[i].norm2();
		normvel += vel[i].norm2();
	}
	
	sumforce = sqrt(sumforce);
	normvel = sqrt(normvel);
	
	for(int i=0; i<ucell.nat; ++i)
    {
        for(int j=0; j<3; ++j)
        {
            vel[i][j] = (1.0-alpha) * vel[i][j] + alpha * force[i][j] / sumforce * normvel;
        }
	}
	
	if(P > 0)
	{
	    negative_count++;
		if(negative_count >= N_min)
        {
			mdp.md_dt = min(mdp.md_dt*finc, dt_max);
			alpha *= f_alpha;
		}
	}
	else
	{
		mdp.md_dt *= fdec;
		negative_count = 0;
		
		for(int i=0; i<ucell.nat; ++i)
        {
			for(int j=0; j<3; ++j)
            {
                vel[i][j] = 0;
            }
		}
		
		alpha = alpha_start;
	}
}

