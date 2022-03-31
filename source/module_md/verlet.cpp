#include "verlet.h"
#include "MD_func.h"
#ifdef __MPI
#include "mpi.h"
#endif
#include "../module_base/timer.h"
#include "module_esolver/esolver.h"

Verlet::Verlet(MD_parameters& MD_para_in, UnitCell_pseudo &unit_in):
    mdp(MD_para_in),
    ucell(unit_in)
{
    stop = false;

    allmass = new double [ucell.nat];
    pos = new ModuleBase::Vector3<double> [ucell.nat];
    vel = new ModuleBase::Vector3<double> [ucell.nat];
    ionmbl = new ModuleBase::Vector3<int> [ucell.nat];
    force = new ModuleBase::Vector3<double> [ucell.nat];
    virial.create(3,3);
    stress.create(3,3);

    // convert to a.u. unit
    mdp.md_dt /= ModuleBase::AU_to_FS;
    mdp.md_tfirst /= ModuleBase::Hartree_to_K;
    mdp.md_tlast /= ModuleBase::Hartree_to_K;

    // LJ parameters
    mdp.lj_rcut *= ModuleBase::ANGSTROM_AU;
	mdp.lj_epsilon /= ModuleBase::Hartree_to_eV;
	mdp.lj_sigma *= ModuleBase::ANGSTROM_AU;

    step_ = 0;
    step_rst_ = 0;
    if(mdp.md_restart) unit_in.init_vel = 1;

    MD_func::InitPos(ucell, pos);
    MD_func::InitVel(ucell, mdp.md_tfirst, allmass, frozen_freedom_, ionmbl, vel);
}

Verlet::~Verlet()
{
    delete []allmass;
    delete []pos;
    delete []vel;
    delete []ionmbl;
    delete []force;
}

void Verlet::setup(ModuleESolver::ESolver *p_esolver)
{
    if(mdp.md_restart)
    {
        restart();
    }

    MD_func::force_virial(p_esolver, step_, mdp, ucell, potential, force, virial);
    MD_func::kinetic_stress(ucell, vel, allmass, kinetic, stress);
    stress += virial;

    temperature_ = 2*kinetic/(double(3*ucell.nat-frozen_freedom_))*ModuleBase::Hartree_to_K;
}

void Verlet::first_half()
{
    for(int i=0; i<ucell.nat; ++i)
    {
        for(int k=0; k<3; ++k)
        {
            if(ionmbl[i][k])
            {
                vel[i][k] += 0.5*force[i][k]*mdp.md_dt/allmass[i];
                pos[i][k] += vel[i][k]*mdp.md_dt;
            }
        }
    }

    ucell.update_pos_tau(pos);
    ucell.periodic_boundary_adjustment();
    MD_func::InitPos(ucell, pos);
}

void Verlet::second_half()
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
}

void Verlet::outputMD()
{
    if(GlobalV::MY_RANK) return;

	std::ofstream ofs;
	if(step_ + step_rst_ == 0)
	{
		ofs.open("run.log", ios::trunc);
	}
	else
	{
		ofs.open("run.log", ios::app);
	}

    temperature_ = 2*kinetic/(double(3*ucell.nat-frozen_freedom_))*ModuleBase::Hartree_to_K;

    const double unit_transform = ModuleBase::HARTREE_SI / pow(ModuleBase::BOHR_RADIUS_SI,3) * 1.0e-8;
    double press = 0.0;
    for(int i=0;i<3;i++)
    {
        press += stress(i,i)/3;
    }

    // std::cout << std::setprecision(6) << std::setiosflags(ios::showpos) << std::setiosflags(ios::fixed);
    std::cout << " " << std::left << std::setw(20) << "Energy" 
            << std::left << std::setw(20) << "Potential" 
            << std::left << std::setw(20) << "Kinetic" 
            << std::left << std::setw(20) << "Temperature" 
            << std::left << std::setw(20) << "Pressure (KBAR)" <<std::endl;
    std::cout << " " << std::left << std::setw(20) << potential+kinetic
            << std::left << std::setw(20) << potential
            << std::left << std::setw(20) << kinetic
            << std::left << std::setw(20) << temperature_
            << std::left << std::setw(20) << press*unit_transform <<std::endl;
	std::cout << " ------------------------------------------------------------" << std::endl;
    MD_func::outStress(virial, stress);

    if(step_ == 0)
	{
		ofs << std::left << std::setw(20) << "Step" 
            << std::left << std::setw(20) << "Energy" 
            << std::left << std::setw(20) << "Potential" 
            << std::left << std::setw(20) << "Kinetic" 
            << std::left << std::setw(20) << "Temperature" 
            << std::left << std::setw(20) << "Pressure (KBAR)" <<std::endl;
	}
    ofs << std::left << std::setw(20) << step_ + step_rst_
        << std::left << std::setw(20) << potential+kinetic
        << std::left << std::setw(20) << potential
        << std::left << std::setw(20) << kinetic
        << std::left << std::setw(20) << temperature_
        << std::left << std::setw(20) << press*unit_transform <<std::endl;
    ofs.close();
}

void Verlet::write_restart()
{
    if(!GlobalV::MY_RANK)
    {
		std::stringstream ssc;
		ssc << GlobalV::global_out_dir << "Restart_md.dat";
		std::ofstream file(ssc.str().c_str());

        file << step_ + step_rst_ << std::endl;
		file.close();
	}
#ifdef __MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif
}

void Verlet::restart()
{
    if(!GlobalV::MY_RANK)
    {
		std::stringstream ssc;
		ssc << GlobalV::global_out_dir << "Restart_md.dat";
		std::ifstream file(ssc.str().c_str());

        if(!file)
		{
			std::cout<< "please ensure whether 'Restart_md.dat' exists!" << std::endl;
            ModuleBase::WARNING_QUIT("verlet", "no Restart_md.dat ！");
		}

		file >> step_rst_;

		file.close();
	}

#ifdef __MPI
	MPI_Bcast(&step_rst_, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
}