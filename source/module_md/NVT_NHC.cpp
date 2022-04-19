#include "NVT_NHC.h"
#include "MD_func.h"
#ifdef __MPI
#include "mpi.h"
#endif
#include "../module_base/timer.h"

NVT_NHC::NVT_NHC(MD_parameters& MD_para_in, UnitCell_pseudo &unit_in) : Verlet(MD_para_in, unit_in)
{
    if(mdp.md_tfirst == 0)
    {
        std::cout << " md_tfirst must be larger than 0 in NHC !!! " << std::endl;
        ModuleBase::WARNING_QUIT("NVT_NHC", " md_tfirst must be larger than 0 in NHC !!! ");
    }

    // init NHC
    Q = new double [mdp.md_mnhc];
	G = new double [mdp.md_mnhc];
	eta = new double [mdp.md_mnhc];
	veta = new double [mdp.md_mnhc];

    for(int i=0; i<mdp.md_mnhc; ++i)
    {
        eta[i] = veta[i] = G[i] = 0;
    }

    //w[0] = 1;

    w[0] = 0.784513610477560;
	w[6] = 0.784513610477560;
	w[1] = 0.235573213359357;
	w[5] = 0.235573213359357;
	w[2] = -1.17767998417887;
	w[4] = -1.17767998417887;
	w[3] = 1-w[0]-w[1]-w[2]-w[4]-w[5]-w[6];
}

NVT_NHC::~NVT_NHC()
{
    delete []Q;
    delete []G;
    delete []eta;
    delete []veta;
}

void NVT_NHC::setup(ModuleESolver::ESolver *p_ensolve)
{
    ModuleBase::TITLE("NVT_NHC", "setup");
    ModuleBase::timer::tick("NVT_NHC", "setup");

    Verlet::setup(p_ensolve);

    temp_target();
    
    update_mass();
    
    for(int m=1; m<mdp.md_mnhc; ++m)
    {
        G[m] = (Q[m-1]*veta[m-1]*veta[m-1]-t_target) / Q[m];
    }

    ModuleBase::timer::tick("NVT_NHC", "setup");
}

void NVT_NHC::first_half()
{
    ModuleBase::TITLE("NVT_NHC", "first_half");
    ModuleBase::timer::tick("NVT_NHC", "first_half");

    // update target T
    temp_target();

    integrate();

    Verlet::first_half();

    ModuleBase::timer::tick("NVT_NHC", "first_half");
}

void NVT_NHC::second_half()
{
    ModuleBase::TITLE("NVT_NHC", "second_half");
    ModuleBase::timer::tick("NVT_NHC", "second_half");

    Verlet::second_half();

    integrate();

    ModuleBase::timer::tick("NVT_NHC", "second_half");
}

void NVT_NHC::outputMD(std::ofstream &ofs)
{
    Verlet::outputMD(ofs);
}

void NVT_NHC::write_restart()
{
    if(!GlobalV::MY_RANK)
    {
		std::stringstream ssc;
		ssc << GlobalV::global_out_dir << "Restart_md.dat";
		std::ofstream file(ssc.str().c_str());

        file << step_ + step_rst_ << std::endl;
        file << mdp.md_mnhc << std::endl;
        for(int i=0; i<mdp.md_mnhc; ++i)
        {
            file << eta[i] << "   ";
        }
        file << std::endl;
        for(int i=0; i<mdp.md_mnhc; ++i)
        {
            file << veta[i] << "   ";
        }
		file.close();
	}
#ifdef __MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif
}

void NVT_NHC::restart()
{
    if(!GlobalV::MY_RANK)
    {
		std::stringstream ssc;
		ssc << GlobalV::global_out_dir << "Restart_md.dat";
		std::ifstream file(ssc.str().c_str());

        if(!file)
		{
			std::cout<< "Please ensure whether 'Restart_md.dat' exists !" << std::endl;
            ModuleBase::WARNING_QUIT("NVT_NHC", "no Restart_md.dat ！");
		}
        double Mnum;
		file >> step_rst_ >> Mnum;
        if(Mnum != mdp.md_mnhc)
		{
			std::cout<< "Num of NHC is not the same !" << std::endl;
            ModuleBase::WARNING_QUIT("NVT_NHC", "no Restart_md.dat ！");
		}
        for(int i=0; i<mdp.md_mnhc; ++i)
        {
            file >> eta[i];
        }
        for(int i=0; i<mdp.md_mnhc; ++i)
        {
            file >> veta[i];
        }

		file.close();
	}

#ifdef __MPI
	MPI_Bcast(&step_rst_, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(eta, mdp.md_mnhc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(veta, mdp.md_mnhc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
}

void NVT_NHC::integrate()
{
    double scale = 1.0;
    kinetic = MD_func::GetAtomKE(ucell.nat, vel, allmass);
    double KE = kinetic;

    // update mass
    update_mass();

    // update force
    if(Q[0] > 0) 
    {
        G[0] = (2*KE - (3*ucell.nat - frozen_freedom_)*t_target) / Q[0];
    }
    else 
    {
        G[0] = 0;
    }

    for(int i=0; i<nc; ++i)
    {
        for(int j=0; j<nsy; ++j)
        {
            double delta = w[j] * mdp.md_dt / nc;

            // propogate veta
            veta[mdp.md_mnhc-1] += G[mdp.md_mnhc-1] * delta /4.0;

            for(int m=mdp.md_mnhc-2; m>=0; --m)
            {
                double aa = exp(-veta[m]*delta/8.0);
                veta[m] = veta[m] * aa * aa + G[m] * aa * delta /4.0;
            }

            scale *= exp(-veta[0]*delta/2.0);
            if(!isfinite(scale))
            {
                std::cout << " Please set a proper md_tfreq !!! " << std::endl;
                ModuleBase::WARNING_QUIT("NVT_NHC", " Please set a proper md_tfreq !!! ");
            }
            
            KE = kinetic * scale * scale;

            // update force
            if(Q[0] > 0) 
            {
                G[0] = (2*KE - (3*ucell.nat - frozen_freedom_)*t_target) / Q[0];
            }
            else 
            {
                G[0] = 0;
            }

            // propogate eta
            for(int m=0; m<mdp.md_mnhc; ++m)
            {
                eta[m] += veta[m] * delta / 2.0;
            }

            // propogate veta
            for(int m=0; m<mdp.md_mnhc-1; ++m)
            {
                double aa = exp(-veta[m+1]*delta/8.0);
                veta[m] = veta[m] * aa * aa + G[m] * aa * delta /4.0;

                G[m+1] = (Q[m]*veta[m]*veta[m]-t_target) / Q[m+1];
            }
            veta[mdp.md_mnhc-1] += G[mdp.md_mnhc-1] * delta /4.0;
        }
    }
    
    for(int i=0; i<ucell.nat; ++i)
    {
        vel[i] *= scale;
    }
}

void NVT_NHC::temp_target()
{
    double delta = (double)(step_ + step_rst_) / GlobalV::MD_NSTEP;
    t_target = mdp.md_tfirst + delta * (mdp.md_tlast - mdp.md_tfirst);
}

void NVT_NHC::update_mass()
{
    Q[0] = (3*ucell.nat - frozen_freedom_) * t_target / mdp.md_tfreq / mdp.md_tfreq;
    for(int m=1; m<mdp.md_mnhc; ++m)
    {
        Q[m] = t_target / mdp.md_tfreq / mdp.md_tfreq;
    }
}