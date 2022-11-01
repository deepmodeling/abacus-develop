#include "Nose_Hoover.h"
#include "MD_func.h"
#ifdef __MPI
#include "mpi.h"
#endif
#include "../module_base/timer.h"

Nose_Hoover::Nose_Hoover(MD_parameters& MD_para_in, UnitCell_pseudo &unit_in) : MDrun(MD_para_in, unit_in)
{
    if(mdp.md_tfirst == 0)
    {
        ModuleBase::WARNING_QUIT("Nose_Hoover", " md_tfirst must be larger than 0 in NHC !!! ");
    }

    // init NHC
    Q = new double [mdp.md_tchain];
	G = new double [mdp.md_tchain];
	eta = new double [mdp.md_tchain];
	veta = new double [mdp.md_tchain];

    for(int i=0; i<mdp.md_tchain; ++i)
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

Nose_Hoover::~Nose_Hoover()
{
    delete []Q;
    delete []G;
    delete []eta;
    delete []veta;
}

void Nose_Hoover::setup(ModuleESolver::ESolver *p_ensolve)
{
    ModuleBase::TITLE("Nose_Hoover", "setup");
    ModuleBase::timer::tick("Nose_Hoover", "setup");

    MDrun::setup(p_ensolve);

    t_target = MD_func::target_temp(step_ + step_rst_, mdp.md_tfirst, mdp.md_tlast);
    
    update_mass();
    
    for(int m=1; m<mdp.md_tchain; ++m)
    {
        G[m] = (Q[m-1]*veta[m-1]*veta[m-1]-t_target) / Q[m];
    }

    ModuleBase::timer::tick("Nose_Hoover", "setup");
}

void Nose_Hoover::first_half()
{
    ModuleBase::TITLE("Nose_Hoover", "first_half");
    ModuleBase::timer::tick("Nose_Hoover", "first_half");

    // update target T
    t_target = MD_func::target_temp(step_ + step_rst_, mdp.md_tfirst, mdp.md_tlast);

    integrate();

    MDrun::first_half();

    ModuleBase::timer::tick("Nose_Hoover", "first_half");
}

void Nose_Hoover::second_half()
{
    ModuleBase::TITLE("Nose_Hoover", "second_half");
    ModuleBase::timer::tick("Nose_Hoover", "second_half");

    MDrun::second_half();

    integrate();

    ModuleBase::timer::tick("Nose_Hoover", "second_half");
}

void Nose_Hoover::outputMD(std::ofstream &ofs, bool cal_stress)
{
    MDrun::outputMD(ofs, cal_stress);
}

void Nose_Hoover::write_restart()
{
    if(!GlobalV::MY_RANK)
    {
		std::stringstream ssc;
		ssc << GlobalV::global_out_dir << "Restart_md.dat";
		std::ofstream file(ssc.str().c_str());

        file << step_ + step_rst_ << std::endl;
        file << mdp.md_tchain << std::endl;
        for(int i=0; i<mdp.md_tchain; ++i)
        {
            file << eta[i] << "   ";
        }
        file << std::endl;
        for(int i=0; i<mdp.md_tchain; ++i)
        {
            file << veta[i] << "   ";
        }
		file.close();
	}
#ifdef __MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif
}

void Nose_Hoover::restart()
{
    bool ok = true;
    bool ok2 = true;

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
            double Mnum;
            file >> step_rst_ >> Mnum;

            if( Mnum != mdp.md_tchain )
            {
                ok2 = false;
            }

            if(ok2)
            {
                for(int i=0; i<mdp.md_tchain; ++i)
                {
                    file >> eta[i];
                }
                for(int i=0; i<mdp.md_tchain; ++i)
                {
                    file >> veta[i];
                }
            }

            file.close();
        }
    }

#ifdef __MPI
    MPI_Bcast(&ok, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&ok2, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

    if(!ok)
    {
        ModuleBase::WARNING_QUIT("mdrun", "no Restart_md.dat !");
    }
    if(!ok2)
    {
        ModuleBase::WARNING_QUIT("mdrun", "Num of NHC is not the same !");
    }

#ifdef __MPI
	MPI_Bcast(&step_rst_, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(eta, mdp.md_tchain, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(veta, mdp.md_tchain, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
}

void Nose_Hoover::integrate()
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
            veta[mdp.md_tchain-1] += G[mdp.md_tchain-1] * delta /4.0;

            for(int m=mdp.md_tchain-2; m>=0; --m)
            {
                double aa = exp(-veta[m]*delta/8.0);
                veta[m] = veta[m] * aa * aa + G[m] * aa * delta /4.0;
            }

            scale *= exp(-veta[0]*delta/2.0);
            if(!isfinite(scale))
            {
                ModuleBase::WARNING_QUIT("Nose_Hoover", " Please set a proper md_tfreq !!! ");
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
            for(int m=0; m<mdp.md_tchain; ++m)
            {
                eta[m] += veta[m] * delta / 2.0;
            }

            // propogate veta
            for(int m=0; m<mdp.md_tchain-1; ++m)
            {
                double aa = exp(-veta[m+1]*delta/8.0);
                veta[m] = veta[m] * aa * aa + G[m] * aa * delta /4.0;

                G[m+1] = (Q[m]*veta[m]*veta[m]-t_target) / Q[m+1];
            }
            veta[mdp.md_tchain-1] += G[mdp.md_tchain-1] * delta /4.0;
        }
    }
    
    for(int i=0; i<ucell.nat; ++i)
    {
        vel[i] *= scale;
    }
}

void Nose_Hoover::update_mass()
{
    Q[0] = (3*ucell.nat - frozen_freedom_) * t_target / mdp.md_tfreq / mdp.md_tfreq;
    for(int m=1; m<mdp.md_tchain; ++m)
    {
        Q[m] = t_target / mdp.md_tfreq / mdp.md_tfreq;
    }
}