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

    // init NPT related variables
    for(int i=0; i<6; ++i)
    {
        pstart[i] = pstop[i] = pfreq[i] = p_target[i] = pflag[i] = 0;
    }

    // determine the NPT methods
    if(mdp.md_pmode == "iso")
    {
        mdp.md_pcouple = "xyz";
        pstart[0] = pstart[1] = pstart[2] = mdp.md_pfirst;
        pstop[0] = pstop[1] = pstop[2] = mdp.md_plast;
        pfreq[0] = pfreq[1] = pfreq[2] = mdp.md_pfreq;
        pflag[0] = pflag[1] = pflag[2] = 1;
    }
    else if(mdp.md_pmode == "aniso")
    {
        if(mdp.md_pcouple == "xyz")
        {
            ModuleBase::WARNING_QUIT("Nose_Hoover", "md_pcouple==xyz will convert aniso to iso!");
        }
        pstart[0] = pstart[1] = pstart[2] = mdp.md_pfirst;
        pstop[0] = pstop[1] = pstop[2] = mdp.md_plast;
        pfreq[0] = pfreq[1] = pfreq[2] = mdp.md_pfreq;
        pflag[0] = pflag[1] = pflag[2] = 1;
    }
    else if(mdp.md_pmode == "tri")
    {
        pstart[0] = pstart[1] = pstart[2] = mdp.md_pfirst;
        pstop[0] = pstop[1] = pstop[2] = mdp.md_plast;
        pfreq[0] = pfreq[1] = pfreq[2] = mdp.md_pfreq;
        pflag[0] = pflag[1] = pflag[2] = 1;

        pstart[3] = pstart[4] = pstart[5] = 0;
        pstop[3] = pstop[4] = pstop[5] = 0;
        pfreq[3] = pfreq[4] = pfreq[5] = mdp.md_pfreq;
        pflag[3] = pflag[4] = pflag[5] = 1;
    }
    else if(mdp.md_pmode != "none")
    {
        ModuleBase::WARNING_QUIT("Nose_Hoover", "No such md_pmode yet!");
    }

    // determine whether NPT ensemble
    npt_flag = 0;
    for(int i=0; i<6; ++i)
    {
        npt_flag += pflag[i];
    }
    pdim = pflag[0] + pflag[1] + pflag[2];


    // allocate thermostats coupled with particles
    mass_eta = new double [mdp.md_tchain];
    eta = new double [mdp.md_tchain];
    v_eta = new double [mdp.md_tchain+1];
    g_eta = new double [mdp.md_tchain];

    v_eta[mdp.md_tchain] = 0;
    for(int i=0; i<mdp.md_tchain; ++i)
    {
        eta[i] = v_eta[i] = g_eta[i] = 0;
    }

    // allocate barostat and thermostats coupled with barostat
    if(npt_flag)
    {
        for(int i=0; i<6; ++i)
        {
            omega[i] = v_omega[i] = mass_omega[i] = 0;
        }

        if(mdp.md_pchain)
        {
            mass_peta = new double [mdp.md_pchain];
            peta = new double [mdp.md_pchain];
            v_peta = new double [mdp.md_pchain+1];
            g_peta = new double [mdp.md_pchain];

            v_peta[mdp.md_pchain] = 0;
            for(int i=0; i<mdp.md_pchain; ++i)
            {
                peta[i] = v_peta[i] = g_peta[i] = 0;
            }
        }
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
    delete []mass_eta;
    delete []eta;
    delete []v_eta;
    delete []g_eta;

    if(npt_flag && mdp.md_pchain)
    {
        delete []mass_peta;
        delete []peta;
        delete []v_peta;
        delete []g_peta;
    }
}

void Nose_Hoover::setup(ModuleESolver::ESolver *p_ensolve)
{
    ModuleBase::TITLE("Nose_Hoover", "setup");
    ModuleBase::timer::tick("Nose_Hoover", "setup");

    MDrun::setup(p_ensolve);

    // determine target temperature
    t_target = MD_func::target_temp(step_ + step_rst_, mdp.md_tfirst, mdp.md_tlast);

    // init thermostats coupled with particles
    mass_eta[0] = (3*ucell.nat - frozen_freedom_) * t_target / mdp.md_tfreq / mdp.md_tfreq;
    for(int m=1; m<mdp.md_tchain; ++m)
    {
        mass_eta[m] = t_target / mdp.md_tfreq / mdp.md_tfreq;
        g_eta[m] = (mass_eta[m-1]*v_eta[m-1]*v_eta[m-1]-t_target) / mass_eta[m];
    }

    // NPT ensemble
    if(npt_flag)  
    {
        // determine target stress 
        target_stress();

        // couple stress component due to md_pcouple
        couple_stress();

        // init barostat
        double nkt = (ucell.nat + 1) * t_target;

        for(int i=0; i<6; ++i)
        {
            if(pflag[i])
            {
                mass_omega[i] = nkt / pfreq[i] / pfreq[i];
            }
        }

        // init thermostats coupled with barostat
        if(mdp.md_pchain)
        {
            mass_peta[0] = t_target / mdp.md_pfreq / mdp.md_pfreq;
            for(int m=1; m<mdp.md_tchain; ++m)
            {
                mass_peta[m] = t_target / mdp.md_pfreq / mdp.md_pfreq;
                g_peta[m] = (mass_peta[m-1]*v_peta[m-1]*v_peta[m-1]-t_target) / mass_peta[m];
            }
        }
    }

    ModuleBase::timer::tick("Nose_Hoover", "setup");
}

void Nose_Hoover::first_half()
{
    ModuleBase::TITLE("Nose_Hoover", "first_half");
    ModuleBase::timer::tick("Nose_Hoover", "first_half");

    // update thermostats coupled with barostat if NPT ensemble
    if(npt_flag && mdp.md_tchain)
    {
        stress_integrate();
    }

    // update target T
    t_target = MD_func::target_temp(step_ + step_rst_, mdp.md_tfirst, mdp.md_tlast);

    // update thermostats coupled with particles
    temp_integrate();

    if(npt_flag)
    {
        // update temperature and stress due to velocity rescaling
        t_current = MD_func::current_temp(kinetic, ucell.nat, frozen_freedom_, allmass, vel);
        MD_func::compute_stress(ucell, vel, allmass, virial, stress);

        // couple stress component due to md_pcouple
        couple_stress();

        // determine target stress 
        target_stress();


    }



    MDrun::first_half();

    ModuleBase::timer::tick("Nose_Hoover", "first_half");
}

void Nose_Hoover::second_half()
{
    ModuleBase::TITLE("Nose_Hoover", "second_half");
    ModuleBase::timer::tick("Nose_Hoover", "second_half");

    MDrun::second_half();

    temp_integrate();

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
            file << v_eta[i] << "   ";
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
                    file >> v_eta[i];
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
    MPI_Bcast(v_eta, mdp.md_tchain, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
}

void Nose_Hoover::temp_integrate()
{
    t_current = MD_func::current_temp(kinetic, ucell.nat, frozen_freedom_, allmass, vel);

    // update mass_eta
    mass_eta[0] = (3*ucell.nat - frozen_freedom_) * t_target / mdp.md_tfreq / mdp.md_tfreq;
    for(int m=1; m<mdp.md_tchain; ++m)
    {
        mass_eta[m] = t_target / mdp.md_tfreq / mdp.md_tfreq;
    }

    // propogate g_eta
    if(mass_eta[0] > 0) 
    {
        g_eta[0] = (2*kinetic - (3*ucell.nat - frozen_freedom_)*t_target) / mass_eta[0];
    }
    else 
    {
        g_eta[0] = 0;
    }

    // integrate loop
    double factor;
    double scale = 1.0;
    double KE = kinetic;
    for(int i=0; i<nc_tchain; ++i)
    {
        for(int j=0; j<nys; ++j)
        {
            double delta = w[j] * mdp.md_dt / nc_tchain;

            // propogate v_eta
            for(int m=mdp.md_tchain-1; m>=0; --m)
            {
                factor = exp(-v_eta[m+1] * delta / 8.0);
                v_eta[m] *= factor;
                v_eta[m] += g_eta[m] * delta /4.0;
                v_eta[m] *= factor;
            }

            // propogate eta
            for(int m=0; m<mdp.md_tchain; ++m)
            {
                eta[m] += v_eta[m] * delta / 2.0;
            }

            // update rescale factor of particle velocity
            scale *= exp(-v_eta[0] * delta / 2.0);
            if(!isfinite(scale))
            {
                ModuleBase::WARNING_QUIT("Nose_Hoover", "Please set a proper md_tfreq!");
            }
            KE = kinetic * scale * scale;

            // propogate g_eta
            if(mass_eta[0] > 0) 
            {
                g_eta[0] = (2*KE - (3*ucell.nat - frozen_freedom_)*t_target) / mass_eta[0];
            }
            else 
            {
                g_eta[0] = 0;
            }

            // propogate v_eta
            v_eta[0] *= factor;
            v_eta[0] += g_eta[0] * delta /4.0;
            v_eta[0] *= factor;

            for(int m=1; m<mdp.md_tchain; ++m)
            {
                factor = exp(-v_eta[m+1] * delta / 8.0);
                v_eta[m] *= factor;
                g_eta[m] = (mass_eta[m-1] * v_eta[m-1] * v_eta[m-1] - t_target) / mass_eta[m];
                v_eta[m] += g_eta[m] * delta / 4.0;
                v_eta[m] *= factor;
            }
        }
    }

    // rescale velocity due to thermostats
    for(int i=0; i<ucell.nat; ++i)
    {
        vel[i] *= scale;
    }
}

void Nose_Hoover::stress_integrate()
{
    // the freedom of lattice
    int pdof = npt_flag;

    // update kenetic energy of lattice
    double ke_omega = 0;
    for(int i=0; i<6; ++i)
    {
        if(pflag[i])
        {
            ke_omega += mass_omega[i] * v_omega[i] * v_omega[i];
        }
    }

    // update force
    double lkt_press = t_target;
    if(mdp.md_pmode != "iso")
    {
        lkt_press *= pdof;
    }
    g_peta[0] = (ke_omega - lkt_press) / mass_peta[0];

    // integrate loop
    double factor;
    double scale = 1.0;
    double kecurrent = ke_omega;
    for(int i=0; i<nc_pchain; ++i)
    {
        for(int j=0; j<nys; ++j)
        {
            double delta = w[j] * mdp.md_dt / nc_pchain;

            // propogate v_peta
            for(int m=mdp.md_pchain-1; m>=0; --m)
            {
                factor = exp(-v_peta[m+1] * delta / 8.0);
                v_peta[m] *= factor;
                v_peta[m] += g_peta[m] * delta /4.0;
                v_peta[m] *= factor;
            }

            // propogate peta
            for(int m=0; m<mdp.md_pchain; ++m)
            {
                peta[m] += v_peta[m] * delta / 2.0;
            }

            // update rescale factor of lattice velocity
            scale *= exp(-v_peta[0] * delta / 2.0);
            kecurrent = ke_omega * scale * scale;

            // propogate g_peta
            g_peta[0] = (kecurrent - lkt_press) / mass_peta[0];

            // propogate v_peta
            v_peta[0] *= factor;
            v_peta[0] += g_peta[0] * delta /4.0;
            v_peta[0] *= factor;

            for(int m=1; m<mdp.md_pchain; ++m)
            {
                factor = exp(-v_peta[m+1] * delta / 8.0);
                v_peta[m] *= factor;
                g_peta[m] = (mass_peta[m-1] * v_peta[m-1] * v_peta[m-1] - t_target) / mass_peta[m];
                v_peta[m] += g_eta[m] * delta / 4.0;
                v_peta[m] *= factor;
            }
        }
    }

    // rescale lattice due to thermostats
    for(int i=0; i<6; ++i)
    {
        if(pflag[i])
        {
            v_omega[i] *= scale;
        }
    }
}

void Nose_Hoover::omega_integrate()
{
    double term_one = 0;


}


void Nose_Hoover::target_stress()
{
    double delta = (double)(step_ + step_rst_) / GlobalV::MD_NSTEP;

    p_hydro = 0;
    for(int i=0; i<3; ++i)
    {
        if(pflag[i])
        {
            p_target[i] = pstart[i] + delta * (pstop[i] - pstart[i]);
            p_hydro += p_target[i];
        }
    }
    if(pdim)
    {
        p_hydro /= pdim;
    }

    for(int i=3; i<6; ++i)
    {
        if(pflag[i])
        {
            p_target[i] = pstart[i] + delta * (pstop[i] - pstart[i]);
        }
    }
}

void Nose_Hoover::couple_stress()
{
    if(mdp.md_pcouple == "xyz")
    {
        double ave = (stress(0,0) + stress(1,1) + stress(2,2)) / 3.0;
        p_current[0] = p_current[1] = p_current[2] = ave;
    }
    else if(mdp.md_pcouple == "xy")
    {
        double ave = (stress(0,0) + stress(1,1)) / 2.0;
        p_current[0] = p_current[1] = ave;
        p_current[2] = stress(2,2);
    }
    else if(mdp.md_pcouple == "yz")
    {
        double ave = (stress(1,1) + stress(2,2)) / 2.0;
        p_current[1] = p_current[2] = ave;
        p_current[0] = stress(0,0);
    }
    else if(mdp.md_pcouple == "xz")
    {
        double ave = (stress(0,0) + stress(2,2)) / 2.0;
        p_current[0] = p_current[2] = ave;
        p_current[1] = stress(1,1);
    }
    else
    {
        p_current[0] = stress(0,0);
        p_current[1] = stress(1,1);
        p_current[2] = stress(2,2);
    }

    if(!isfinite(p_current[0]) || !isfinite(p_current[1]) || !isfinite(p_current[2]))
    {
        ModuleBase::WARNING_QUIT("Nose_Hoover", "Non-numeric stress component!");
    }

    p_current[3] = stress(1,2);
    p_current[4] = stress(0,2);
    p_current[5] = stress(0,1);

    if(!isfinite(p_current[3]) || !isfinite(p_current[4]) || !isfinite(p_current[5]))
    {
        ModuleBase::WARNING_QUIT("Nose_Hoover", "Non-numeric stress component!");
    }
}