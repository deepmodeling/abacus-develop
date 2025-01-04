#include "fire.h"

#include "md_func.h"
#ifdef __MPI
#include "mpi.h"
#endif
#include "module_base/timer.h"

FIRE::FIRE(const Parameter& param_in, UnitCell& unit_in) : MD_base(param_in, unit_in)
{
    force_thr = param_in.inp.force_thr;
    dt_max = -1.0;
    alpha_start = 0.10;
    alpha = alpha_start;

    finc = 1.1;
    fdec = 0.5;
    f_alpha = 0.99;
    n_min = 4;
    negative_count = 0;
    max = 0.0;
    force_thr = 1e-3;
}

FIRE::~FIRE()
{
}

void FIRE::setup(ModuleESolver::ESolver* p_esolver, const std::string& global_readin_dir)
{
    ModuleBase::TITLE("FIRE", "setup");
    ModuleBase::timer::tick("FIRE", "setup");

    MD_base::setup(p_esolver, global_readin_dir);

    check_force();

    ModuleBase::timer::tick("FIRE", "setup");

    return;
}

void FIRE::first_half(std::ofstream& ofs)
{
    ModuleBase::TITLE("FIRE", "first_half");
    ModuleBase::timer::tick("FIRE", "first_half");

    MD_base::update_vel(force);

    check_fire();

    MD_base::update_pos();

    ModuleBase::timer::tick("FIRE", "first_half");

    return;
}


void FIRE::second_half(void)
{
    ModuleBase::TITLE("FIRE", "second_half");
    ModuleBase::timer::tick("FIRE", "second_half");

    MD_base::update_vel(force);

    check_force();

    ModuleBase::timer::tick("FIRE", "second_half");

    return;
}


void FIRE::print_md(std::ofstream& ofs, const bool& cal_stress)
{
    MD_base::print_md(ofs, cal_stress);

    ofs << "\n Largest gradient in force is " << max * ModuleBase::Hartree_to_eV * ModuleBase::ANGSTROM_AU << " eV/A." << std::endl;
    ofs << " Threshold is " << PARAM.inp.force_thr_ev << " eV/A." << std::endl;
    std::cout << " LARGEST GRAD (eV/A)  : " << max * ModuleBase::Hartree_to_eV * ModuleBase::ANGSTROM_AU << std::endl;

    return;
}


void FIRE::write_restart(const std::string& global_out_dir)
{
    if (!my_rank)
    {
        std::stringstream ssc;
        ssc << global_out_dir << "Restart_md.dat";
        std::ofstream file(ssc.str().c_str());

        file << step_ + step_rst_ << std::endl;
        file << md_tfirst << std::endl;
        file << alpha << std::endl;
        file << negative_count << std::endl;
        file << dt_max << std::endl;
        file << md_dt << std::endl;
        file.close();
    }
#ifdef __MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    return;
}


void FIRE::restart(const std::string& global_readin_dir)
{
    bool ok = true;

    if (!my_rank)
    {
        std::stringstream ssc;
        ssc << global_readin_dir << "Restart_md.dat";
        std::ifstream file(ssc.str().c_str());

        if (!file)
        {
            ok = false;
        }

        if (ok)
        {
            file >> step_rst_ >> md_tfirst >> alpha >> negative_count >> dt_max >> md_dt;
            file.close();
        }
    }

#ifdef __MPI
    MPI_Bcast(&ok, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

    if (!ok)
    {
        ModuleBase::WARNING_QUIT("mdrun", "no Restart_md.dat !");
    }

#ifdef __MPI
    MPI_Bcast(&step_rst_, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&md_tfirst, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&alpha, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&negative_count, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&dt_max, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&md_dt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

    return;
}


void FIRE::check_force(void)
{
    max = 0;

    for (int i = 0; i < ucell.nat; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            if (max < std::abs(force[i][j]))
            {
                max = std::abs(force[i][j]);
            }
        }
    }

    if (2.0 * max < force_thr)
    {
        stop = true;
    }

    return;
}


void FIRE::check_fire(void)
{
    double P = 0.0;
    double sumforce = 0.0;
    double normvel = 0.0;

    /// initial dt_max
    if (dt_max < 0)
    {
        dt_max = 2.5 * md_dt;
    }

    for (int i = 0; i < ucell.nat; ++i)
    {
        P += vel[i].x * force[i].x + vel[i].y * force[i].y + vel[i].z * force[i].z;
        sumforce += force[i].norm2();
        normvel += vel[i].norm2();
    }

    sumforce = sqrt(sumforce);
    normvel = sqrt(normvel);

    for (int i = 0; i < ucell.nat; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            vel[i][j] = (1.0 - alpha) * vel[i][j] + alpha * force[i][j] / sumforce * normvel;
        }
    }

    if (P > 0)
    {
        negative_count++;
        if (negative_count >= n_min)
        {
            md_dt = std::min(md_dt * finc, dt_max);
            alpha *= f_alpha;
        }
    }
    else
    {
        md_dt *= fdec;
        negative_count = 0;

        for (int i = 0; i < ucell.nat; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                vel[i][j] = 0;
            }
        }

        alpha = alpha_start;
    }
    
    return;
}
