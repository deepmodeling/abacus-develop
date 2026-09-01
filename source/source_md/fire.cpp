#include "fire.h"

#include "md_func.h"
#ifdef __MPI
#include "mpi.h"
#endif
#include "source_base/timer.h"

FIRE::FIRE(const Parameter& param_in, MDCell& mdcell_in) : MD_base(param_in, mdcell_in)
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

    // BUGFIX:
    // Do not override the force convergence threshold read from INPUT.
    // force_thr is stored in internal force unit, Hartree/Bohr.
    // force_thr = 1e-3;
}

FIRE::~FIRE()
{
}

void FIRE::setup(ModuleESolver::ESolver* p_esolver, const std::string& global_readin_dir)
{
    ModuleBase::TITLE("FIRE", "setup");
    ModuleBase::timer::start("FIRE", "setup");

    MD_base::setup(p_esolver, global_readin_dir);

    check_force();

    ModuleBase::timer::end("FIRE", "setup");

    return;
}

void FIRE::first_half(std::ofstream& ofs)
{
    ModuleBase::TITLE("FIRE", "first_half");
    ModuleBase::timer::start("FIRE", "first_half");

    MD_base::update_vel();

    check_fire();

    MD_base::update_pos();

    ModuleBase::timer::end("FIRE", "first_half");

    return;
}


void FIRE::second_half(void)
{
    ModuleBase::TITLE("FIRE", "second_half");
    ModuleBase::timer::start("FIRE", "second_half");

    MD_base::update_vel();

    check_force();

    ModuleBase::timer::end("FIRE", "second_half");

    return;
}


void FIRE::print_md(std::ofstream& ofs, const bool& cal_stress)
{
    MD_base::print_md(ofs, cal_stress);

    const double max_force = max * ModuleBase::Hartree_to_eV * ModuleBase::ANGSTROM_AU;

	ofs << " LARGEST FORCE (eV/A)      : " << max_force << std::endl;
	std::cout << " LARGEST FORCE (eV/A)  : " << max_force << std::endl;

    return;
}


void FIRE::write_restart(const std::string& global_out_dir)
{
    if (!my_rank)
    {
        std::stringstream ssc;
        ssc << global_out_dir << "Restart_md.txt";
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
    MPI_Barrier(mdcell.communicator());
#endif

    return;
}


void FIRE::restart(const std::string& global_readin_dir)
{
    bool ok = true;

    if (!my_rank)
    {
        std::stringstream ssc;
        ssc << global_readin_dir << "Restart_md.txt";
        std::ifstream file(ssc.str().c_str());

        if (!file)
        {
            ok = false;
        }

        if (ok)
        {
            file >> step_rst_ >> md_tfirst >> alpha >> negative_count >> dt_max >> md_dt;
            if(!file)
            {
                ok = false;
            }    
            file.close();
        }
    }

#ifdef __MPI
    MPI_Bcast(&ok, 1, MPI_C_BOOL, 0, mdcell.communicator());
#endif

    if (!ok)
    {
        ModuleBase::WARNING_QUIT("mdrun", "no Restart_md.txt !");
    }

#ifdef __MPI
    MPI_Bcast(&step_rst_, 1, MPI_INT, 0, mdcell.communicator());
    MPI_Bcast(&md_tfirst, 1, MPI_DOUBLE, 0, mdcell.communicator());
    MPI_Bcast(&alpha, 1, MPI_DOUBLE, 0, mdcell.communicator());
    MPI_Bcast(&negative_count, 1, MPI_INT, 0, mdcell.communicator());
    MPI_Bcast(&dt_max, 1, MPI_DOUBLE, 0, mdcell.communicator());
    MPI_Bcast(&md_dt, 1, MPI_DOUBLE, 0, mdcell.communicator());
#endif

    return;
}

void FIRE::check_force(void)
{
    max = 0.0;

    std::int64_t movable_dof = 0;

    for (const LocalAtom& atom : mdcell.owned_atoms())
    {
        for (int j = 0; j < 3; ++j)
        {
            // Only movable degrees of freedom should be used
            // in the FIRE convergence criterion.
            //
            // For example:
            // m 1 1 1 -> x/y/z are included.
            // m 1 0 1 -> y is excluded.
            // m 0 0 0 -> this atom contributes no DOF to convergence.
            if (!atom.mbl[j])
            {
                continue;
            }

            ++movable_dof;

            if (max < std::abs(atom.force[j]))
            {
                max = std::abs(atom.force[j]);
            }
        }
    }

 #ifdef __MPI
    if (mdcell.mpi_size() > 1)
    {
        double global_max = 0.0;
        std::int64_t global_movable_dof = 0;
        MPI_Allreduce(&max, &global_max, 1, MPI_DOUBLE, MPI_MAX, mdcell.communicator());
        MPI_Allreduce(&movable_dof, &global_movable_dof, 1, MPI_INT64_T, MPI_SUM, mdcell.communicator());
        max = global_max;
        movable_dof = global_movable_dof;
    }
#endif

    // If there are no movable degrees of freedom, there is nothing to optimize.
    if (movable_dof == 0)
    {
        stop = true;
        return;
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

    std::int64_t movable_dof = 0;

    // Compute P, |F| and |v| only on movable degrees of freedom.
    // Fixed atoms/directions may have non-zero raw forces, but they should not
    // affect the FIRE velocity projection or adaptive time-step control.
    for (LocalAtom& atom : mdcell.mutable_owned_atoms())
    {
        for (int j = 0; j < 3; ++j)
        {
            if (!atom.mbl[j])
            {
                // Keep frozen components clean.
                atom.vel[j] = 0.0;
                continue;
            }

            ++movable_dof;

            P += atom.vel[j] * atom.force[j];
            sumforce += atom.force[j] * atom.force[j];
            normvel += atom.vel[j] * atom.vel[j];
        }
    }

 #ifdef __MPI
    if (mdcell.mpi_size() > 1)
    {
        double local_values[3] = {P, sumforce, normvel};
        double global_values[3] = {0.0, 0.0, 0.0};
        std::int64_t global_movable_dof = 0;
        MPI_Allreduce(local_values, global_values, 3, MPI_DOUBLE, MPI_SUM, mdcell.communicator());
        MPI_Allreduce(&movable_dof, &global_movable_dof, 1, MPI_INT64_T, MPI_SUM, mdcell.communicator());
        P = global_values[0];
        sumforce = global_values[1];
        normvel = global_values[2];
        movable_dof = global_movable_dof;
    }
#endif

    // No movable degrees of freedom: nothing to update.
    if (movable_dof == 0)
    {
        return;
    }

    sumforce = std::sqrt(sumforce);
    normvel = std::sqrt(normvel);

    // If force or velocity norm is zero, the velocity projection is undefined.
    // Avoid 0/0. In a truly converged case check_force() should stop the run.
    if (sumforce > 0.0 && normvel > 0.0)
    {
        for (LocalAtom& atom : mdcell.mutable_owned_atoms())
        {
            for (int j = 0; j < 3; ++j)
            {
                if (!atom.mbl[j])
                {
                    atom.vel[j] = 0.0;
                    continue;
                }

                atom.vel[j] = (1.0 - alpha) * atom.vel[j]
                            + alpha * atom.force[j] / sumforce * normvel;
            }
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

        for (LocalAtom& atom : mdcell.mutable_owned_atoms())
        {
            for (int j = 0; j < 3; ++j)
            {
                atom.vel[j] = 0;
            }
        }

        alpha = alpha_start;
    }

    return;
}
