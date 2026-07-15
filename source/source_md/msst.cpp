#include "msst.h"

#include "md_func.h"
#ifdef __MPI
#include "mpi.h"
#include "source_base/parallel_reduce.h"
#endif
#include "source_base/timer.h"

MSST::MSST(const Parameter& param_in, MdCell& mdcell_in) : MD_base(param_in, mdcell_in)
{
    msst_qmass = mdp.msst_qmass / pow(ModuleBase::ANGSTROM_AU, 4) / pow(ModuleBase::AU_to_MASS, 2);
    msst_vel = mdp.msst_vel * ModuleBase::ANGSTROM_AU * ModuleBase::AU_to_FS;
    msst_vis = mdp.msst_vis / ModuleBase::AU_to_MASS / ModuleBase::ANGSTROM_AU * ModuleBase::AU_to_FS;

    assert(mdcell.nat() > 0);

    dilation.set(1, 1, 1);
    omega.set(0, 0, 0);
    p0 = 0;
    e0 = 0;
    v0 = 1;
    totmass = 0;
    lag_pos = 0;
    vsum = 0;

    for (int i = 0; i < mdcell.nlocal(); ++i)
    {
        totmass += allmass[i];
    }
#ifdef __MPI
    Parallel_Reduce::reduce_all(&totmass, 1);
#endif
}

MSST::~MSST()
{
}

void MSST::setup(ModuleESolver::ESolver* p_esolver, const std::string& global_readin_dir)
{
    ModuleBase::TITLE("MSST", "setup");
    ModuleBase::timer::start("MSST", "setup");

    MD_base::setup(p_esolver, global_readin_dir);

    int sd = mdp.msst_direction;

    if (!mdp.md_restart)
    {
        lag_pos = 0;
        v0 = mdcell.omega();
        p0 = stress(sd, sd);
        e0 = potential + kinetic;

        if (kinetic > 0 && mdp.msst_tscale > 0)
        {
            double fac1 = mdp.msst_tscale * totmass * 2.0 * kinetic / msst_qmass;
            omega[sd] = -1.0 * sqrt(fac1);
            double fac2 = omega[sd] / v0;

            std::cout << "initial strain rate = " << fac2 << "    msst_tscale = " << mdp.msst_tscale << std::endl;

            for (int i = 0; i < state_.size(); ++i)
            {
                vel[i] *= sqrt(1.0 - mdp.msst_tscale);
            }
            sync_velocity_buffer_to_state();
        }

        MD_func::compute_stress(mdcell, vel, allmass, cal_stress, virial, stress);
        t_current = MD_func::current_temp(kinetic, mdcell, frozen_freedom_, allmass, vel);
    }

    ModuleBase::timer::end("MSST", "setup");

    return;
}

void MSST::first_half(std::ofstream& ofs)
{
    ModuleBase::TITLE("MSST", "first_half");
    ModuleBase::timer::start("MSST", "first_half");

    const int sd = mdp.msst_direction;
    const double dthalf = 0.5 * md_dt;
    double vol = 0.0;
    energy_ = potential + kinetic;

    /// propagate the time derivative of volume 1/2 step
    propagate_voldot();

    vsum = vel_sum();

    /// save the velocities
    old_v.resize(static_cast<std::size_t>(state_.size()));
    for (int i = 0; i < state_.size(); ++i)
    {
        old_v[static_cast<std::size_t>(i)] = vel[i];
    }

    /// propagate velocity sum 1/2 step by temporarily propagating the velocities
    propagate_vel();

    vsum = vel_sum();

    /// reset the velocities
    for (int i = 0; i < state_.size(); ++i)
    {
        vel[i] = old_v[static_cast<std::size_t>(i)];
    }
    sync_velocity_buffer_to_state();

    /// propagate velocities 1/2 step using the new velocity sum
    propagate_vel();

    /// propagate volume 1/2 step
    vol = mdcell.omega() + omega[sd] * dthalf;

    /// rescale positions and change box size
    rescale(ofs, vol);

    /// propagate atom positions 1 time step
    MD_base::update_pos();

    /// propagate volume 1/2 step
    vol = mdcell.omega() + omega[sd] * dthalf;

    /// rescale positions and change box size
    rescale(ofs, vol);

    ModuleBase::timer::end("MSST", "first_half");

    return;
}

void MSST::second_half()
{
    ModuleBase::TITLE("MSST", "second_half");
    ModuleBase::timer::start("MSST", "second_half");

    const int sd = mdp.msst_direction;
    const double dthalf = 0.5 * md_dt;
    energy_ = potential + kinetic;

    /// propagate velocities 1/2 step
    propagate_vel();

    vsum = vel_sum();
    MD_func::compute_stress(mdcell, vel, allmass, cal_stress, virial, stress);
    t_current = MD_func::current_temp(kinetic, mdcell, frozen_freedom_, allmass, vel);

    /// propagate the time derivative of volume 1/2 step
    propagate_voldot();

    /// calculate Lagrangian position
    lag_pos -= msst_vel * mdcell.omega() / v0 * md_dt;

    ModuleBase::timer::end("MSST", "second_half");

    return;
}


void MSST::print_md(std::ofstream& ofs, const bool& cal_stress)
{
    MD_base::print_md(ofs, cal_stress);

    return;
}

void MSST::write_restart(const std::string& global_out_dir)
{
    if (!my_rank)
    {
        std::stringstream ssc;
        ssc << global_out_dir << "Restart_md.txt";
        std::ofstream file(ssc.str().c_str());

        file << step_ + step_rst_ << std::endl;
        file << md_tfirst << std::endl;
        file << omega[mdp.msst_direction] << std::endl;
        file << e0 << std::endl;
        file << v0 << std::endl;
        file << p0 << std::endl;
        file << lag_pos << std::endl;

        file.close();
    }
#ifdef __MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    return;
}


void MSST::restart(const std::string& global_readin_dir)
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
            file >> step_rst_ >> md_tfirst >> omega[mdp.msst_direction] >> e0 >> v0 >> p0 >> lag_pos;
            file.close();
        }
    }

#ifdef __MPI
    MPI_Bcast(&ok, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

    if (!ok)
    {
        ModuleBase::WARNING_QUIT("mdrun", "no Restart_md.txt !");
    }

#ifdef __MPI
    MPI_Bcast(&step_rst_, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&md_tfirst, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&omega[mdp.msst_direction], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&e0, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&v0, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&p0, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&lag_pos, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

    return;
}

double MSST::vel_sum() const
{
    double vsum = 0;

    for (int i = 0; i < state_.size(); ++i)
    {
        vsum += vel[i].norm2();
    }
#ifdef __MPI
    Parallel_Reduce::reduce_all(&vsum, 1);
#endif
    return vsum;
}

void MSST::rescale(std::ofstream& ofs, const double& volume)
{
    int sd = mdp.msst_direction;

    assert(mdcell.omega() > 0.0);

    dilation[sd] = volume / mdcell.omega();
    ModuleBase::Matrix3 latvec = mdcell.latvec();
    latvec.e11 *= dilation[0];
    latvec.e22 *= dilation[1];
    latvec.e33 *= dilation[2];
    mdcell.set_lattice_vectors(latvec);
    mdcell.refresh_cart_from_frac();

    /// rescale velocity
    for (int i = 0; i < state_.size(); ++i)
    {
        vel[i][sd] *= dilation[sd];
    }
    sync_velocity_buffer_to_state();
    refresh_runtime_storage_from_mdcell();
    static_cast<void>(ofs);
}


void MSST::propagate_vel()
{
    const int sd = mdp.msst_direction;
    const double dthalf = 0.5 * md_dt;
    const double fac = msst_vis * pow(omega[sd], 2) / (vsum * mdcell.omega());

    for (int i = 0; i < state_.size(); ++i)
    {
        ModuleBase::Vector3<double> const_C = force[i] / allmass[i];
        ModuleBase::Vector3<double> const_D;
        const_D.set(fac / allmass[i], fac / allmass[i], fac / allmass[i]);
        const_D[sd] -= 2 * omega[sd] / mdcell.omega();

        for (int k = 0; k < 3; ++k)
        {
            if (fabs(dthalf * const_D[k]) > 1e-6)
            {
                double expd = exp(dthalf * const_D[k]);
                vel[i][k] = expd * (const_C[k] + const_D[k] * vel[i][k] - const_C[k] / expd) / const_D[k];
            }
            else
            {
                vel[i][k]
                    += (const_C[k] + const_D[k] * vel[i][k]) * dthalf
                       + 0.5 * (const_D[k] * const_D[k] * vel[i][k] + const_C[k] * const_D[k]) * dthalf * dthalf;
            }
        }
    }
    sync_velocity_buffer_to_state();

    return;
}


void MSST::propagate_voldot()
{
    const int sd = mdp.msst_direction;
    const double dthalf = 0.5 * md_dt;
    double p_current = stress(sd, sd);
    double p_msst = msst_vel * msst_vel * totmass * (v0 - mdcell.omega()) / (v0 * v0);
    double const_A = totmass * (p_current - p0 - p_msst) / msst_qmass;
    double const_B = totmass * msst_vis / (msst_qmass * mdcell.omega());

    /// prevent the increase of volume
    if (mdcell.omega() > v0 && const_A > 0)
    {
        const_A = -const_A;
    }

    /// avoid singularity at B = 0 with Taylor expansion
    double fac = const_B * dthalf;
    if (fac > 1e-6)
    {
        omega[sd] = (omega[sd] + const_A * (exp(fac) - 1) / const_B) * exp(-fac);
    }
    else
    {
        omega[sd] += (const_A - const_B * omega[sd]) * dthalf
                     + 0.5 * (const_B * const_B * omega[sd] - const_A * const_B) * dthalf * dthalf;
    }

    return;
}
