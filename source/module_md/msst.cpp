#include "msst.h"

#include "md_func.h"
#ifdef __MPI
#include "mpi.h"
#endif
#include "module_base/timer.h"

MSST::MSST(MD_para& MD_para_in, UnitCell& unit_in) : MD_base(MD_para_in, unit_in)
{
    mdp.msst_qmass = mdp.msst_qmass / pow(ModuleBase::ANGSTROM_AU, 4) / pow(ModuleBase::AU_to_MASS, 2);
    mdp.msst_vel = mdp.msst_vel * ModuleBase::ANGSTROM_AU * ModuleBase::AU_to_FS;
    mdp.msst_vis = mdp.msst_vis / ModuleBase::AU_to_MASS / ModuleBase::ANGSTROM_AU * ModuleBase::AU_to_FS;

    old_v = new ModuleBase::Vector3<double>[ucell.nat];
    dilation.set(1, 1, 1);
    omega.set(0, 0, 0);
    p0 = 0;
    e0 = 0;
    v0 = 1;
    totmass = 0;

    for (int i = 0; i < ucell.nat; ++i)
    {
        totmass += allmass[i];
    }
}

MSST::~MSST()
{
    delete[] old_v;
}

void MSST::setup(ModuleESolver::ESolver* p_esolver, const std::string& global_readin_dir)
{
    ModuleBase::TITLE("MSST", "setup");
    ModuleBase::timer::tick("MSST", "setup");

    MD_base::setup(p_esolver, global_readin_dir);
    ucell.cell_parameter_updated = true;

    int sd = mdp.msst_direction;

    if (!mdp.md_restart)
    {
        lag_pos = 0;
        v0 = ucell.omega;
        p0 = stress(sd, sd);
        e0 = potential + kinetic;

        if (kinetic > 0 && mdp.msst_tscale > 0)
        {
            double fac1 = mdp.msst_tscale * totmass * 2.0 * kinetic / mdp.msst_qmass;
            omega[sd] = -1.0 * sqrt(fac1);
            double fac2 = omega[sd] / v0;

            std::cout << "initial strain rate = " << fac2 << "    msst_tscale = " << mdp.msst_tscale << std::endl;

            for (int i = 0; i < ucell.nat; ++i)
            {
                vel[i] *= sqrt(1.0 - mdp.msst_tscale);
            }
        }

        MD_func::compute_stress(ucell, vel, allmass, mdp.cal_stress, virial, stress);
        t_current = MD_func::current_temp(kinetic, ucell.nat, frozen_freedom_, allmass, vel);
    }

    ModuleBase::timer::tick("MSST", "setup");
}

void MSST::first_half(std::ofstream& ofs)
{
    ModuleBase::TITLE("MSST", "first_half");
    ModuleBase::timer::tick("MSST", "first_half");

    const int sd = mdp.msst_direction;
    const double dthalf = 0.5 * mdp.md_dt;
    double vol;
    energy_ = potential + kinetic;

    /// propagate the time derivative of volume 1/2 step
    propagate_voldot();

    vsum = vel_sum();

    /// save the velocities
    for (int i = 0; i < ucell.nat; ++i)
    {
        old_v[i] = vel[i];
    }

    /// propagate velocity sum 1/2 step by temporarily propagating the velocities
    propagate_vel();

    vsum = vel_sum();

    /// reset the velocities
    for (int i = 0; i < ucell.nat; ++i)
    {
        vel[i] = old_v[i];
    }

    /// propagate velocities 1/2 step using the new velocity sum
    propagate_vel();

    /// propagate volume 1/2 step
    vol = ucell.omega + omega[sd] * dthalf;

    /// rescale positions and change box size
    rescale(ofs, vol);

    /// propagate atom positions 1 time step
    MD_base::update_pos();

    /// propagate volume 1/2 step
    vol = ucell.omega + omega[sd] * dthalf;

    /// rescale positions and change box size
    rescale(ofs, vol);

    ModuleBase::timer::tick("MSST", "first_half");
}

void MSST::second_half()
{
    ModuleBase::TITLE("MSST", "second_half");
    ModuleBase::timer::tick("MSST", "second_half");

    const int sd = mdp.msst_direction;
    const double dthalf = 0.5 * mdp.md_dt;
    energy_ = potential + kinetic;

    /// propagate velocities 1/2 step
    propagate_vel();

    vsum = vel_sum();
    MD_func::compute_stress(ucell, vel, allmass, mdp.cal_stress, virial, stress);
    t_current = MD_func::current_temp(kinetic, ucell.nat, frozen_freedom_, allmass, vel);

    /// propagate the time derivative of volume 1/2 step
    propagate_voldot();

    /// calculate Lagrangian position
    lag_pos -= mdp.msst_vel * ucell.omega / v0 * mdp.md_dt;

    ModuleBase::timer::tick("MSST", "second_half");
}

void MSST::print_md(std::ofstream& ofs, const bool& cal_stress)
{
    MD_base::print_md(ofs, cal_stress);
}

void MSST::write_restart(const std::string& global_out_dir)
{
    if (!mdp.my_rank)
    {
        std::stringstream ssc;
        ssc << global_out_dir << "Restart_md.dat";
        std::ofstream file(ssc.str().c_str());

        file << step_ + step_rst_ << std::endl;
        file << mdp.md_tfirst << std::endl;
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
}

void MSST::restart(const std::string& global_readin_dir)
{
    bool ok = true;

    if (!mdp.my_rank)
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
            file >> step_rst_ >> mdp.md_tfirst >> omega[mdp.msst_direction] >> e0 >> v0 >> p0 >> lag_pos;
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
    MPI_Bcast(&mdp.md_tfirst, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&omega[mdp.msst_direction], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&e0, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&v0, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&p0, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&lag_pos, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
}

double MSST::vel_sum()
{
    double vsum = 0;

    for (int i = 0; i < ucell.nat; ++i)
    {
        vsum += vel[i].norm2();
    }

    return vsum;
}

void MSST::rescale(std::ofstream& ofs, const double& volume)
{
    int sd = mdp.msst_direction;

    dilation[sd] = volume / ucell.omega;
    ucell.latvec.e11 *= dilation[0];
    ucell.latvec.e22 *= dilation[1];
    ucell.latvec.e33 *= dilation[2];

    ucell.setup_cell_after_vc(ofs);

    /// rescale velocity
    for (int i = 0; i < ucell.nat; ++i)
    {
        vel[i][sd] *= dilation[sd];
    }
}

void MSST::propagate_vel()
{
    if (mdp.my_rank == 0)
    {
        const int sd = mdp.msst_direction;
        const double dthalf = 0.5 * mdp.md_dt;
        const double fac = mdp.msst_vis * pow(omega[sd], 2) / (vsum * ucell.omega);

        for (int i = 0; i < ucell.nat; ++i)
        {
            ModuleBase::Vector3<double> const_C = force[i] / allmass[i];
            ModuleBase::Vector3<double> const_D;
            const_D.set(fac / allmass[i], fac / allmass[i], fac / allmass[i]);
            const_D[sd] -= 2 * omega[sd] / ucell.omega;

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
    }

#ifdef __MPI
    MPI_Bcast(vel, ucell.nat * 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
}

void MSST::propagate_voldot()
{
    const int sd = mdp.msst_direction;
    const double dthalf = 0.5 * mdp.md_dt;
    double p_current = stress(sd, sd);
    double p_msst = mdp.msst_vel * mdp.msst_vel * totmass * (v0 - ucell.omega) / (v0 * v0);
    double const_A = totmass * (p_current - p0 - p_msst) / mdp.msst_qmass;
    double const_B = totmass * mdp.msst_vis / (mdp.msst_qmass * ucell.omega);

    /// prevent the increase of volume
    if (ucell.omega > v0 && const_A > 0)
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
}
