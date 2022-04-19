#include "MSST.h"
#include "MD_func.h"
#ifdef __MPI
#include "mpi.h"
#endif
#include "../module_base/timer.h"
#include "module_esolver/esolver.h"

MSST::MSST(MD_parameters& MD_para_in, UnitCell_pseudo &unit_in) : Verlet(MD_para_in, unit_in)
{
    std::cout << "MSST" << std::endl;

    mdp.msst_qmass = mdp.msst_qmass / pow(ModuleBase::ANGSTROM_AU, 4) / pow(ModuleBase::AU_to_MASS, 2);
    mdp.msst_vel = mdp.msst_vel * ModuleBase::ANGSTROM_AU * ModuleBase::AU_to_FS;
    mdp.msst_vis = mdp.msst_vis / ModuleBase::AU_to_MASS / ModuleBase::ANGSTROM_AU * ModuleBase::AU_to_FS;

    old_v = new ModuleBase::Vector3<double> [ucell.nat];
    dilation.set(1,1,1);
    omega.set(0,0,0);
    p0 = 0;
    e0 = 0;
    v0 = 1;
    totmass = 0;

    for(int i=0; i<ucell.nat; ++i)
    {
        totmass += allmass[i];
    }
}

MSST::~MSST()
{
    delete []old_v;
}

void MSST::setup(ModuleESolver::ESolver *p_esolver)
{
    ModuleBase::TITLE("MSST", "setup");
    ModuleBase::timer::tick("MSST", "setup");

    int sd = mdp.msst_direction;

    MD_func::force_virial(p_esolver, step_, mdp, ucell, potential, force, virial);
    MD_func::kinetic_stress(ucell, vel, allmass, kinetic, stress);
    stress += virial;

    if(mdp.md_restart)
    {
        restart();
    }
    else
    {
        lag_pos = 0;
        v0 = ucell.omega;
        p0 = stress(sd, sd);
        e0 = potential + kinetic;

        if(kinetic > 0 && mdp.msst_tscale > 0)
        {
            double fac1 = mdp.msst_tscale * totmass * 2.0 * kinetic / mdp.msst_qmass;
            omega[sd] = -1.0 * sqrt(fac1);
            double fac2 = omega[sd] / v0;

            std::cout << "initial strain rate = " << fac2 << "    msst_tscale = " << mdp.msst_tscale << std::endl;

            for(int i=0; i<ucell.nat; ++i)
            {
                vel[i] *= sqrt(1.0 - mdp.msst_tscale);
            }
        }

        MD_func::kinetic_stress(ucell, vel, allmass, kinetic, stress);
        stress += virial;
    }

    ModuleBase::timer::tick("MSST", "setup");
}

void MSST::first_half()
{
    ModuleBase::TITLE("MSST", "first_half");
    ModuleBase::timer::tick("MSST", "first_half");

    const int sd = mdp.msst_direction;
    const double dthalf = 0.5 * mdp.md_dt;

    energy_ = potential + kinetic;

    // propagate the time derivative of volume 1/2 step
    propagate_voldot();

    vsum = vel_sum();

    // save the velocities
    for(int i=0; i<ucell.nat; ++i)
    {
        old_v[i] = vel[i];
    }

    // propagate velocity sum 1/2 step by temporarily propagating the velocities
    propagate_vel();

    vsum = vel_sum();

    // reset the velocities
    for(int i=0; i<ucell.nat; ++i)
    {
        vel[i] = old_v[i];
    }

    // propagate velocities 1/2 step using the new velocity sum
    propagate_vel();

    // propagate volume 1/2 step
    double vol = ucell.omega + omega[sd] * dthalf;

    // rescale positions and change box size
    rescale(vol);

    // propagate atom positions 1 time step
    for(int i=0; i<ucell.nat; ++i)
    {
        pos[i] += vel[i] * mdp.md_dt;
    }
    ucell.update_pos_tau(pos);
    ucell.periodic_boundary_adjustment();

    // propagate volume 1/2 step
    vol = ucell.omega + omega[sd] * dthalf;

    // rescale positions and change box size
    rescale(vol);

    ModuleBase::timer::tick("MSST", "first_half");
}

void MSST::second_half()
{
    ModuleBase::TITLE("MSST", "second_half");
    ModuleBase::timer::tick("MSST", "second_half");

    const int sd = mdp.msst_direction;
    const double dthalf = 0.5 * mdp.md_dt;

    energy_ = potential + kinetic;

    // propagate velocities 1/2 step
    propagate_vel();

    vsum = vel_sum();
    MD_func::kinetic_stress(ucell, vel, allmass, kinetic, stress);
    stress += virial;

    // propagate the time derivative of volume 1/2 step
    propagate_voldot();

    // calculate Lagrangian position
    lag_pos -= mdp.msst_vel * ucell.omega / v0 * mdp.md_dt;

    ModuleBase::timer::tick("MSST", "second_half");
}

void MSST::outputMD(std::ofstream &ofs)
{
    Verlet::outputMD(ofs);
}

void MSST::write_restart()
{
    if(!GlobalV::MY_RANK)
    {
		std::stringstream ssc;
		ssc << GlobalV::global_out_dir << "Restart_md.dat";
		std::ofstream file(ssc.str().c_str());

        file << step_ + step_rst_ << std::endl;
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

void MSST::restart()
{
    if(!GlobalV::MY_RANK)
    {
		std::stringstream ssc;
		ssc << GlobalV::global_out_dir << "Restart_md.dat";
		std::ifstream file(ssc.str().c_str());

        if(!file)
		{
			std::cout<< "please ensure whether 'Restart_md.dat' exists!" << std::endl;
            ModuleBase::WARNING_QUIT("MSST", "no Restart_md.dat ！");
		}

		file >> step_rst_ >> omega[mdp.msst_direction] >> e0 >> v0 >> p0 >> lag_pos;

		file.close();
	}

#ifdef __MPI
	MPI_Bcast(&step_rst_, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&omega[mdp.msst_direction], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&e0, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&v0, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&p0, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&lag_pos, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
}

double MSST::extra_term()
{
    return 0;
}

double MSST::vel_sum()
{
    double vsum = 0;

    for(int i=0; i<ucell.nat; ++i)
    {
        vsum += vel[i].norm2();
    }

    return vsum;
}

void MSST::rescale(double volume)
{
    int sd = mdp.msst_direction;

    dilation[sd] = volume/ucell.omega;
    ucell.latvec.e11 *= dilation[0];
    ucell.latvec.e22 *= dilation[1];
    ucell.latvec.e33 *= dilation[2];

    ucell.setup_cell_after_vc(GlobalV::ofs_running);
    MD_func::InitPos(ucell, pos);

    // rescale velocity
    for(int i=0; i<ucell.nat; ++i)
    {
        vel[i][sd] *= dilation[sd];
    }
}

void MSST::propagate_vel()
{
    const int sd = mdp.msst_direction;
    const double dthalf = 0.5 * mdp.md_dt;
    const double fac = mdp.msst_vis * pow(omega[sd], 2) / (vsum * ucell.omega);

    for(int i=0; i<ucell.nat; ++i)
    {
        ModuleBase::Vector3<double> const_C = force[i] / allmass[i];
        ModuleBase::Vector3<double> const_D;
        const_D.set(fac/allmass[i], fac/allmass[i], fac/allmass[i]);
        const_D[sd] -= 2 * omega[sd] / ucell.omega;

        for(int k=0; k<3; ++k)
        {
            if( fabs(dthalf*const_D[k]) > 1e-6 )
            {
                double expd = exp(dthalf*const_D[k]);
                vel[i][k] = expd * ( const_C[k] + const_D[k] * vel[i][k] - const_C[k] / expd ) / const_D[k];
            }
            else
            {
                vel[i][k] += ( const_C[k] + const_D[k] * vel[i][k] ) * dthalf +
                    0.5 * (const_D[k] * const_D[k] * vel[i][k] + const_C[k] * const_D[k] ) * dthalf * dthalf;
            }
        }
    }
}

void MSST::propagate_voldot()
{
    const int sd = mdp.msst_direction;
    const double dthalf = 0.5 * mdp.md_dt;
    double p_current = stress(sd, sd);
    double p_msst = mdp.msst_vel * mdp.msst_vel * totmass * (v0 - ucell.omega) / (v0 * v0);
    double const_A = totmass * (p_current - p0 - p_msst) / mdp.msst_qmass;
    double const_B = totmass * mdp.msst_vis / (mdp.msst_qmass * ucell.omega);

    // prevent the increase of volume
    if(ucell.omega > v0 && const_A > 0)
    {
        const_A = -const_A;
    }

    // avoid singularity at B = 0 with Taylor expansion
    double fac = const_B * dthalf;
    if(fac > 1e-6)
    {
        omega[sd] = (omega[sd] + const_A * (exp(fac) - 1) / const_B) * exp(-fac);
    }
    else
    {
        omega[sd] += (const_A - const_B * omega[sd]) * dthalf +
            0.5 * (const_B * const_B * omega[sd] - const_A * const_B) * dthalf * dthalf;
    }
}
