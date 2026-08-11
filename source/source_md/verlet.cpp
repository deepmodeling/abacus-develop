#include "verlet.h"

#include "md_func.h"
#include "source_base/timer.h"

#ifdef __MPI
#include <mpi.h>
#endif

Verlet::Verlet(const Parameter& param_in, MDCell& mdcell_in) : MD_base(param_in, mdcell_in)
{
}

Verlet::~Verlet()
{
}


void Verlet::setup(ModuleESolver::ESolver* p_esolver, const std::string& global_readin_dir)
{
    ModuleBase::TITLE("Verlet", "setup");
    ModuleBase::timer::start("Verlet", "setup");

    MD_base::setup(p_esolver, global_readin_dir);

    ModuleBase::timer::end("Verlet", "setup");
}


void Verlet::first_half(std::ofstream& ofs)
{
    ModuleBase::TITLE("Verlet", "first_half");
    ModuleBase::timer::start("Verlet", "first_half");

    MD_base::update_vel();
    MD_base::update_pos();

    ModuleBase::timer::end("Verlet", "first_half");
}


void Verlet::second_half()
{
    ModuleBase::TITLE("Verlet", "second_half");
    ModuleBase::timer::start("Verlet", "second_half");

    MD_base::update_vel();
    apply_thermostat();

    ModuleBase::timer::end("Verlet", "second_half");
}


void Verlet::apply_thermostat(void)
{
    double t_target = 0.0;
    t_current = MD_func::current_temp(kinetic, mdcell, frozen_freedom_);

    if (mdp.md_type == "nve")
    {
    }
    else if (mdp.md_thermostat == "rescaling")
    {
        t_target = MD_func::target_temp(step_ + step_rst_, mdp.md_nstep, md_tfirst, md_tlast);
        if (std::abs(t_target - t_current) * ModuleBase::Hartree_to_K > mdp.md_tolerance)
        {
            thermalize(0, t_current, t_target);
        }
    }
    else if (mdp.md_thermostat == "rescale_v")
    {
        if ((step_ + step_rst_) % mdp.md_nraise == 0)
        {
            t_target = MD_func::target_temp(step_ + step_rst_, mdp.md_nstep, md_tfirst, md_tlast);
            thermalize(0, t_current, t_target);
        }
    }
    else if (mdp.md_thermostat == "anderson")
    {
        for (LocalAtom& atom : mdcell.mutable_owned_atoms())
        {
            if (static_cast<double>(std::rand()) / RAND_MAX <= 1.0 / mdp.md_nraise)
            {
                const double deviation = sqrt(md_tlast / atom.mass);
                for (int k = 0; k < 3; ++k)
                {
                    if (atom.mbl[k])
                    {
                        atom.vel[k] = deviation * MD_func::gaussrand();
                    }
                }
            }
        }
    }
    else if (mdp.md_thermostat == "berendsen")
    {
        t_target = MD_func::target_temp(step_ + step_rst_, mdp.md_nstep, md_tfirst, md_tlast);
        thermalize(mdp.md_nraise, t_current, t_target);
    }
    else if (mdp.md_thermostat == "csvr")
    {
        t_target = MD_func::target_temp(step_ + step_rst_, mdp.md_nstep, md_tfirst, md_tlast);
        apply_csvr(t_current, t_target);
    }
    else
    {
        ModuleBase::WARNING_QUIT("Verlet", "No such thermostat!");
    }
}


void Verlet::thermalize(const int& nraise, const double& current_temp, const double& target_temp)
{
    double fac = 0.0;
    if (nraise > 0 && current_temp > 0 && target_temp > 0)
    {
        fac = sqrt(1 + (target_temp / current_temp - 1) / nraise);
    }
    else if (nraise == 0 && current_temp > 0 && target_temp > 0)
    {
        fac = sqrt(target_temp / current_temp);
    }

    for (LocalAtom& atom : mdcell.mutable_owned_atoms()) atom.vel *= fac;
}


void Verlet::apply_csvr(const double& current_temp, const double& target_temp)
{
    // CSVR thermostat: Canonical Sampling through Velocity Rescaling
    // Reference: G. Bussi, D. Donadio, M. Parrinello, J. Chem. Phys. 126, 014101 (2007)

    if (current_temp <= 0.0 || target_temp <= 0.0)
    {
        return;
    }

    // Get degrees of freedom (3N - frozen)
    std::int64_t ndeg = MD_func::global_dof(mdcell, frozen_freedom_);

    // Calculate kinetic energies
    double kin_energy = current_temp * static_cast<double>(ndeg) * 0.5;  // in Hartree
    double kin_target = target_temp * ndeg * 0.5;   // in Hartree

    // Calculate tau parameter (characteristic time scale / dt)
    double taut = mdp.md_csvr_tau / mdp.md_dt;

    // Calculate decay factor
    double factor = 0.0;
    if (taut > 0.1)
    {
        factor = exp(-1.0 / taut);
    }

    double scale = 1.0;
#ifdef __MPI
    if (mdcell.mpi_size() > 1)
    {
        if (mdcell.mpi_rank() == 0)
        {
            const double rr = MD_func::gaussrand();
            double sumnoises = 0.0;
            for (int i = 0; i < ndeg - 1; ++i)
            {
                const double random_value = MD_func::gaussrand();
                sumnoises += random_value * random_value;
            }
            const double factor2 = (1.0 - factor) * kin_target / kin_energy / ndeg;
            const double resample = std::max(0.0,
                                             factor + factor2 * (rr * rr + sumnoises)
                                                 + 2.0 * rr * sqrt(factor * factor2));
            scale = sqrt(resample);
        }
        MPI_Bcast(&scale, 1, MPI_DOUBLE, 0, mdcell.communicator());
    }
    else
#endif
    {
        const double rr = MD_func::gaussrand();
        double sumnoises = 0.0;
        for (int i = 0; i < ndeg - 1; ++i)
        {
            const double random_value = MD_func::gaussrand();
            sumnoises += random_value * random_value;
        }
        const double factor2 = (1.0 - factor) * kin_target / kin_energy / ndeg;
        const double resample = std::max(0.0,
                                         factor + factor2 * (rr * rr + sumnoises)
                                             + 2.0 * rr * sqrt(factor * factor2));
        scale = sqrt(resample);
    }

    // Apply velocity scaling
    for (LocalAtom& atom : mdcell.mutable_owned_atoms()) atom.vel *= scale;
}


void Verlet::print_md(std::ofstream& ofs, const bool& cal_stress)
{
    MD_base::print_md(ofs, cal_stress);
    return;
}


void Verlet::write_restart(const std::string& global_out_dir)
{
    MD_base::write_restart(global_out_dir);
    return;
}


void Verlet::restart(const std::string& global_readin_dir)
{
    MD_base::restart(global_readin_dir);
    return;
}
