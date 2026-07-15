#include "md_base.h"
#include "md_func.h"
#ifdef __MPI
#include "mpi.h"
#endif
#include "source_io/module_output/print_info.h"
#include <algorithm>

void MD_base::refresh_runtime_storage_from_mdcell()
{
    delete[] allmass;
    delete[] pos;
    delete[] vel;
    delete[] ionmbl;
    delete[] force;

    allmass = new double[mdcell.nlocal()];
    pos = new ModuleBase::Vector3<double>[mdcell.nlocal()];
    vel = new ModuleBase::Vector3<double>[mdcell.nlocal()];
    ionmbl = new ModuleBase::Vector3<int>[mdcell.nlocal()];
    force = new ModuleBase::Vector3<double>[mdcell.nlocal()];
    state_ = MdStateView::from_mdcell(mdcell);
    for (int i = 0; i < mdcell.nlocal(); ++i)
    {
        allmass[i] = state_.mass(i);
        vel[i] = state_.vel(i);
        ionmbl[i] = state_.mbl(i);
        force[i] = state_.force(i);
        pos[i].set(0.0, 0.0, 0.0);
    }
}

void MD_base::sync_velocity_buffer_to_state()
{
    for (int i = 0; i < state_.size(); ++i)
    {
        state_.vel(i) = vel[i];
    }
}

MD_base::MD_base(const Parameter& param_in, MdCell& mdcell_in)
: mdp(param_in.mdp), mdcell(mdcell_in)
{
    my_rank = param_in.globalv.myrank;
    cal_stress = param_in.inp.cal_stress;
    if (mdp.md_seed >= 0)
    {
        srand(mdp.md_seed);
    }

    stop = false;

    assert(mdcell.nlocal() > 0);

    allmass = nullptr;
    pos = nullptr;
    vel = nullptr;
    ionmbl = nullptr;
    force = nullptr;
    virial.create(3, 3);
    stress.create(3, 3);

    assert(ModuleBase::AU_to_FS!=0.0);
    assert(ModuleBase::Hartree_to_K!=0.0);

    /// convert to a.u. unit
    md_dt = mdp.md_dt / ModuleBase::AU_to_FS;
    md_tfirst = mdp.md_tfirst / ModuleBase::Hartree_to_K;
    md_tlast = mdp.md_tlast / ModuleBase::Hartree_to_K;

    step_ = 0;
    step_rst_ = 0;

    refresh_runtime_storage_from_mdcell();
    MD_func::init_vel(mdcell, my_rank, mdp.md_restart, md_tfirst, allmass, frozen_freedom_, ionmbl, vel);
    for (int i = 0; i < mdcell.nlocal(); ++i)
    {
        state_.vel(i) = vel[i];
        state_.mbl(i) = ionmbl[i];
    }
    t_current = MD_func::current_temp(kinetic, mdcell, frozen_freedom_, allmass, vel);
}


MD_base::~MD_base()
{
    delete[] allmass;
    delete[] pos;
    delete[] vel;
    delete[] ionmbl;
    delete[] force;
}


void MD_base::setup(ModuleESolver::ESolver* p_esolver, const std::string& global_readin_dir)
{
    if (mdp.md_restart)
    {
        restart(global_readin_dir);
        refresh_runtime_storage_from_mdcell();
    }

    // mohan add 2026-01-04
    const int stress_step = 0;
    const int force_step = 0;
    const int istep_print = step_ + step_rst_ + 1;

	ModuleIO::print_screen(stress_step, force_step, istep_print);

    MD_func::force_virial(p_esolver, step_, mdcell, potential, force, cal_stress, virial);
    MD_func::compute_stress(mdcell, vel, allmass, cal_stress, virial, stress);

    return;
}


void MD_base::first_half(std::ofstream& ofs)
{
    update_vel(force);
    update_pos();

    return;
}


void MD_base::second_half()
{
    update_vel(force);

    return;
}


void MD_base::update_pos()
{
    for (int i = 0; i < state_.size(); ++i)
    {
        for (int k = 0; k < 3; ++k)
        {
            if (state_.mbl(i)[k])
            {
                pos[i][k] = state_.vel(i)[k] * md_dt / mdcell.lat0();
            }
            else
            {
                pos[i][k] = 0;
            }
        }
        pos[i] = pos[i] * mdcell.GT();
        state_.frac(i) += pos[i];
        state_.frac(i).x -= std::floor(state_.frac(i).x);
        state_.frac(i).y -= std::floor(state_.frac(i).y);
        state_.frac(i).z -= std::floor(state_.frac(i).z);
        state_.cart(i) = state_.frac(i) * mdcell.latvec();
    }

#ifdef __MPI
    mdcell.migrate_owned_atoms();
    refresh_runtime_storage_from_mdcell();
#endif

    return;
}


void MD_base::update_vel(const ModuleBase::Vector3<double>* force)
{
    for (int i = 0; i < state_.size(); ++i)
    {
        for (int k = 0; k < 3; ++k)
        {
            if (state_.mbl(i)[k])
            {
                state_.vel(i)[k] += 0.5 * force[i][k] * md_dt / state_.mass(i);
            }
        }
    }
    for (int i = 0; i < state_.size(); ++i)
    {
        vel[i] = state_.vel(i);
    }
    return;
}


void MD_base::print_md(std::ofstream& ofs, const bool& cal_stress)
{
    t_current = MD_func::current_temp(kinetic, mdcell, frozen_freedom_, allmass, vel);

    if (my_rank!=0)
    {
        return;
    }

    assert(ModuleBase::BOHR_RADIUS_SI>0.0);

    const double unit_transform = ModuleBase::HARTREE_SI / pow(ModuleBase::BOHR_RADIUS_SI, 3) * 1.0e-8;
    double press = 0.0;
    for (int i = 0; i < 3; i++)
    {
        press += stress(i, i) / 3;
    }

    // screen output
    std::cout << std::setprecision(8);
    std::cout << " ------------------------------------------------------------------------------------------------"
              << std::endl;
    std::cout << " " << std::left << std::setw(20) << "Energy (Ry)" << std::left << std::setw(20) << "Potential (Ry)"
              << std::left << std::setw(20) << "Kinetic (Ry)" << std::left << std::setw(20) << "Temperature (K)";

    if (cal_stress)
    {
        std::cout << std::left << std::setw(20) << "Pressure (kbar)";
    }

    std::cout << std::endl;
    std::cout << " " << std::left << std::setw(20) << 2 * (potential + kinetic) << std::left << std::setw(20)
              << 2 * potential << std::left << std::setw(20) << 2 * kinetic << std::left << std::setw(20)
              << t_current * ModuleBase::Hartree_to_K;

    if (cal_stress)
    {
        std::cout << std::left << std::setw(20) << press * unit_transform;
    }

    std::cout << std::endl;
    std::cout << " ------------------------------------------------------------------------------------------------"
              << std::endl;

    // running_log output
    ofs.unsetf(std::ios::fixed);
    ofs << std::setprecision(8);

    if (cal_stress)
    {
        MD_func::print_stress(ofs, virial, stress);
	    ofs << std::endl;
    }

    ofs << " ------------------------------------------------------------------------------------------------"
        << std::endl;
    ofs << " " << std::left << std::setw(20) << "Energy (Ry)" << std::left << std::setw(20) << "Potential (Ry)"
        << std::left << std::setw(20) << "Kinetic (Ry)" << std::left << std::setw(20) << "Temperature (K)";

    if (cal_stress)
    {
        ofs << std::left << std::setw(20) << "Pressure (kbar)";
    }

    ofs << std::endl;
    ofs << " " << std::left << std::setw(20) << 2 * (potential + kinetic) << std::left << std::setw(20) << 2 * potential
        << std::left << std::setw(20) << 2 * kinetic << std::left << std::setw(20)
        << t_current * ModuleBase::Hartree_to_K;

    if (cal_stress)
    {
        ofs << std::left << std::setw(20) << press * unit_transform;
    }

    ofs << std::endl;
    ofs << " ------------------------------------------------------------------------------------------------"
        << std::endl;
    ofs << std::endl;
    return;
}


void MD_base::write_restart(const std::string& global_out_dir)
{
    if (!my_rank)
    {
        std::stringstream ssc;
        ssc << global_out_dir << "Restart_md.txt";
        std::ofstream file(ssc.str().c_str());

        file << step_ + step_rst_ << std::endl;
        file << md_tfirst << std::endl;
        file.close();
    }
#ifdef __MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    return;
}


void MD_base::restart(const std::string& global_readin_dir)
{
    MD_func::current_md_info(my_rank, global_readin_dir, step_rst_, md_tfirst);

    return;
}
