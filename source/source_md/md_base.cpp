#include "md_base.h"
#include "md_func.h"
#include "source_cell/unitcell.h"
#ifdef __MPI
#include "mpi.h"
#endif
#include "source_io/module_output/print_info.h"
#include <algorithm>
#include <iomanip>

MD_base::MD_base(const Parameter& param_in, MDCell& mdcell_in)
: mdp(param_in.mdp), mdcell(mdcell_in)
{
#ifdef __MPI
    my_rank = mdcell.mpi_rank();
#else
    my_rank = param_in.globalv.myrank;
#endif
    cal_stress = param_in.inp.cal_stress;
    if (mdp.md_seed >= 0)
    {
        srand(mdp.md_seed + my_rank * 104729);
    }

    stop = false;

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

    MD_func::init_vel(mdcell, mdp.md_restart, md_tfirst, frozen_freedom_);
    t_current = MD_func::current_temp(kinetic, mdcell, frozen_freedom_);
}


MD_base::~MD_base() {}


void MD_base::setup(ModuleESolver::ESolver* p_esolver, const std::string& global_readin_dir)
{
    if (mdp.md_restart)
    {
        restart(global_readin_dir);
    }

    // mohan add 2026-01-04
    const int stress_step = 0;
    const int force_step = 0;
    const int istep_print = step_ + step_rst_ + 1;

	ModuleIO::print_screen(stress_step, force_step, istep_print);

    MD_func::force_virial(p_esolver, step_, mdcell, potential, cal_stress, virial);
    MD_func::compute_stress(mdcell, cal_stress, virial, stress);
    if (mdcell.has_backing_unitcell())
    {
        mdcell.backing_unitcell().ionic_position_updated = true;
    }

    return;
}


void MD_base::first_half(std::ofstream& ofs)
{
    update_vel();
    update_pos();

    return;
}


void MD_base::second_half()
{
    update_vel();

    return;
}


void MD_base::update_pos()
{
    std::vector<LocalAtom>& atoms = mdcell.mutable_owned_atoms();
    for (std::size_t i = 0; i < atoms.size(); ++i)
    {
        LocalAtom& atom = atoms[i];
        ModuleBase::Vector3<double> pos;
        for (int k = 0; k < 3; ++k)
        {
            if (atom.mbl[k])
            {
                pos[k] = atom.vel[k] * md_dt / mdcell.lat0();
            }
            else
            {
                pos[k] = 0;
            }
        }
        pos = pos * mdcell.GT();
        atom.frac += pos;
        atom.frac.x -= std::floor(atom.frac.x);
        atom.frac.y -= std::floor(atom.frac.y);
        atom.frac.z -= std::floor(atom.frac.z);
        atom.cart = atom.frac * mdcell.latvec();
    }

    mdcell.migrate_owned_atoms();

    return;
}


void MD_base::update_vel()
{
    std::vector<LocalAtom>& atoms = mdcell.mutable_owned_atoms();
    for (std::size_t i = 0; i < atoms.size(); ++i)
    {
        LocalAtom& atom = atoms[i];
        for (int k = 0; k < 3; ++k)
        {
            if (atom.mbl[k])
            {
                atom.vel[k] += 0.5 * atom.force[k] * md_dt / atom.mass;
            }
        }
    }
    return;
}


void MD_base::print_md(std::ofstream& ofs, const bool& cal_stress)
{
    t_current = MD_func::current_temp(kinetic, mdcell, frozen_freedom_);

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
    MPI_Barrier(mdcell.communicator());
#endif

    return;
}


void MD_base::restart(const std::string& global_readin_dir)
{
    MD_func::current_md_info(mdcell, global_readin_dir, step_rst_, md_tfirst);

    return;
}
