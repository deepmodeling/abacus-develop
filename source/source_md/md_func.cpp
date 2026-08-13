#include "md_func.h"

#include "source_base/global_variable.h"
#include "source_base/timer.h"
#include "source_io/module_parameter/parameter.h"

#include <cerrno>
#include <cstring>
#include <fcntl.h>
#include <unistd.h>


namespace MD_func
{
#ifdef __MPI
namespace
{
bool write_dump_at(const int file, const std::string& data, const MPI_Offset offset)
{
    std::size_t written = 0;
    while (written < data.size())
    {
        const ssize_t count = pwrite(file, data.data() + written, data.size() - written,
                                     static_cast<off_t>(offset + written));
        if (count <= 0) return false;
        written += static_cast<std::size_t>(count);
    }
    return true;
}
}
#endif
double gaussrand()
{
    static double v1=0.0;
    static double v2=0.0;
    static double S=0.0;
    static int phase = 0;
    double xx=0.0;

    if (phase == 0)
    {
        do
        {
            double U1 = static_cast<double>(std::rand()) / RAND_MAX;
            double U2 = static_cast<double>(std::rand()) / RAND_MAX;

            v1 = 2.0 * U1 - 1.0;
            v2 = 2.0 * U2 - 1.0;
            S = v1 * v1 + v2 * v2;
        } while (S >= 1 || S == 0);

        xx = v1 * sqrt(-2.0 * log(S) / S);
    }
    else
    {
        xx = v2 * sqrt(-2.0 * log(S) / S);
    }

    phase = 1 - phase;

    return xx;
}

double kinetic_energy(const int& natom, const ModuleBase::Vector3<double>* vel, const double* allmass)
{
    double ke = 0;

#pragma omp parallel for reduction(+:ke) schedule(static) if (natom >= 256)
    for (int ion = 0; ion < natom; ++ion)
    {
        ke += 0.5 * allmass[ion] * vel[ion].norm2();
    }

    return ke;
}

void compute_stress(const UnitCell& unit_in,
                    const ModuleBase::Vector3<double>* vel,
                    const double* allmass,
                    const bool& cal_stress,
                    const ModuleBase::matrix& virial,
                    ModuleBase::matrix& stress)
{
    if (cal_stress)
    {
        ModuleBase::matrix temperature_tensor(3, 3);
        temp_vector(unit_in.nat, vel, allmass, temperature_tensor);
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                stress(i, j) = virial(i, j) + temperature_tensor(i, j) / unit_in.omega;
            }
        }
    }

    return;
}

void compute_stress(const MDCell& mdcell,
                    const bool& cal_stress,
                    const ModuleBase::matrix& virial,
                    ModuleBase::matrix& stress)
{
    if (!cal_stress) return;
    ModuleBase::matrix t_vector(3, 3);
    for (std::size_t i = 0; i < mdcell.owned_atoms().size(); ++i)
    {
        const LocalAtom& atom = mdcell.owned_atoms()[i];
        for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b)
            t_vector(a, b) += atom.mass * atom.vel[a] * atom.vel[b];
    }
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, t_vector.c, 9, MPI_DOUBLE, MPI_SUM, mdcell.communicator());
#endif
    for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j)
        stress(i, j) = virial(i, j) + t_vector(i, j) / mdcell.omega();
}

void read_vel(const UnitCell& unit_in, ModuleBase::Vector3<double>* vel)
{
    int iat = 0;
    for (int it = 0; it < unit_in.ntype; ++it)
    {
        for (int ia = 0; ia < unit_in.atoms[it].na; ++ia)
        {
            vel[iat] = unit_in.atoms[it].vel[ia];
            if (unit_in.atoms[it].mbl[ia].x == 0)
            {
                vel[iat].x = 0;
            }
            if (unit_in.atoms[it].mbl[ia].y == 0)
            {
                vel[iat].y = 0;
            }
            if (unit_in.atoms[it].mbl[ia].z == 0)
            {
                vel[iat].z = 0;
            }
            ++iat;
        }
    }
    assert(iat == unit_in.nat);

    return;
}

void rescale_vel(const int& natom,
                 const double& temperature,
                 const double* allmass,
                 const int& frozen_freedom,
                 ModuleBase::Vector3<double>* vel)
{
    double factor = 0.0;
    if (3 * natom == frozen_freedom || temperature == 0)
    {
        factor = 0;
    }
    else
    {
        factor = 0.5 * (3 * natom - frozen_freedom) * temperature / kinetic_energy(natom, vel, allmass);
    }

#pragma omp parallel for schedule(static) if (natom >= 256)
    for (int i = 0; i < natom; i++)
    {
        vel[i] = vel[i] * sqrt(factor);
    }
}

void rand_vel(const int& natom,
              const double& temperature,
              const double* allmass,
              const int& frozen_freedom,
              const ModuleBase::Vector3<int> frozen,
              const ModuleBase::Vector3<int>* ionmbl,
              const int& my_rank,
              ModuleBase::Vector3<double>* vel)
{
    if (!my_rank)
    {
        double tot_mass = 0;
        ModuleBase::Vector3<double> tot_momentum;
        for (int i = 0; i < natom; i++)
        {
            tot_mass += allmass[i];
            double sigma = sqrt(temperature / allmass[i]);
            for (int k = 0; k < 3; ++k)
            {
                if (ionmbl[i][k] == 0)
                {
                    vel[i][k] = 0;
                }
                else
                {
                    vel[i][k] = gaussrand() * sigma;
                }

                if (frozen[k] == 0)
                {
                    tot_momentum[k] += allmass[i] * vel[i][k];
                }
            }
        }

        for (int k = 0; k < 3; ++k)
        {
            if (frozen[k] == 0)
            {
                for (int i = 0; i < natom; i++)
                {
                    vel[i][k] -= tot_momentum[k] / tot_mass;
                }
            }
        }

        // rescale the velocity to the target temperature
        rescale_vel(natom, temperature, allmass, frozen_freedom, vel);
    }

#ifdef __MPI
    MPI_Bcast(vel, natom * 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

    return;
}

void init_vel(const UnitCell& unit_in,
              const int& my_rank,
              const bool& restart,
              double& temperature,
              double* allmass,
              int& frozen_freedom,
              ModuleBase::Vector3<int>* ionmbl,
              ModuleBase::Vector3<double>* vel)
{
    ModuleBase::Vector3<int> frozen;
    get_mass_mbl(unit_in, allmass, frozen, ionmbl);
    frozen_freedom = frozen.x + frozen.y + frozen.z;
    if (frozen.x == 0)
    {
        ++frozen_freedom;
    }
    if (frozen.y == 0)
    {
        ++frozen_freedom;
    }
    if (frozen.z == 0)
    {
        ++frozen_freedom;
    }

    if (unit_in.init_vel)
    {
        std::cout << " Reading velocities from STRU file" << std::endl;
        read_vel(unit_in, vel);
        double kinetic = 0.0;
        double t_current = MD_func::current_temp(kinetic, unit_in.nat, frozen_freedom, allmass, vel);
        if (restart)
        {
            std::cout << " Restart MD, current temperature is " << t_current * ModuleBase::Hartree_to_K << " K"
                      << std::endl;
        }
        else if (temperature < 0)
        {
            std::cout << " Autoset the initial tempearture to " << t_current * ModuleBase::Hartree_to_K << " K"
                      << std::endl;
            temperature = t_current;
        }
        else
        {
            std::cout << " Initial temeprature from INPUT is " << temperature * ModuleBase::Hartree_to_K << " K"
                      << std::endl;
            std::cout << " Reading temperature from STRU is " << t_current * ModuleBase::Hartree_to_K << " K"
                      << std::endl;
            std::cout << " Rescale velocties to initial temperature" << std::endl;
            rescale_vel(unit_in.nat, temperature, allmass, frozen_freedom, vel);
        }
    }
    else
    {
        std::cout << " Random velocities according to initial temperature " 
                  << temperature * ModuleBase::Hartree_to_K << " K"
                  << std::endl;
        rand_vel(unit_in.nat, temperature, allmass, frozen_freedom, frozen, ionmbl, my_rank, vel);
    }
}

void init_vel(MDCell& mdcell,
              const bool& restart,
              double& temperature,
              std::int64_t& frozen_freedom)
{
    std::vector<LocalAtom>& atoms = mdcell.mutable_owned_atoms();
    ModuleBase::Vector3<std::int64_t> frozen(0, 0, 0);
    for (std::size_t i = 0; i < atoms.size(); ++i)
        for (int k = 0; k < 3; ++k) if (!atoms[i].mbl[k]) ++frozen[k];
#ifdef __MPI
    if (mdcell.mpi_size() > 1)
    {
        MPI_Allreduce(MPI_IN_PLACE, &frozen.x, 3, MPI_INT64_T, MPI_SUM, mdcell.communicator());
    }
#endif
    frozen_freedom = frozen.x + frozen.y + frozen.z;
    if (!frozen.x) ++frozen_freedom;
    if (!frozen.y) ++frozen_freedom;
    if (!frozen.z) ++frozen_freedom;
    if (mdcell.init_vel())
    {
        double kinetic = 0.0;
        const double current = current_temp(kinetic, mdcell, frozen_freedom);
        if (!restart && current > 0.0 && temperature > 0.0)
        {
            const double factor = sqrt(temperature / current);
            for (LocalAtom& atom : atoms) atom.vel *= factor;
        }
        return;
    }

    double local_mass = 0.0;
    ModuleBase::Vector3<double> momentum(0.0, 0.0, 0.0);
    for (LocalAtom& atom : atoms)
    {
        local_mass += atom.mass;
        for (int k = 0; k < 3; ++k)
        {
            atom.vel[k] = atom.mbl[k] ? gaussrand() * sqrt(temperature / atom.mass) : 0.0;
            if (frozen[k] == 0) momentum[k] += atom.mass * atom.vel[k];
        }
    }
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, &local_mass, 1, MPI_DOUBLE, MPI_SUM, mdcell.communicator());
    MPI_Allreduce(MPI_IN_PLACE, &momentum.x, 3, MPI_DOUBLE, MPI_SUM, mdcell.communicator());
#endif
    for (int k = 0; k < 3; ++k)
    {
        if (frozen[k] == 0 && local_mass > 0.0)
        {
            for (LocalAtom& atom : atoms) atom.vel[k] -= momentum[k] / local_mass;
        }
    }
    double kinetic = 0.0;
    const double current = current_temp(kinetic, mdcell, frozen_freedom);
    if (current > 0.0 && temperature > 0.0)
    {
        const double factor = sqrt(temperature / current);
        for (LocalAtom& atom : atoms) atom.vel *= factor;
    }
    static_cast<void>(restart);
}

void force_virial(ModuleESolver::ESolver* p_esolver,
                  const int& istep,
                  UnitCell& unit_in,
                  double& potential,
                  ModuleBase::Vector3<double>* force,
                  const bool& cal_stress,
                  ModuleBase::matrix& virial)
{
    ModuleBase::TITLE("MD_func", "force_virial");
    ModuleBase::timer::start("MD_func", "force_virial");

    p_esolver->runner(unit_in, istep);

    potential = p_esolver->cal_energy();

    ModuleBase::matrix force_temp(unit_in.nat, 3);
    p_esolver->cal_force(unit_in, force_temp);

    if (cal_stress)
    {
        p_esolver->cal_stress(unit_in, virial);
    }

    /// convert Rydberg to Hartree
    potential *= 0.5;
    force_temp *= 0.5;
    virial *= 0.5;

    const int natom = unit_in.nat;
#pragma omp parallel for schedule(static) if (natom >= 256)
    for (int i = 0; i < natom; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            force[i][j] = force_temp(i, j);
        }
    }

    ModuleBase::timer::end("MD_func", "force_virial");

    return;
}

void force_virial(ModuleESolver::ESolver* p_esolver,
                  const int& istep,
                  MDCell& mdcell,
                  double& potential,
                  const bool& cal_stress,
                  ModuleBase::matrix& virial)
{
    ModuleBase::TITLE("MD_func", "force_virial");
    ModuleBase::timer::start("MD_func", "force_virial");
    if (p_esolver->supports_mdcell())
    {
        p_esolver->runner(static_cast<BaseCell&>(mdcell), istep);
        potential = 0.5 * p_esolver->cal_energy();
        ModuleBase::matrix local_force;
        p_esolver->cal_force(static_cast<BaseCell&>(mdcell), local_force);
        for (LocalAtom& atom : mdcell.mutable_owned_atoms()) atom.force *= 0.5;
        if (cal_stress) { p_esolver->cal_stress(static_cast<BaseCell&>(mdcell), virial); virial *= 0.5; }
    }
    else
    {
        if (!mdcell.has_backing_unitcell()) ModuleBase::WARNING_QUIT("MD_func::force_virial", "This ESolver requires UnitCell, but MDCell has no backing UnitCell.");
        UnitCell& ucell = mdcell.backing_unitcell();
        std::vector<std::vector<ModuleBase::Vector3<double>>> backing_velocities(
            static_cast<std::size_t>(ucell.ntype));
        for (int it = 0; it < ucell.ntype; ++it)
        {
            backing_velocities[static_cast<std::size_t>(it)] = ucell.atoms[it].vel;
        }
        mdcell.sync_backing_unitcell();
        for (int it = 0; it < ucell.ntype; ++it)
        {
            std::copy(backing_velocities[static_cast<std::size_t>(it)].begin(),
                      backing_velocities[static_cast<std::size_t>(it)].end(),
                      ucell.atoms[it].vel.begin());
        }
        p_esolver->runner(ucell, istep); potential = 0.5 * p_esolver->cal_energy();
        ModuleBase::matrix full_force(ucell.nat, 3); p_esolver->cal_force(ucell, full_force); full_force *= 0.5;
        if (cal_stress) { p_esolver->cal_stress(ucell, virial); virial *= 0.5; }
        std::vector<int> offsets(ucell.ntype + 1, 0); for (int it=0; it<ucell.ntype; ++it) offsets[it+1]=offsets[it]+ucell.atoms[it].na;
        for (LocalAtom& atom : mdcell.mutable_owned_atoms()) { const int iat=offsets[atom.type]+atom.type_index; atom.force.set(full_force(iat,0),full_force(iat,1),full_force(iat,2)); }
    }
    ModuleBase::timer::end("MD_func", "force_virial");
}


void print_stress(std::ofstream& ofs, const ModuleBase::matrix& virial, const ModuleBase::matrix& stress)
{
    double stress_scalar = 0.0;
    double virial_scalar = 0.0;

    for (int i = 0; i < 3; i++)
    {
        stress_scalar += stress(i, i) / 3.0;
        virial_scalar += virial(i, i) / 3.0;
    }

    const double unit_transform = ModuleBase::HARTREE_SI / pow(ModuleBase::BOHR_RADIUS_SI, 3) * 1.0e-8;

    
    ofs << " ELECTRONIC      PART OF STRESS: " << virial_scalar * unit_transform << " kbar" << std::endl;
    ofs << " IONIC (KINETIC) PART OF STRESS: " << (stress_scalar - virial_scalar) * unit_transform << " kbar" << std::endl;
    ofs << " MD PRESSURE (ELECTRONS+IONS)  : " << stress_scalar * unit_transform << " kbar" << std::endl;

    // one should use 'print_stress' function in ../source/source_io/output_log.cpp
/*
    ofs.unsetf(std::ios::fixed);
    ofs << std::setprecision(8) << std::endl;
    ModuleBase::GlobalFunc::NEW_PART("MD STRESS (kbar)");
    for (int i = 0; i < 3; i++)
    {
        ofs << std::setw(15) << stress(i, 0) * unit_transform << std::setw(15) << stress(i, 1) * unit_transform
            << std::setw(15) << stress(i, 2) * unit_transform << std::endl;
    }
    ofs << std::setiosflags(std::ios::left);
*/

    return;
}

void dump_info(const int& step,
               const std::string& global_out_dir,
               const UnitCell& unit_in,
               const Parameter& param_in,
               const ModuleBase::matrix& virial,
               const ModuleBase::Vector3<double>* force,
               const ModuleBase::Vector3<double>* vel)
{
    if (param_in.globalv.myrank)
    {
        return;
    }

    std::stringstream file;
    file << global_out_dir << "MD_dump";
    std::ofstream ofs;
    if (step == 0)
    {
        ofs.open(file.str(), std::ios::trunc);
    }
    else
    {
        ofs.open(file.str(), std::ios::app);
    }

    const double unit_pos = unit_in.lat0 / ModuleBase::ANGSTROM_AU;                                  ///< Angstrom
    const double unit_vel = 1.0 / ModuleBase::ANGSTROM_AU / ModuleBase::AU_to_FS;                    ///< Angstrom/fs
    const double unit_virial = ModuleBase::HARTREE_SI / pow(ModuleBase::BOHR_RADIUS_SI, 3) * 1.0e-8; ///< kBar
    const double unit_force = ModuleBase::Hartree_to_eV * ModuleBase::ANGSTROM_AU;                   ///< eV/Angstrom

    ofs << "MDSTEP:  " << step << std::endl;
    ofs << std::setprecision(12) << std::setiosflags(std::ios::fixed);

    ofs << "LATTICE_CONSTANT: " << unit_in.lat0_angstrom << " Angstrom" << std::endl;

    ofs << "LATTICE_VECTORS" << std::endl;
    ofs << "  " << unit_in.latvec.e11 << "  " << unit_in.latvec.e12 << "  " << unit_in.latvec.e13 << std::endl;
    ofs << "  " << unit_in.latvec.e21 << "  " << unit_in.latvec.e22 << "  " << unit_in.latvec.e23 << std::endl;
    ofs << "  " << unit_in.latvec.e31 << "  " << unit_in.latvec.e32 << "  " << unit_in.latvec.e33 << std::endl;

    if (param_in.inp.cal_stress && param_in.mdp.dump_virial)
    {
        ofs << "VIRIAL (kbar)" << std::endl;
        for (int i = 0; i < 3; ++i)
        {
            ofs << "  " << virial(i, 0) * unit_virial << "  " << virial(i, 1) * unit_virial << "  "
                << virial(i, 2) * unit_virial << std::endl;
        }
    }

    ofs << "INDEX    LABEL    POSITION (Angstrom)";
    if (param_in.mdp.dump_force)
    {
        ofs << "    FORCE (eV/Angstrom)";
    }
    if (param_in.mdp.dump_vel)
    {
        ofs << "    VELOCITY (Angstrom/fs)";
    }
    ofs << std::endl;

    int index = 0;
    for (int it = 0; it < unit_in.ntype; ++it)
    {
        for (int ia = 0; ia < unit_in.atoms[it].na; ++ia)
        {
            ofs << "  " << index << "  " << unit_in.atoms[it].label << "  " << unit_in.atoms[it].tau[ia].x * unit_pos
                << "  " << unit_in.atoms[it].tau[ia].y * unit_pos << "  " << unit_in.atoms[it].tau[ia].z * unit_pos;

            if (param_in.mdp.dump_force)
            {
                ofs << "  " << force[index].x * unit_force << "  " << force[index].y * unit_force << "  "
                    << force[index].z * unit_force;
            }

            if (param_in.mdp.dump_vel)
            {
                ofs << "  " << vel[index].x * unit_vel << "  " << vel[index].y * unit_vel << "  "
                    << vel[index].z * unit_vel;
            }
            ofs << std::endl;
            index++;
        }
    }

    ofs << std::endl;
    ofs << std::endl;
    ofs.close();

    return;
}

void dump_info(const int& step,
               const std::string& global_out_dir,
               const MDCell& mdcell,
               const Parameter& param_in,
               const ModuleBase::matrix& virial)
{
    std::stringstream file;
    file << global_out_dir << "MD_dump";
    const double unit_pos = mdcell.lat0() / ModuleBase::ANGSTROM_AU;
    const double unit_vel = 1.0 / ModuleBase::ANGSTROM_AU / ModuleBase::AU_to_FS;
    const double unit_virial = ModuleBase::HARTREE_SI / pow(ModuleBase::BOHR_RADIUS_SI, 3) * 1.0e-8;
    const double unit_force = ModuleBase::Hartree_to_eV * ModuleBase::ANGSTROM_AU;
    std::ostringstream header;
    header << "MDSTEP:  " << step << "\n" << std::fixed << std::setprecision(12);
    header << "LATTICE_CONSTANT: " << mdcell.lat0() * ModuleBase::BOHR_TO_A << " Angstrom\nLATTICE_VECTORS\n";
    header << "  " << mdcell.latvec().e11 << "  " << mdcell.latvec().e12 << "  " << mdcell.latvec().e13 << "\n";
    header << "  " << mdcell.latvec().e21 << "  " << mdcell.latvec().e22 << "  " << mdcell.latvec().e23 << "\n";
    header << "  " << mdcell.latvec().e31 << "  " << mdcell.latvec().e32 << "  " << mdcell.latvec().e33 << "\n";
    if (param_in.inp.cal_stress && param_in.mdp.dump_virial)
    {
        header << "VIRIAL (kbar)\n";
        for (int i = 0; i < 3; ++i)
            header << "  " << virial(i, 0) * unit_virial << "  " << virial(i, 1) * unit_virial << "  " << virial(i, 2) * unit_virial << "\n";
    }
    header << "INDEX    LABEL    POSITION (Angstrom)";
    if (param_in.mdp.dump_force) header << "    FORCE (eV/Angstrom)";
    if (param_in.mdp.dump_vel) header << "    VELOCITY (Angstrom/fs)";
    header << "\n";
    std::vector<std::int64_t> type_offsets(mdcell.type_atom_counts().size() + 1, 0);
    for (std::size_t it = 0; it < mdcell.type_atom_counts().size(); ++it)
    {
        type_offsets[it + 1] = type_offsets[it] + mdcell.type_atom_counts()[it];
    }
    std::ostringstream local;
    local << std::fixed << std::setprecision(12);
    for (int i = 0; i < mdcell.nlocal(); ++i)
    {
        const LocalAtom& atom = mdcell.owned_atoms()[static_cast<std::size_t>(i)];
        local << "  " << type_offsets[static_cast<std::size_t>(atom.type)] + atom.type_index
              << "  " << mdcell.type_labels()[static_cast<std::size_t>(atom.type)]
              << "  " << atom.cart.x * unit_pos << "  " << atom.cart.y * unit_pos << "  " << atom.cart.z * unit_pos;
        if (param_in.mdp.dump_force)
            local << "  " << atom.force.x * unit_force << "  " << atom.force.y * unit_force << "  " << atom.force.z * unit_force;
        if (param_in.mdp.dump_vel)
            local << "  " << atom.vel.x * unit_vel << "  " << atom.vel.y * unit_vel << "  " << atom.vel.z * unit_vel;
        local << "\n";
    }
#ifdef __MPI
    const MPI_Comm comm = mdcell.communicator();
    int rank = 0; MPI_Comm_rank(comm, &rank);
    MPI_Offset base = 0;
    if (rank == 0)
    {
        const int file_descriptor = open(file.str().c_str(), O_CREAT | O_WRONLY | (step == 0 ? O_TRUNC : O_APPEND), 0666);
        if (file_descriptor < 0) ModuleBase::WARNING_QUIT("MD_func::dump_info", "cannot open MD_dump.");
        base = lseek(file_descriptor, 0, SEEK_END);
        if (!write_dump_at(file_descriptor, header.str(), base) || close(file_descriptor) != 0) ModuleBase::WARNING_QUIT("MD_func::dump_info", "cannot write MD_dump header.");
    }
    MPI_Bcast(&base, 1, MPI_OFFSET, 0, comm);
    const MPI_Offset atom_offset = base + static_cast<MPI_Offset>(header.str().size());
    const MPI_Offset local_size = static_cast<MPI_Offset>(local.str().size());
    MPI_Offset rank_offset = 0;
    MPI_Exscan(&local_size, &rank_offset, 1, MPI_OFFSET, MPI_SUM, comm);
    if (rank == 0) rank_offset = 0;
    const int file_descriptor = open(file.str().c_str(), O_WRONLY);
    if (file_descriptor < 0 || !write_dump_at(file_descriptor, local.str(), atom_offset + rank_offset) || close(file_descriptor) != 0) ModuleBase::WARNING_QUIT("MD_func::dump_info", "cannot write MD_dump atoms.");
#else
    std::ofstream ofs(file.str().c_str(), step == 0 ? std::ios::trunc : std::ios::app);
    ofs << header.str() << local.str();
#endif
}

void get_mass_mbl(const UnitCell& unit_in,
                  double* allmass,
                  ModuleBase::Vector3<int>& frozen,
                  ModuleBase::Vector3<int>* ionmbl)
{
    int ion = 0;
    frozen.set(0, 0, 0);
    for (int it = 0; it < unit_in.ntype; it++)
    {
        for (int i = 0; i < unit_in.atoms[it].na; i++)
        {
            allmass[ion] = unit_in.atoms[it].mass / ModuleBase::AU_to_MASS;
            ionmbl[ion] = unit_in.atoms[it].mbl[i];
            if (ionmbl[ion].x == 0) {
                ++frozen.x;
}
            if (ionmbl[ion].y == 0) {
                ++frozen.y;
}
            if (ionmbl[ion].z == 0) {
                ++frozen.z;
}

            ion++;
        }
    }

    return;
}

double target_temp(const int& istep, const int& nstep, const double& tfirst, const double& tlast)
{
    assert(nstep>0);
    double delta = static_cast<double>(istep) / nstep;
    return tfirst + delta * (tlast - tfirst);
}

double current_temp(double& kinetic,
                    const int& natom,
                    const int& frozen_freedom,
                    const double* allmass,
                    const ModuleBase::Vector3<double>* vel)
{
    if (3 * natom == frozen_freedom)
    {
        kinetic = 0.0;
        return 0.0;
    }
    kinetic = kinetic_energy(natom, vel, allmass);
    return 2 * kinetic / (3 * natom - frozen_freedom);
}

double current_temp(double& kinetic,
                    const MDCell& mdcell,
                    const std::int64_t& frozen_freedom)
{
    kinetic = 0.0;
    for (std::size_t i = 0; i < mdcell.owned_atoms().size(); ++i)
    {
        const LocalAtom& atom = mdcell.owned_atoms()[i];
        kinetic += 0.5 * atom.mass * atom.vel.norm2();
    }
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, &kinetic, 1, MPI_DOUBLE, MPI_SUM, mdcell.communicator());
#endif
    const std::int64_t dof = 3 * mdcell.nat() - frozen_freedom;
    if (dof == 0) return 0.0;
    return 2.0 * kinetic / static_cast<double>(dof);
}

std::int64_t global_dof(const MDCell& mdcell, const std::int64_t& frozen_freedom)
{
    static_cast<void>(frozen_freedom);
    std::int64_t natom = mdcell.nlocal();
    std::int64_t local_frozen[3] = {0, 0, 0};
    for (int i = 0; i < mdcell.nlocal(); ++i)
    {
        const ModuleBase::Vector3<int>& mbl = mdcell.owned_atoms()[static_cast<std::size_t>(i)].mbl;
        if (mbl.x == 0) ++local_frozen[0];
        if (mbl.y == 0) ++local_frozen[1];
        if (mbl.z == 0) ++local_frozen[2];
    }
    std::int64_t global_frozen[3] = {local_frozen[0], local_frozen[1], local_frozen[2]};
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, &natom, 1, MPI_INT64_T, MPI_SUM, mdcell.communicator());
    MPI_Allreduce(MPI_IN_PLACE, global_frozen, 3, MPI_INT64_T, MPI_SUM, mdcell.communicator());
#endif
    std::int64_t total_frozen = global_frozen[0] + global_frozen[1] + global_frozen[2];
    if (global_frozen[0] == 0) ++total_frozen;
    if (global_frozen[1] == 0) ++total_frozen;
    if (global_frozen[2] == 0) ++total_frozen;
    return 3 * natom - total_frozen;
}

void temp_vector(const int& natom,
                 const ModuleBase::Vector3<double>* vel,
                 const double* allmass,
                 ModuleBase::matrix& t_vector)
{
    t_vector.create(3, 3);

    double t00 = 0.0;
    double t01 = 0.0;
    double t02 = 0.0;
    double t10 = 0.0;
    double t11 = 0.0;
    double t12 = 0.0;
    double t20 = 0.0;
    double t21 = 0.0;
    double t22 = 0.0;

#pragma omp parallel for reduction(+:t00, t01, t02, t10, t11, t12, t20, t21, t22) schedule(static) if (natom >= 256)
    for (int ion = 0; ion < natom; ++ion)
    {
        const double mass = allmass[ion];
        const double vx = vel[ion].x;
        const double vy = vel[ion].y;
        const double vz = vel[ion].z;

        t00 += mass * vx * vx;
        t01 += mass * vx * vy;
        t02 += mass * vx * vz;
        t10 += mass * vy * vx;
        t11 += mass * vy * vy;
        t12 += mass * vy * vz;
        t20 += mass * vz * vx;
        t21 += mass * vz * vy;
        t22 += mass * vz * vz;
    }

    t_vector(0, 0) = t00;
    t_vector(0, 1) = t01;
    t_vector(0, 2) = t02;
    t_vector(1, 0) = t10;
    t_vector(1, 1) = t11;
    t_vector(1, 2) = t12;
    t_vector(2, 0) = t20;
    t_vector(2, 1) = t21;
    t_vector(2, 2) = t22;

    return;
}


void current_md_info(const int& my_rank, const std::string& file_dir, int& md_step, double& temperature)
{
    bool ok = true;

    if (my_rank == 0)
    {
        std::stringstream ssc;
        ssc << file_dir << "Restart_md.txt";
        std::ifstream file(ssc.str().c_str());

        if (!file)
        {
            ok = false;
        }

        if (ok)
        {
            file >> md_step >> temperature;
            file.close();
        }
    }

#ifdef __MPI
    MPI_Bcast(&ok, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
#endif

    if (!ok)
    {
        ModuleBase::WARNING_QUIT("current_md_info", "no Restart_md.txt!");
    }

#ifdef __MPI
    MPI_Bcast(&md_step, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&temperature, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

    return;
}

void current_md_info(const MDCell& mdcell, const std::string& file_dir, int& md_step, double& temperature)
{
    bool ok = true;

#ifdef __MPI
    const int rank = mdcell.mpi_rank();
#else
    const int rank = 0;
#endif
    if (rank == 0)
    {
        std::stringstream ssc;
        ssc << file_dir << "Restart_md.txt";
        std::ifstream file(ssc.str().c_str());
        if (!file)
        {
            ok = false;
        }
        if (ok)
        {
            file >> md_step >> temperature;
        }
    }

#ifdef __MPI
    MPI_Bcast(&ok, 1, MPI_C_BOOL, 0, mdcell.communicator());
#endif
    if (!ok)
    {
        ModuleBase::WARNING_QUIT("current_md_info", "no Restart_md.txt!");
    }
#ifdef __MPI
    MPI_Bcast(&md_step, 1, MPI_INT, 0, mdcell.communicator());
    MPI_Bcast(&temperature, 1, MPI_DOUBLE, 0, mdcell.communicator());
#endif
}

} // namespace MD_func
