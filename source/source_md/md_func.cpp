#include "md_func.h"
#include "source_cell/domain_decomposition.h"

#include "source_base/global_variable.h"
#include "source_base/parallel_cell.h"
#include "source_base/timer.h"
#include "source_io/module_output/output_log.h"

#include <algorithm>
#include <fstream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <utility>


namespace MD_func
{

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

void init_vel(MDCell& mdcell,
              const bool& init_vel,
              const bool& restart,
              double& temperature,
              std::int64_t& frozen_freedom)
{
    if (mdcell.has_backing_unitcell())
    {
        const UnitCell& unit = mdcell.backing_unitcell();
        std::vector<double> masses(static_cast<std::size_t>(unit.nat), 0.0);
        std::vector<ModuleBase::Vector3<int> > mobile(static_cast<std::size_t>(unit.nat));
        std::vector<ModuleBase::Vector3<double> > velocities(static_cast<std::size_t>(unit.nat));
        ModuleBase::Vector3<int> frozen;
        get_mass_mbl(unit, masses.data(), frozen, mobile.data());
        frozen_freedom = frozen.x + frozen.y + frozen.z;
        if (frozen.x == 0) ++frozen_freedom;
        if (frozen.y == 0) ++frozen_freedom;
        if (frozen.z == 0) ++frozen_freedom;

        if (init_vel)
        {
            read_vel(unit, velocities.data());
            double kinetic = 0.0;
            const double current = current_temp(kinetic, unit.nat, frozen_freedom, masses.data(), velocities.data());
            if (!restart && temperature >= 0.0 && current > 0.0)
            {
                rescale_vel(unit.nat, temperature, masses.data(), frozen_freedom, velocities.data());
            }
            else if (!restart && temperature < 0.0)
            {
                temperature = current;
            }
        }
        else
        {
#ifdef __MPI
            const int rank = mdcell.mpi_rank();
#else
            const int rank = 0;
#endif
            rand_vel(unit.nat,
                     temperature,
                     masses.data(),
                     static_cast<int>(frozen_freedom),
                     frozen,
                     mobile.data(),
                     rank,
                     velocities.data());
        }

        std::vector<int> offsets(static_cast<std::size_t>(unit.ntype + 1), 0);
        for (int it = 0; it < unit.ntype; ++it)
        {
            offsets[static_cast<std::size_t>(it + 1)] = offsets[static_cast<std::size_t>(it)] + unit.atoms[it].na;
        }
        for (LocalAtom& atom : mdcell.owned_atoms())
        {
            atom.vel = velocities[static_cast<std::size_t>(offsets[static_cast<std::size_t>(atom.type)] + atom.type_index)];
        }
        return;
    }

    std::vector<LocalAtom>& atoms = mdcell.owned_atoms();
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
    if (init_vel)
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
}

void force_virial(ModuleESolver::ESolver* p_esolver,
                  const int& istep,
                  MDCell& mdcell,
                  DomainDecomposition& decomp,
                  double& potential,
                  const bool& cal_stress,
                  ModuleBase::matrix& virial,
                  const bool& md_out_force)
{
    ModuleBase::TITLE("MD_func", "force_virial");
    ModuleBase::timer::start("MD_func", "force_virial");
    if (!mdcell.has_backing_unitcell())
    {
        decomp.prepare_neighbors(mdcell);
        p_esolver->runner(static_cast<BaseCell&>(mdcell), istep);
        decomp.accumulate_ghost_forces(mdcell);
        potential = 0.5 * p_esolver->cal_energy();
        for (LocalAtom& atom : mdcell.owned_atoms()) atom.force *= 0.5;
        if (md_out_force)
        {
            ModuleIO::print_force(GlobalV::ofs_running, mdcell, "TOTAL-FORCE (eV/Angstrom)");
        }
        if (cal_stress) { p_esolver->cal_stress(static_cast<BaseCell&>(mdcell), virial); virial *= 0.5; }
    }
    else
    {
        UnitCell& ucell = mdcell.backing_unitcell();
        std::vector<std::vector<ModuleBase::Vector3<double>>> backing_velocities(
            static_cast<std::size_t>(ucell.ntype));
        for (int it = 0; it < ucell.ntype; ++it)
        {
            backing_velocities[static_cast<std::size_t>(it)] = ucell.atoms[it].vel;
        }
        // Step zero evaluates the initial coordinates before the first MD drift.
        const bool positions_updated = istep > 0;
        mdcell.sync_backing_unitcell(positions_updated);
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
        for (LocalAtom& atom : mdcell.owned_atoms()) { const int iat=offsets[atom.type]+atom.type_index; atom.force.set(full_force(iat,0),full_force(iat,1),full_force(iat,2)); }
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

namespace
{
void write_dump_rows(std::ostream& output,
                     const std::string& text,
                     const std::int64_t begin,
                     const std::int64_t end)
{
    // Spatial ownership is unrelated to the stable type-major atom index.
    std::map<std::int64_t, std::string> rows;
    std::istringstream input(text);
    std::string line;
    while (std::getline(input, line))
    {
        std::istringstream record(line);
        std::int64_t index = 0;
        if (!(record >> index) || index < begin || index >= end)
        {
            ModuleBase::WARNING_QUIT("MD_func::dump_info", "invalid MD_dump atom index.");
        }
        const bool inserted = rows.emplace(index, line).second;
        if (!inserted)
        {
            ModuleBase::WARNING_QUIT("MD_func::dump_info", "duplicate MD_dump atom index.");
        }
    }
    const std::size_t expected_count = static_cast<std::size_t>(end - begin);
    if (rows.size() != expected_count)
    {
        ModuleBase::WARNING_QUIT("MD_func::dump_info", "missing MD_dump atoms.");
    }
    for (const auto& row : rows)
    {
        output << row.second << "\n";
    }
}
} // namespace

void dump_info(const int& step,
               const std::string& global_out_dir,
               const MDCell& mdcell,
               const MD_para& mdp,
               const bool cal_stress,
               const ModuleBase::matrix& virial)
{
    const int max_batch_bytes = 1024 * 1024;
    dump_info(step, global_out_dir, mdcell, mdp, cal_stress, virial, max_batch_bytes);
}

void dump_info(const int& step,
               const std::string& global_out_dir,
               const MDCell& mdcell,
               const MD_para& mdp,
               const bool cal_stress,
               const ModuleBase::matrix& virial,
               const int max_batch_bytes)
{
    if (max_batch_bytes <= 0)
    {
        throw std::invalid_argument("MD_dump batch budget must be positive.");
    }
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
    if (cal_stress && mdp.dump_virial)
    {
        header << "VIRIAL (kbar)\n";
        for (int i = 0; i < 3; ++i)
            header << "  " << virial(i, 0) * unit_virial << "  " << virial(i, 1) * unit_virial << "  " << virial(i, 2) * unit_virial << "\n";
    }
    header << "INDEX    LABEL    POSITION (Angstrom)";
    if (mdp.dump_force) header << "    FORCE (eV/Angstrom)";
    if (mdp.dump_vel) header << "    VELOCITY (Angstrom/fs)";
    header << "\n";
    std::vector<std::int64_t> type_offsets(mdcell.type_atom_counts().size() + 1, 0);
    for (std::size_t it = 0; it < mdcell.type_atom_counts().size(); ++it)
    {
        type_offsets[it + 1] = type_offsets[it] + mdcell.type_atom_counts()[it];
    }
    const auto format_atom = [&](const LocalAtom& atom) -> std::string
    {
        std::ostringstream row;
        row << std::fixed << std::setprecision(12);
        row << "  " << type_offsets[static_cast<std::size_t>(atom.type)] + atom.type_index
            << "  " << mdcell.type_labels()[static_cast<std::size_t>(atom.type)]
            << "  " << atom.cart.x * unit_pos << "  " << atom.cart.y * unit_pos << "  " << atom.cart.z * unit_pos;
        if (mdp.dump_force)
            row << "  " << atom.force.x * unit_force << "  " << atom.force.y * unit_force << "  " << atom.force.z * unit_force;
        if (mdp.dump_vel)
            row << "  " << atom.vel.x * unit_vel << "  " << atom.vel.y * unit_vel << "  " << atom.vel.z * unit_vel;
        row << "\n";
        return row.str();
    };
    std::vector<std::pair<std::int64_t, const LocalAtom*>> ordered_atoms;
    const std::size_t local_atom_count = mdcell.owned_atoms().size();
    ordered_atoms.reserve(local_atom_count);
    int local_row_bytes = 0;
    int local_overflow = 0;
    const std::size_t maximum_count = std::numeric_limits<int>::max();
    for (const LocalAtom& atom : mdcell.owned_atoms())
    {
        const std::int64_t index = type_offsets[static_cast<std::size_t>(atom.type)] + atom.type_index;
        const LocalAtom* atom_pointer = &atom;
        ordered_atoms.emplace_back(index, atom_pointer);
        const std::string row = format_atom(atom);
        if (row.size() > maximum_count)
        {
            local_overflow = 1;
        }
        else
        {
            const int row_bytes = static_cast<int>(row.size());
            local_row_bytes = std::max(local_row_bytes, row_bytes);
        }
    }
    std::sort(ordered_atoms.begin(), ordered_atoms.end(),
              [](const std::pair<std::int64_t, const LocalAtom*>& left,
                 const std::pair<std::int64_t, const LocalAtom*>& right)
              { return left.first < right.first; });
    ModuleBase::CommunicationDomain domain;
#ifdef __MPI
    domain.initialize(mdcell.communicator());
#endif
    if (domain.max(local_overflow) != 0)
    {
        throw std::length_error("A single MD_dump row exceeds the MPI count range.");
    }
    const int maximum_row_bytes = domain.max(local_row_bytes);
    const int row_bytes = std::max(maximum_row_bytes, 1);
    const int batch_bytes = std::max(max_batch_bytes, row_bytes);
    const std::int64_t batch_atoms = batch_bytes / row_bytes;
    std::ofstream ofs;
    if (domain.rank() == 0)
    {
        const std::string filename = file.str();
        const std::ios_base::openmode mode = step == 0 ? std::ios::trunc : std::ios::app;
        ofs.open(filename.c_str(), mode);
        if (!ofs)
        {
            ModuleBase::WARNING_QUIT("MD_func::dump_info", "cannot open MD_dump.");
        }
        ofs << header.str();
    }
    std::size_t local_index = 0;
    for (std::int64_t begin = 0; begin < type_offsets.back();)
    {
        // Each global-index range contains at most batch_atoms rows across all ranks.
        // Counts and displacements are bounded per batch, independent of frame size.
        const std::int64_t remaining = type_offsets.back() - begin;
        const std::int64_t count = std::min(batch_atoms, remaining);
        const std::int64_t end = begin + count;
        std::ostringstream local;
        while (local_index < ordered_atoms.size() && ordered_atoms[local_index].first < end)
        {
            const LocalAtom& atom = *ordered_atoms[local_index].second;
            local << format_atom(atom);
            ++local_index;
        }
        const std::string local_text = local.str();
        const std::string gathered = domain.gather_text(local_text, batch_bytes);
        if (domain.rank() == 0)
        {
            write_dump_rows(ofs, gathered, begin, end);
        }
        begin = end;
    }
    if (domain.rank() == 0)
    {
        ofs << "\n\n";
        ofs.close();
        if (!ofs)
        {
            ModuleBase::WARNING_QUIT("MD_func::dump_info", "cannot write MD_dump.");
        }
    }
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

std::int64_t global_dof(const MDCell& mdcell)
{
    std::int64_t local_frozen[3] = {0, 0, 0};
    for (int i = 0; i < mdcell.owned_atoms().size(); ++i)
    {
        const ModuleBase::Vector3<int>& mbl = mdcell.owned_atoms()[static_cast<std::size_t>(i)].mbl;
        if (mbl.x == 0) ++local_frozen[0];
        if (mbl.y == 0) ++local_frozen[1];
        if (mbl.z == 0) ++local_frozen[2];
    }
    std::int64_t global_frozen[3] = {local_frozen[0], local_frozen[1], local_frozen[2]};
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, global_frozen, 3, MPI_INT64_T, MPI_SUM, mdcell.communicator());
#endif
    std::int64_t total_frozen = global_frozen[0] + global_frozen[1] + global_frozen[2];
    if (global_frozen[0] == 0) ++total_frozen;
    if (global_frozen[1] == 0) ++total_frozen;
    if (global_frozen[2] == 0) ++total_frozen;
    return 3 * mdcell.nat() - total_frozen;
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
