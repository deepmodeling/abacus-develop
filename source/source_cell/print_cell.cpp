#include <regex>
#include <cassert>
#include <algorithm>
#include <cstdint>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>

#include "print_cell.h"
#include "source_cell/mdcell.h"
#include "source_base/formatter.h"
#include "source_base/parallel_cell.h"
#include "source_base/tool_title.h"
#include "source_base/tool_quit.h"
#include "source_base/global_variable.h"
#include "source_base/output.h"

#ifdef __MPI
#include <mpi.h>
#endif

namespace unitcell
{
    void print_tau(Atom* atoms,
                   const std::string& Coordinate,
                   const int ntype,
                   const double lat0,
                   std::ofstream &ofs)
    {
        ModuleBase::TITLE("UnitCell", "print_tau");
        // assert (direct || Coordinate == "Cartesian" || Coordinate == "Cartesian_angstrom"); // this line causes abort in unittest ReadAtomPositionsCACXY.
        // previously there are two if-statements, the first is `if(Coordinate == "Direct")` and the second is `if(Coordinate == "Cartesian" || Coordiante == "Cartesian_angstrom")`
        // however the Coordinate can also be value among Cartesian_angstrom_center_xy, Cartesian_angstrom_center_xz, Cartesian_angstrom_center_yz and Cartesian_angstrom_center_xyz

        // if Coordinate has value one of them, this print_tau will not print anything.
        std::regex pattern("Direct|Cartesian(_angstrom)?(_center_(xy|xz|yz|xyz))?");
        assert(std::regex_search(Coordinate, pattern));
        bool direct = (Coordinate == "Direct");

        //----------------------
        // print atom positions
        //----------------------
        std::string table;
        table += direct? " DIRECT COORDINATES\n": FmtCore::format(" CARTESIAN COORDINATES ( UNIT = %15.8f Bohr )\n", lat0);
        table += FmtCore::format("%5s%19s%19s%19s%8s\n", "atom", "x", "y", "z", "mag");
        for(int it = 0; it < ntype; it++)
        {
            for (int ia = 0; ia < atoms[it].na; ia++)
            {
                const double& x = direct? atoms[it].taud[ia].x: atoms[it].tau[ia].x;
                const double& y = direct? atoms[it].taud[ia].y: atoms[it].tau[ia].y;
                const double& z = direct? atoms[it].taud[ia].z: atoms[it].tau[ia].z;
                table += FmtCore::format("%5s%19.12f%19.12f%19.12f%8.4f\n", 
                                        atoms[it].label, 
                                        x, 
                                        y, 
                                        z, 
                                        atoms[it].mag[ia]); 
            }
        }
        table += "\n";
        ofs << table; 


        // print velocities
        ofs << " ATOMIC VELOCITIES" << std::endl;
        ofs << std::setprecision(12);
        ofs << std::setw(5) << "atom" 
            << std::setw(19) << "vx" 
            << std::setw(19) << "vy" 
            << std::setw(19) << "vz"
            << std::endl;
 
        for(int it = 0; it < ntype; it++)
        {
            for (int ia = 0; ia < atoms[it].na; ia++)
            {
                ofs << std::setw(5) << atoms[it].label;
                ofs << " " << std::setw(18) << atoms[it].vel[ia].x;
                ofs << " " << std::setw(18) << atoms[it].vel[ia].y;
                ofs << " " << std::setw(18) << atoms[it].vel[ia].z;
                ofs << std::endl;
            }
        }
        ofs << std::endl;
        ofs << std::setprecision(6); // return to 6, as original


        return;
    }

    void print_stru_file(const UnitCell& ucell,
                         const Atom*     atoms,
                         const ModuleBase::Matrix3& latvec,
                         const std::string& fn,
                         const std::string& header,
                         const int& nspin,
                         const bool& vel,
                         const bool& magmom,
                         const bool& orb,
                         const bool& dpks_desc,
                         const int& iproc,
                         const ModuleBase::matrix& force,
                         const bool& has_force,
                         const bool  append)
    {
        ModuleBase::TITLE("UnitCell","print_stru_file");
        if (iproc != 0)
        {
            return; // old: if(GlobalV::MY_RANK != 0) return;
        }
        // optional header comments
        std::string str;
        if (!header.empty())
        {
            str = header;
        }
        // ATOMIC_SPECIES
        str += "ATOMIC_SPECIES\n";
        for(int it=0; it<ucell.ntype; it++)
        { 
            str += FmtCore::format("%s %8.4f %s %s\n", 
                                    ucell.atoms[it].label, 
                                    ucell.atoms[it].mass, 
                                    ucell.pseudo_fn[it], 
                                    ucell.pseudo_type[it]); 
        }
        // NUMERICAL_ORBITAL
        if(orb)
        {
            str += "\nNUMERICAL_ORBITAL\n";
            for(int it = 0; it < ucell.ntype; it++) 
            { 
                str += ucell.orbital_fn[it] + "\n"; 
            }
        }
        // NUMERICAL_DESCRIPTOR
        if(dpks_desc) 
        { 
            str += "\nNUMERICAL_DESCRIPTOR\n" + ucell.descriptor_file + "\n"; 
        }
        // LATTICE_CONSTANT: fixed to one Angstrom expressed in Bohr, so that the
        // lattice vectors below can be written directly in Angstrom.
        const double lat0_angstrom = 1.0 / ModuleBase::BOHR_TO_A;
        str += "\nLATTICE_CONSTANT\n"
             + FmtCore::format("%-.10f", lat0_angstrom)
             + " # in Bohr (= 1 Angstrom); lattice vectors below are in Angstrom\n";
        // LATTICE_VECTORS: internal vectors are dimensionless multiples of ucell.lat0;
        // multiply by ucell.lat0 * BOHR_TO_A to get the physical vectors in Angstrom.
        const double lat_scale = ucell.lat0 * ModuleBase::BOHR_TO_A;
        str += "\nLATTICE_VECTORS # in Angstrom\n";
        str += FmtCore::format("%.16f %.16f %.16f\n",
                               latvec.e11 * lat_scale, latvec.e12 * lat_scale, latvec.e13 * lat_scale);
        str += FmtCore::format("%.16f %.16f %.16f\n",
                               latvec.e21 * lat_scale, latvec.e22 * lat_scale, latvec.e23 * lat_scale);
        str += FmtCore::format("%.16f %.16f %.16f\n",
                               latvec.e31 * lat_scale, latvec.e32 * lat_scale, latvec.e33 * lat_scale);
        // ATOMIC_POSITIONS
        str += "\nATOMIC_POSITIONS\n";
        int nat_ = 0; // counter iat, for printing out Mulliken magmom who is indexed by iat
        // Force output is controlled by the explicit has_force flag rather than
        // inferred from matrix dimensions, so a zero-allocated matrix is not
        // mistaken for computed forces.
        if (has_force && (force.nr != ucell.nat || force.nc != 3))
        {
            ModuleBase::WARNING_QUIT("print_stru_file",
                                     "has_force is true but force matrix dimensions do not match nat x 3");
        }
        // Positions are always written as Cartesian_angstrom so that the
        // coordinate format is independent of force availability; has_force
        // only controls the per-atom f fields and the unit annotation.
        const std::string unit_note = has_force ? " # positions in Angstrom, forces in eV/Angstrom\n"
                                                : " # positions in Angstrom\n";
        str += "Cartesian_angstrom" + unit_note;
        // Internal Cartesian tau is in units of lat0 (Bohr); convert to Angstrom.
        const double pos_conv = ucell.lat0 * ModuleBase::BOHR_TO_A;
        const double force_conv = ModuleBase::Ry_to_eV / ModuleBase::BOHR_TO_A; // Ry/Bohr to eV/Angstrom
        for(int it = 0; it < ucell.ntype; it++)
        {
            str += "\n" + ucell.atoms[it].label + " #label\n";
            // Output real initial magnetism: for nspin=2 use mag[0], for nspin=4 use norm of m_loc_[0]
            double start_mag = ucell.magnet.start_mag[it];
            if (atoms[it].na > 0) {
                if (nspin == 2) {
                    start_mag = atoms[it].mag[0];
                } else if (nspin == 4) {
                    start_mag = std::sqrt(std::pow(atoms[it].m_loc_[0].x, 2)
                                        + std::pow(atoms[it].m_loc_[0].y, 2)
                                        + std::pow(atoms[it].m_loc_[0].z, 2));
                }
            }
            str += FmtCore::format("%.4f #magnetism (default, overridden by per-atom mag below)\n", start_mag);
            str += FmtCore::format("%d #number of atoms\n", atoms[it].na);
            for(int ia = 0; ia < atoms[it].na; ia++)
            {
                // output position
                const double& x = atoms[it].tau[ia].x;
                const double& y = atoms[it].tau[ia].y;
                const double& z = atoms[it].tau[ia].z;
                str += FmtCore::format("%.10f %.10f %.10f", x*pos_conv, y*pos_conv, z*pos_conv);
                str += FmtCore::format(" m%2d%2d%2d", atoms[it].mbl[ia].x, atoms[it].mbl[ia].y, atoms[it].mbl[ia].z);
                if (vel) // output velocity
                {
                    str += FmtCore::format(" v %.10f %.10f %.10f", atoms[it].vel[ia].x, atoms[it].vel[ia].y, atoms[it].vel[ia].z);
                }
                if (has_force) // output force
                {
                    str += FmtCore::format(" f %.6f %.6f %.6f",
                                           force(nat_, 0)*force_conv,
                                           force(nat_, 1)*force_conv,
                                           force(nat_, 2)*force_conv);
                }
                if (nspin == 2) // output magnetic information
                {
                    if (magmom && !ucell.atom_mulliken.empty()) {
                        str += FmtCore::format(" mag %.4f", ucell.atom_mulliken[nat_][1]);
                    } else {
                        str += FmtCore::format(" mag %.4f", atoms[it].mag[ia]);
                    }
                }
                else if (nspin == 4) // output magnetic information
                {
                    if (magmom && !ucell.atom_mulliken.empty()) {
                        str += FmtCore::format(" mag %.4f %.4f %.4f",
                                                ucell.atom_mulliken[nat_][1],
                                                ucell.atom_mulliken[nat_][2],
                                                ucell.atom_mulliken[nat_][3]);
                    } else {
                        str += FmtCore::format(" mag %.4f %.4f %.4f",
                                                atoms[it].m_loc_[ia].x,
                                                atoms[it].m_loc_[ia].y,
                                                atoms[it].m_loc_[ia].z);
                    }
                }
                str += "\n";
                nat_++;
            }
        }
        std::ios_base::openmode mode = std::ios::out;
        if (append)
        {
            mode |= std::ios::app;
        }
        std::ofstream ofs(fn.c_str(), mode);
        ofs << str;
        ofs.close();
        return;
    }

    void print_cell(const UnitCell& ucell, std::ofstream& ofs)
    {
        ModuleBase::GlobalFunc::OUT(ofs, "print_unitcell()");

        ModuleBase::GlobalFunc::OUT(ofs, "latName", ucell.latName);
        ModuleBase::GlobalFunc::OUT(ofs, "ntype", ucell.ntype);
        ModuleBase::GlobalFunc::OUT(ofs, "nat", ucell.nat);
        ModuleBase::GlobalFunc::OUT(ofs, "lat0", ucell.lat0);
        ModuleBase::GlobalFunc::OUT(ofs, "lat0_angstrom", ucell.lat0_angstrom);
        ModuleBase::GlobalFunc::OUT(ofs, "tpiba", ucell.tpiba);
        ModuleBase::GlobalFunc::OUT(ofs, "omega", ucell.omega);

        output::printM3(ofs, "Lattices Vector (R) : ", ucell.latvec);
        output::printM3(ofs, "Supercell lattice vector : ", ucell.latvec_supercell);
        output::printM3(ofs, "Reciprocal lattice Vector (G): ", ucell.G);
        output::printM3(ofs, "GGT : ", ucell.GGT);

        ofs << std::endl;
        return;
    }
}

namespace
{
std::string mdcell_stru_header(const MDCell& cell, const StruMeta& metadata)
{
    std::ostringstream output;
    output << std::fixed << std::setprecision(10);
    output << "ATOMIC_SPECIES\n";
    for (std::size_t it = 0; it < metadata.species.size(); ++it)
    {
        const StruSpecies& species = metadata.species[it];
        output << cell.type_labels()[it] << " " << std::setprecision(4) << cell.type_masses()[it] << std::setprecision(10);
        if (!species.pseudo_file.empty()) output << " " << species.pseudo_file;
        if (!species.pseudo_type.empty()) output << " " << species.pseudo_type;
        output << "\n";
    }
    bool has_orbitals = false;
    for (std::size_t it = 0; it < metadata.species.size(); ++it)
        has_orbitals = has_orbitals || !metadata.species[it].orbital_file.empty();
    if (has_orbitals)
    {
        output << "\nNUMERICAL_ORBITAL\n";
        for (std::size_t it = 0; it < metadata.species.size(); ++it)
            output << metadata.species[it].orbital_file << "\n";
    }
    if (!metadata.descriptor_file.empty()) output << "\nNUMERICAL_DESCRIPTOR\n" << metadata.descriptor_file << "\n";
    const double lat0_angstrom = 1.0 / ModuleBase::BOHR_TO_A;
    const double lat_scale = cell.lat0() * ModuleBase::BOHR_TO_A;
    output << "\nLATTICE_CONSTANT\n" << lat0_angstrom << "\n\nLATTICE_VECTORS\n";
    const ModuleBase::Matrix3& lattice = cell.latvec();
    output << lattice.e11 * lat_scale << " " << lattice.e12 * lat_scale << " " << lattice.e13 * lat_scale << "\n";
    output << lattice.e21 * lat_scale << " " << lattice.e22 * lat_scale << " " << lattice.e23 * lat_scale << "\n";
    output << lattice.e31 * lat_scale << " " << lattice.e32 * lat_scale << " " << lattice.e33 * lat_scale << "\n";
    output << "\nATOMIC_POSITIONS\nCartesian_angstrom\n";
    return output.str();
}

std::string mdcell_type_header(const MDCell& cell, const StruMeta& metadata, const std::size_t it)
{
    const StruSpecies& species = metadata.species[it];
    std::ostringstream output;
    output << "\n" << cell.type_labels()[it] << " #label\n";
    output << std::fixed << std::setprecision(4) << species.start_mag << " #magnetism\n";
    output << cell.type_atom_counts()[it] << " #number of atoms\n";
    return output.str();
}

std::string format_mdcell_atom(const MDCell& cell, const LocalAtom& atom, const bool include_force)
{
    std::ostringstream output;
    const double position_scale = cell.lat0() * ModuleBase::BOHR_TO_A;
    output << std::fixed << std::setprecision(10)
           << atom.cart.x * position_scale << " " << atom.cart.y * position_scale << " "
           << atom.cart.z * position_scale
           << " m " << atom.mbl.x << " " << atom.mbl.y << " " << atom.mbl.z;
    output << " v " << atom.vel.x << " " << atom.vel.y << " " << atom.vel.z;
    if (include_force)
    {
        const double force_scale = ModuleBase::Hartree_to_eV * ModuleBase::ANGSTROM_AU;
        output << " f " << atom.force.x * force_scale << " " << atom.force.y * force_scale
               << " " << atom.force.z * force_scale;
    }
    return output.str();
}

void write_mdcell_rows(std::ostream& output, const std::string& text,
                       const std::int64_t begin, const std::int64_t end)
{
    std::map<std::int64_t, std::string> rows;
    std::istringstream input(text);
    std::string line;
    while (std::getline(input, line))
    {
        std::istringstream record(line);
        std::int64_t index = -1;
        if (!(record >> index) || index < begin || index >= end)
        {
            throw std::runtime_error("Invalid MDCell STRU atom index.");
        }
        const std::size_t separator = line.find('\t');
        if (separator == std::string::npos || !rows.emplace(index, line.substr(separator + 1)).second)
        {
            throw std::runtime_error("Duplicate or malformed MDCell STRU atom row.");
        }
    }
    if (rows.size() != static_cast<std::size_t>(end - begin))
    {
        throw std::runtime_error("Missing MDCell STRU atoms.");
    }
    for (const std::pair<const std::int64_t, std::string>& row : rows)
    {
        output << row.second << '\n';
    }
}
}

namespace unitcell
{
StruMeta make_stru_meta(const UnitCell& ucell)
{
    StruMeta metadata;
    metadata.species.resize(static_cast<std::size_t>(ucell.ntype));
    for (int it = 0; it < ucell.ntype; ++it)
    {
        StruSpecies& species = metadata.species[it];
        if (static_cast<std::size_t>(it) < ucell.pseudo_fn.size()) species.pseudo_file = ucell.pseudo_fn[it];
        if (static_cast<std::size_t>(it) < ucell.pseudo_type.size()) species.pseudo_type = ucell.pseudo_type[it];
        if (static_cast<std::size_t>(it) < ucell.orbital_fn.size()) species.orbital_file = ucell.orbital_fn[it];
        if (static_cast<std::size_t>(it) < ucell.magnet.start_mag.size()) species.start_mag = ucell.magnet.start_mag[it];
    }
    metadata.descriptor_file = ucell.descriptor_file;
    return metadata;
}
}

namespace mdcell
{
void print_stru_file(const MDCell& cell, const StruMeta& stru_meta, const std::string& fn,
                     const std::string& frame_comment, const bool append, const bool include_force)
{
    if (stru_meta.species.size() != cell.type_labels().size()
        || stru_meta.species.size() != cell.type_masses().size()
        || stru_meta.species.size() != cell.type_atom_counts().size())
    {
        throw std::runtime_error("MDCell STRU metadata does not match the MDCell type data.");
    }
    const std::string header = frame_comment + mdcell_stru_header(cell, stru_meta);
    ModuleBase::CommunicationDomain domain;
#ifdef __MPI
    domain.initialize(cell.communicator());
#endif
    std::vector<std::vector<std::pair<std::int64_t, std::string> > > local_rows(stru_meta.species.size());
    int invalid_identity = 0;
    int local_row_bytes = 0;
    for (const LocalAtom& atom : cell.owned_atoms())
    {
        if (atom.type < 0 || static_cast<std::size_t>(atom.type) >= stru_meta.species.size()
            || atom.type_index < 0 || atom.type_index >= cell.type_atom_counts()[atom.type])
        {
            invalid_identity = 1;
            continue;
        }
        const std::string row = std::to_string(atom.type_index) + "\t"
                                + format_mdcell_atom(cell, atom, include_force) + "\n";
        if (row.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        {
            invalid_identity = 1;
            continue;
        }
        local_row_bytes = std::max(local_row_bytes, static_cast<int>(row.size()));
        local_rows[atom.type].push_back(std::make_pair(atom.type_index, row));
    }
    if (domain.max(invalid_identity) != 0)
    {
        throw std::runtime_error("Invalid MDCell STRU atom identity or row size.");
    }
    for (auto& rows : local_rows)
    {
        std::sort(rows.begin(), rows.end(),
                  [](const std::pair<std::int64_t, std::string>& left,
                     const std::pair<std::int64_t, std::string>& right) { return left.first < right.first; });
    }
    const int largest_row_bytes = domain.max(local_row_bytes);
    const int row_bytes = std::max(largest_row_bytes, 1);
    const int batch_bytes = std::max(1024 * 1024, row_bytes);
    const std::int64_t batch_atoms = batch_bytes / row_bytes;
    std::ofstream output;
    if (domain.rank() == 0)
    {
        const std::ios_base::openmode mode = append ? std::ios::out | std::ios::app : std::ios::out;
        output.open(fn.c_str(), mode);
        if (output)
        {
            output << header;
        }
    }
    const int open_failed = domain.rank() == 0 && !output ? 1 : 0;
    if (domain.max(open_failed) != 0)
    {
        throw std::runtime_error("Unable to open MDCell STRU file: " + fn);
    }
    for (std::size_t it = 0; it < stru_meta.species.size(); ++it)
    {
        const std::string type_header = mdcell_type_header(cell, stru_meta, it);
        if (domain.rank() == 0)
        {
            output << type_header;
        }
        std::size_t local_index = 0;
        const std::vector<std::pair<std::int64_t, std::string> >& rows = local_rows[it];
        for (std::int64_t begin = 0; begin < cell.type_atom_counts()[it];)
        {
            const std::int64_t remaining = cell.type_atom_counts()[it] - begin;
            const std::int64_t count = std::min(batch_atoms, remaining);
            const std::int64_t end = begin + count;
            std::string local_text;
            while (local_index < rows.size() && rows[local_index].first < end)
            {
                local_text += rows[local_index].second;
                ++local_index;
            }
            const std::string gathered = domain.gather_text(local_text, batch_bytes);
            int invalid_rows = 0;
            if (domain.rank() == 0)
            {
                try
                {
                    write_mdcell_rows(output, gathered, begin, end);
                }
                catch (const std::runtime_error&)
                {
                    invalid_rows = 1;
                }
            }
            if (domain.max(invalid_rows) != 0)
            {
                throw std::runtime_error("Duplicate, missing or malformed MDCell STRU atom index.");
            }
            begin = end;
        }
    }
    if (domain.rank() == 0)
    {
        output.close();
    }
    const int write_failed = domain.rank() == 0 && !output ? 1 : 0;
    if (domain.max(write_failed) != 0)
    {
        throw std::runtime_error("Unable to write MDCell STRU file: " + fn);
    }
}
}
