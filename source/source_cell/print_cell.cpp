#include <regex>
#include <cassert>
#include <cerrno>
#include <cstring>
#include <fcntl.h>
#include <sstream>
#include <stdexcept>
#include <unistd.h>

#include "print_cell.h"
#include "source_cell/mdcell.h"
#include "source_base/formatter.h"
#include "source_base/tool_title.h"
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
                         const bool& direct,
                         const bool& vel,
                         const bool& magmom,
                         const bool& orb,
                         const bool& dpks_desc,
                         const int& iproc,
                         const ModuleBase::matrix& force)
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
        // If force is provided, output positions in Angstrom and forces in eV/Angstrom.
        // Fractional (Direct) positions are only emitted when no force is needed.
        const bool has_force = (force.nr == ucell.nat && force.nc == 3);
        const bool use_cartesian = has_force || !direct;
        const std::string scale = use_cartesian ? "Cartesian_angstrom" : "Direct";
        std::string unit_note = "\n";
        if (use_cartesian)
        {
            unit_note = has_force ? " # positions in Angstrom, forces in eV/Angstrom\n"
                                  : " # positions in Angstrom\n";
        }
        str += scale + unit_note;
        // Internal Cartesian tau is in units of lat0 (Bohr); convert to Angstrom.
        const double pos_conv = use_cartesian ? ucell.lat0 * ModuleBase::BOHR_TO_A : 1.0;
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
                const double& x = use_cartesian ? atoms[it].tau[ia].x : atoms[it].taud[ia].x;
                const double& y = use_cartesian ? atoms[it].tau[ia].y : atoms[it].taud[ia].y;
                const double& z = use_cartesian ? atoms[it].tau[ia].z : atoms[it].taud[ia].z;
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
        std::ofstream ofs(fn.c_str());
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
    output << "\nLATTICE_CONSTANT\n" << cell.lat0() << "\n\nLATTICE_VECTORS\n";
    const ModuleBase::Matrix3& lattice = cell.latvec();
    output << lattice.e11 << " " << lattice.e12 << " " << lattice.e13 << "\n";
    output << lattice.e21 << " " << lattice.e22 << " " << lattice.e23 << "\n";
    output << lattice.e31 << " " << lattice.e32 << " " << lattice.e33 << "\n";
    output << "\nATOMIC_POSITIONS\nCartesian\n";
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

std::string local_mdcell_atoms(const MDCell& cell, const std::size_t type)
{
    std::string output;
    for (std::size_t iat = 0; iat < cell.owned_atoms().size(); ++iat)
    {
        const LocalAtom& atom = cell.owned_atoms()[iat];
        if (atom.type == static_cast<int>(type))
        {
            std::ostringstream atom_output;
            atom_output << std::fixed << std::setprecision(10)
                        << atom.cart.x << " " << atom.cart.y << " " << atom.cart.z
                        << " m " << atom.mbl.x << " " << atom.mbl.y << " " << atom.mbl.z
                        << " v " << atom.vel.x << " " << atom.vel.y << " " << atom.vel.z << "\n";
            output += atom_output.str();
        }
    }
    return output;
}

#ifdef __MPI
bool write_at(const int file, const std::string& data, MPI_Offset offset)
{
    std::size_t written = 0;
    while (written < data.size())
    {
        const ssize_t count = pwrite(file,
                                     data.data() + written,
                                     data.size() - written,
                                     static_cast<off_t>(offset + written));
        if (count <= 0)
        {
            return false;
        }
        written += static_cast<std::size_t>(count);
    }
    return true;
}
#endif
}

namespace unitcell
{
StruMeta make_stru_meta(const UnitCell& ucell)
{
    StruMeta metadata;
    metadata.species.resize(static_cast<std::size_t>(ucell.ntype));
    for (int it = 0; it < ucell.ntype; ++it)
    {
        StruSpecies& species = metadata.species[static_cast<std::size_t>(it)];
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
void print_stru_file(const MDCell& cell, const StruMeta& stru_meta, const std::string& fn)
{
    if (stru_meta.species.size() != cell.type_labels().size()
        || stru_meta.species.size() != cell.type_masses().size()
        || stru_meta.species.size() != cell.type_atom_counts().size())
    {
        throw std::runtime_error("MDCell STRU metadata does not match the MDCell type data.");
    }
    const std::string header = mdcell_stru_header(cell, stru_meta);
#ifdef __MPI
    int rank = 0;
    const MPI_Comm comm = cell.communicator();
    MPI_Comm_rank(comm, &rank);
    int header_ok = 1;
    if (rank == 0)
    {
        const int header_file = open(fn.c_str(), O_CREAT | O_TRUNC | O_WRONLY, 0666);
        if (header_file < 0)
        {
            header_ok = 0;
        }
        else
        {
            const bool header_written = write_at(header_file, header, 0);
            const bool header_closed = close(header_file) == 0;
            header_ok = header_written && header_closed;
        }
    }
    MPI_Bcast(&header_ok, 1, MPI_INT, 0, comm);
    if (header_ok == 0)
    {
        throw std::runtime_error("Unable to create MDCell restart STRU file: " + fn + ": " + std::strerror(errno));
    }
    MPI_Barrier(comm);

    const int file = open(fn.c_str(), O_WRONLY);
    int file_ok = file >= 0 ? 1 : 0;
    int all_files_ok = 0;
    MPI_Allreduce(&file_ok, &all_files_ok, 1, MPI_INT, MPI_MIN, comm);
    if (all_files_ok == 0)
    {
        if (file >= 0) close(file);
        throw std::runtime_error("Unable to open MDCell restart STRU file: " + fn + ": " + std::strerror(errno));
    }

    MPI_Offset offset = static_cast<MPI_Offset>(header.size());
    for (std::size_t it = 0; it < stru_meta.species.size(); ++it)
    {
        const std::string type_header = mdcell_type_header(cell, stru_meta, it);
        int type_header_ok = 1;
        if (rank == 0) type_header_ok = write_at(file, type_header, offset) ? 1 : 0;
        MPI_Bcast(&type_header_ok, 1, MPI_INT, 0, comm);
        if (type_header_ok == 0)
        {
            close(file);
            throw std::runtime_error("Unable to write MDCell restart STRU type header: " + fn + ": " + std::strerror(errno));
        }
        offset += static_cast<MPI_Offset>(type_header.size());
        const std::string local_atoms = local_mdcell_atoms(cell, it);
        const MPI_Offset local_size = static_cast<MPI_Offset>(local_atoms.size());
        MPI_Offset type_size = 0;
        MPI_Offset rank_offset = 0;
        MPI_Allreduce(&local_size, &type_size, 1, MPI_OFFSET, MPI_SUM, comm);
        MPI_Exscan(&local_size, &rank_offset, 1, MPI_OFFSET, MPI_SUM, comm);
        if (rank == 0) rank_offset = 0;
        int atom_data_ok = 1;
        if (local_size > 0)
        {
            atom_data_ok = write_at(file, local_atoms, offset + rank_offset) ? 1 : 0;
        }
        int all_atom_data_ok = 0;
        MPI_Allreduce(&atom_data_ok, &all_atom_data_ok, 1, MPI_INT, MPI_MIN, comm);
        if (all_atom_data_ok == 0)
        {
            close(file);
            throw std::runtime_error("Unable to write MDCell restart STRU atom data: " + fn + ": " + std::strerror(errno));
        }
        MPI_Barrier(comm);
        offset += type_size;
    }
    if (close(file) != 0)
    {
        throw std::runtime_error("Unable to close MDCell restart STRU file: " + fn + ": " + std::strerror(errno));
    }
#else
    std::ofstream output(fn.c_str());
    output << header;
    for (std::size_t it = 0; it < stru_meta.species.size(); ++it)
        output << mdcell_type_header(cell, stru_meta, it) << local_mdcell_atoms(cell, it);
#endif
}
}
