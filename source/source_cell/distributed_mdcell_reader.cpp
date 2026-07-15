#include "source_cell/distributed_mdcell_reader.h"

#ifdef __MPI

#include "source_base/constants.h"
#include "source_base/vector3.h"
#include "source_cell/module_neighlist/domain_decomposition.h"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace
{
struct RootMetadata
{
    double lat0;
    double omega;
    ModuleBase::Matrix3 latvec;
    ModuleBase::Matrix3 gt;
    std::vector<std::string> labels;
    std::vector<double> masses;
    int nat;
};

struct PackedMdAtom
{
    double cart[3];
    double frac[3];
    double vel[3];
    int mbl[3];
    double mass;
    int type;
    int type_index;
    ModuleNeighList::GlobalAtomId global_id;
};

std::string trim_copy(const std::string& value)
{
    std::size_t begin = 0;
    while (begin < value.size() && std::isspace(static_cast<unsigned char>(value[begin])))
    {
        ++begin;
    }
    std::size_t end = value.size();
    while (end > begin && std::isspace(static_cast<unsigned char>(value[end - 1])))
    {
        --end;
    }
    return value.substr(begin, end - begin);
}

std::string strip_comment(const std::string& line)
{
    const std::size_t pos = line.find('#');
    return trim_copy(pos == std::string::npos ? line : line.substr(0, pos));
}

bool is_blank_or_comment(const std::string& line)
{
    return strip_comment(line).empty();
}

std::string next_data_line(std::ifstream& ifs, const char* context)
{
    std::string line;
    while (std::getline(ifs, line))
    {
        line = strip_comment(line);
        if (!line.empty())
        {
            return line;
        }
    }
    throw std::runtime_error(std::string("Unexpected EOF while reading ") + context + ".");
}

void expect_keyword(std::ifstream& ifs, const char* keyword)
{
    const std::string line = next_data_line(ifs, keyword);
    if (line != keyword)
    {
        throw std::runtime_error(std::string("Expected keyword '") + keyword + "', got '" + line + "'.");
    }
}

double parse_double(const std::string& token, const char* context)
{
    char* end = NULL;
    const double value = std::strtod(token.c_str(), &end);
    if (end == token.c_str() || *end != '\0')
    {
        throw std::runtime_error(std::string("Failed to parse double for ") + context + ": " + token);
    }
    return value;
}

int parse_int(const std::string& token, const char* context)
{
    char* end = NULL;
    const long value = std::strtol(token.c_str(), &end, 10);
    if (end == token.c_str() || *end != '\0')
    {
        throw std::runtime_error(std::string("Failed to parse int for ") + context + ": " + token);
    }
    return static_cast<int>(value);
}

ModuleBase::Vector3<double> wrap_fractional(const ModuleBase::Vector3<double>& frac)
{
    ModuleBase::Vector3<double> wrapped = frac;
    wrapped.x -= std::floor(wrapped.x);
    wrapped.y -= std::floor(wrapped.y);
    wrapped.z -= std::floor(wrapped.z);
    if (wrapped.x >= 1.0 - 1.0e-12 || wrapped.x < 1.0e-12) wrapped.x = 0.0;
    if (wrapped.y >= 1.0 - 1.0e-12 || wrapped.y < 1.0e-12) wrapped.y = 0.0;
    if (wrapped.z >= 1.0 - 1.0e-12 || wrapped.z < 1.0e-12) wrapped.z = 0.0;
    return wrapped;
}

PackedMdAtom pack_atom(const LocalAtom& atom)
{
    PackedMdAtom packed;
    packed.cart[0] = atom.cart.x;
    packed.cart[1] = atom.cart.y;
    packed.cart[2] = atom.cart.z;
    packed.frac[0] = atom.frac.x;
    packed.frac[1] = atom.frac.y;
    packed.frac[2] = atom.frac.z;
    packed.vel[0] = atom.vel.x;
    packed.vel[1] = atom.vel.y;
    packed.vel[2] = atom.vel.z;
    packed.mbl[0] = atom.mbl.x;
    packed.mbl[1] = atom.mbl.y;
    packed.mbl[2] = atom.mbl.z;
    packed.mass = atom.mass;
    packed.type = atom.type;
    packed.type_index = atom.type_index;
    packed.global_id = atom.global_id;
    return packed;
}

LocalAtom unpack_atom(const PackedMdAtom& packed, int owner_rank)
{
    return LocalAtom(ModuleBase::Vector3<double>(packed.cart[0], packed.cart[1], packed.cart[2]),
                     ModuleBase::Vector3<double>(packed.frac[0], packed.frac[1], packed.frac[2]),
                     ModuleBase::Vector3<double>(packed.vel[0], packed.vel[1], packed.vel[2]),
                     ModuleBase::Vector3<double>(0.0, 0.0, 0.0),
                     ModuleBase::Vector3<int>(packed.mbl[0], packed.mbl[1], packed.mbl[2]),
                     packed.mass,
                     packed.type,
                     packed.type_index,
                     packed.global_id,
                     owner_rank,
                     false);
}

void broadcast_metadata(RootMetadata& metadata, MPI_Comm comm)
{
    MPI_Bcast(&metadata.lat0, 1, MPI_DOUBLE, 0, comm);
    MPI_Bcast(&metadata.omega, 1, MPI_DOUBLE, 0, comm);
    MPI_Bcast(&metadata.latvec.e11, 9, MPI_DOUBLE, 0, comm);
    MPI_Bcast(&metadata.gt.e11, 9, MPI_DOUBLE, 0, comm);
    MPI_Bcast(&metadata.nat, 1, MPI_INT, 0, comm);

    int ntype = static_cast<int>(metadata.labels.size());
    MPI_Bcast(&ntype, 1, MPI_INT, 0, comm);
    if (metadata.labels.size() != static_cast<std::size_t>(ntype))
    {
        metadata.labels.resize(static_cast<std::size_t>(ntype));
    }
    if (metadata.masses.size() != static_cast<std::size_t>(ntype))
    {
        metadata.masses.resize(static_cast<std::size_t>(ntype));
    }
    if (ntype > 0)
    {
        MPI_Bcast(&metadata.masses[0], ntype, MPI_DOUBLE, 0, comm);
    }
    for (int it = 0; it < ntype; ++it)
    {
        int len = static_cast<int>(metadata.labels[static_cast<std::size_t>(it)].size());
        MPI_Bcast(&len, 1, MPI_INT, 0, comm);
        if (metadata.labels[static_cast<std::size_t>(it)].size() != static_cast<std::size_t>(len))
        {
            metadata.labels[static_cast<std::size_t>(it)].assign(static_cast<std::size_t>(len), '\0');
        }
        if (len > 0)
        {
            MPI_Bcast(&metadata.labels[static_cast<std::size_t>(it)][0], len, MPI_CHAR, 0, comm);
        }
    }
}

RootMetadata parse_root_metadata(std::ifstream& ifs)
{
    RootMetadata metadata;
    metadata.lat0 = 1.0;
    metadata.omega = 0.0;
    metadata.nat = 0;

    expect_keyword(ifs, "ATOMIC_SPECIES");
    while (true)
    {
        const std::streampos mark = ifs.tellg();
        std::string line = next_data_line(ifs, "ATOMIC_SPECIES body");
        if (line == "LATTICE_CONSTANT")
        {
            ifs.seekg(mark);
            break;
        }
        std::istringstream iss(line);
        std::string label;
        std::string mass_token;
        iss >> label >> mass_token;
        if (label.empty() || mass_token.empty())
        {
            throw std::runtime_error("Invalid ATOMIC_SPECIES line: " + line);
        }
        metadata.labels.push_back(label);
        metadata.masses.push_back(parse_double(mass_token, "atomic mass"));
    }

    expect_keyword(ifs, "LATTICE_CONSTANT");
    metadata.lat0 = parse_double(next_data_line(ifs, "LATTICE_CONSTANT value"), "lattice constant");

    expect_keyword(ifs, "LATTICE_VECTORS");
    for (int row = 0; row < 3; ++row)
    {
        std::istringstream iss(next_data_line(ifs, "LATTICE_VECTORS row"));
        double x = 0.0;
        double y = 0.0;
        double z = 0.0;
        iss >> x >> y >> z;
        if (!iss)
        {
            throw std::runtime_error("Invalid LATTICE_VECTORS row.");
        }
        if (row == 0) { metadata.latvec.e11 = x; metadata.latvec.e12 = y; metadata.latvec.e13 = z; }
        if (row == 1) { metadata.latvec.e21 = x; metadata.latvec.e22 = y; metadata.latvec.e23 = z; }
        if (row == 2) { metadata.latvec.e31 = x; metadata.latvec.e32 = y; metadata.latvec.e33 = z; }
    }
    metadata.gt = metadata.latvec.Inverse();
    metadata.omega = std::abs(metadata.latvec.Det()) * metadata.lat0 * metadata.lat0 * metadata.lat0;
    return metadata;
}

void distribute_atoms_from_root(std::ifstream& ifs,
                                const RootMetadata& metadata,
                                MPI_Comm comm,
                                double cutoff_bohr,
                                double skin_bohr,
                                std::vector<LocalAtom>& owned_atoms)
{
    int rank = 0;
    int size = 1;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    DomainDecomposition decomp;
    decomp.init(comm, metadata.latvec, metadata.lat0, cutoff_bohr, skin_bohr);

    if (rank == 0)
    {
        expect_keyword(ifs, "ATOMIC_POSITIONS");
        const std::string coord_type = next_data_line(ifs, "ATOMIC_POSITIONS type");
        const bool is_cartesian = coord_type == "Cartesian";
        const bool is_direct = coord_type == "Direct";
        if (!is_cartesian && !is_direct)
        {
            throw std::runtime_error("Only Direct and Cartesian ATOMIC_POSITIONS are supported for distributed LJ MD.");
        }

        std::vector<std::vector<PackedMdAtom> > per_rank(static_cast<std::size_t>(size));
        ModuleNeighList::GlobalAtomId global_id = 0;
        for (std::size_t it = 0; it < metadata.labels.size(); ++it)
        {
            const std::string label = next_data_line(ifs, "atom label");
            if (label != metadata.labels[it])
            {
                throw std::runtime_error("ATOMIC_POSITIONS label order does not match ATOMIC_SPECIES.");
            }
            next_data_line(ifs, "magnetism");
            const int nat_type = parse_int(next_data_line(ifs, "atom count"), "atom count");
            for (int ia = 0; ia < nat_type; ++ia)
            {
                std::istringstream iss(next_data_line(ifs, "atom line"));
                double c1 = 0.0;
                double c2 = 0.0;
                double c3 = 0.0;
                iss >> c1 >> c2 >> c3;
                if (!iss)
                {
                    throw std::runtime_error("Invalid atomic coordinate line.");
                }

                ModuleBase::Vector3<double> frac;
                ModuleBase::Vector3<double> cart;
                if (is_cartesian)
                {
                    cart.set(c1, c2, c3);
                    frac = wrap_fractional(cart * metadata.gt);
                    cart = frac * metadata.latvec;
                }
                else
                {
                    frac = wrap_fractional(ModuleBase::Vector3<double>(c1, c2, c3));
                    cart = frac * metadata.latvec;
                }

                ModuleBase::Vector3<int> mbl(1, 1, 1);
                ModuleBase::Vector3<double> vel(0.0, 0.0, 0.0);
                std::string token;
                while (iss >> token)
                {
                    if (token == "m")
                    {
                        std::string mx;
                        std::string my;
                        std::string mz;
                        iss >> mx >> my >> mz;
                        if (!iss)
                        {
                            throw std::runtime_error("Invalid move flag record in STRU.");
                        }
                        mbl.set(parse_int(mx, "move flag x"),
                                parse_int(my, "move flag y"),
                                parse_int(mz, "move flag z"));
                    }
                    else if (token == "v" || token == "vel" || token == "velocity")
                    {
                        std::string vx;
                        std::string vy;
                        std::string vz;
                        iss >> vx >> vy >> vz;
                        if (!iss)
                        {
                            throw std::runtime_error("Invalid velocity record in STRU.");
                        }
                        vel.set(parse_double(vx, "velocity x"),
                                parse_double(vy, "velocity y"),
                                parse_double(vz, "velocity z"));
                    }
                }

                const int owner = decomp.owner_rank_from_frac(frac);
                LocalAtom atom(cart,
                               frac,
                               vel,
                               ModuleBase::Vector3<double>(0.0, 0.0, 0.0),
                               mbl,
                               metadata.masses[it] / ModuleBase::AU_to_MASS,
                               static_cast<int>(it),
                               ia,
                               global_id,
                               owner,
                               false);
                per_rank[static_cast<std::size_t>(owner)].push_back(pack_atom(atom));
                ++global_id;
            }
        }

        owned_atoms.clear();
        owned_atoms.reserve(per_rank[0].size());
        for (std::size_t i = 0; i < per_rank[0].size(); ++i)
        {
            owned_atoms.push_back(unpack_atom(per_rank[0][i], 0));
        }
        for (int irank = 1; irank < size; ++irank)
        {
            const int count = static_cast<int>(per_rank[static_cast<std::size_t>(irank)].size());
            MPI_Send(&count, 1, MPI_INT, irank, 8300, comm);
            if (count > 0)
            {
                MPI_Send(&per_rank[static_cast<std::size_t>(irank)][0],
                         count * static_cast<int>(sizeof(PackedMdAtom)),
                         MPI_BYTE,
                         irank,
                         8301,
                         comm);
            }
        }
    }
    else
    {
        int count = 0;
        MPI_Recv(&count, 1, MPI_INT, 0, 8300, comm, MPI_STATUS_IGNORE);
        std::vector<PackedMdAtom> buffer(static_cast<std::size_t>(count));
        if (count > 0)
        {
            MPI_Recv(&buffer[0],
                     count * static_cast<int>(sizeof(PackedMdAtom)),
                     MPI_BYTE,
                     0,
                     8301,
                     comm,
                     MPI_STATUS_IGNORE);
        }
        owned_atoms.clear();
        owned_atoms.reserve(buffer.size());
        for (std::size_t i = 0; i < buffer.size(); ++i)
        {
            owned_atoms.push_back(unpack_atom(buffer[i], rank));
        }
    }
}
} // namespace

MdCell DistributedMdCellReader::read_lj_stru(const std::string& stru_file,
                                             MPI_Comm comm,
                                             double cutoff_bohr,
                                             double skin_bohr)
{
    int rank = 0;
    MPI_Comm_rank(comm, &rank);

    RootMetadata metadata;
    metadata.lat0 = 0.0;
    metadata.omega = 0.0;
    metadata.nat = 0;
    std::vector<LocalAtom> owned_atoms;

    std::ifstream ifs;
    if (rank == 0)
    {
        ifs.open(stru_file.c_str(), std::ios::in);
        if (!ifs)
        {
            throw std::runtime_error("Failed to open STRU file: " + stru_file);
        }
        metadata = parse_root_metadata(ifs);
    }

    broadcast_metadata(metadata, comm);
    distribute_atoms_from_root(ifs, metadata, comm, cutoff_bohr, skin_bohr, owned_atoms);

    int local_max_global_id = -1;
    for (std::size_t i = 0; i < owned_atoms.size(); ++i)
    {
        local_max_global_id = std::max(local_max_global_id, static_cast<int>(owned_atoms[i].global_id));
    }
    int global_max_global_id = local_max_global_id;
    MPI_Allreduce(MPI_IN_PLACE, &global_max_global_id, 1, MPI_INT, MPI_MAX, comm);
    const int nat = global_max_global_id + 1;

    MdCell mdcell(metadata.latvec,
                  metadata.gt,
                  metadata.lat0,
                  metadata.omega,
                  owned_atoms,
                  metadata.labels,
                  metadata.masses,
                  comm,
                  cutoff_bohr,
                  skin_bohr);
    if (mdcell.nat() != nat)
    {
        throw std::runtime_error("DistributedMdCellReader global atom count mismatch.");
    }
    return mdcell;
}

#endif
