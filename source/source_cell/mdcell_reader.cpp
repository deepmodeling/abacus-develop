#include "source_cell/mdcell_reader.h"

#include "source_base/constants.h"
#include "source_base/parallel_cell.h"
#include "source_base/vector3.h"
#include "source_cell/mdcell.h"

#ifdef __MPI
#include "source_cell/module_neighlist/domain_decomposition.h"
#endif

#include <cctype>
#include <cstdint>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace
{
struct StruMetadata
{
    double lat0;
    double omega;
    ModuleBase::Matrix3 latvec;
    ModuleBase::Matrix3 gt;
    std::vector<std::string> labels;
    std::vector<double> masses;
    std::vector<std::int64_t> type_atom_counts;
    MdStruFileMetadata stru_file_metadata;
};

std::string next_data_line(std::ifstream& ifs, const char* context)
{
    std::string line;
    while (std::getline(ifs, line))
    {
        const std::size_t comment = line.find('#');
        if (comment != std::string::npos)
        {
            line.erase(comment);
        }
        std::size_t begin = 0;
        while (begin < line.size() && std::isspace(static_cast<unsigned char>(line[begin])))
        {
            ++begin;
        }
        std::size_t end = line.size();
        while (end > begin && std::isspace(static_cast<unsigned char>(line[end - 1])))
        {
            --end;
        }
        if (begin != end) return line.substr(begin, end - begin);
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

std::int64_t parse_int64(const std::string& token, const char* context)
{
    char* end = NULL;
    const long long value = std::strtoll(token.c_str(), &end, 10);
    if (end == token.c_str() || *end != '\0')
    {
        throw std::runtime_error(std::string("Failed to parse int64 for ") + context + ": " + token);
    }
    return static_cast<std::int64_t>(value);
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

StruMetadata parse_stru_metadata(std::ifstream& ifs)
{
    StruMetadata metadata;
    metadata.lat0 = 1.0;
    metadata.omega = 0.0;

    expect_keyword(ifs, "ATOMIC_SPECIES");
    while (true)
    {
        const std::streampos mark = ifs.tellg();
        const std::string line = next_data_line(ifs, "ATOMIC_SPECIES body");
        if (line == "LATTICE_CONSTANT")
        {
            ifs.seekg(mark);
            break;
        }
        if (line == "NUMERICAL_ORBITAL")
        {
            for (std::size_t it = 0; it < metadata.stru_file_metadata.species.size(); ++it)
            {
                metadata.stru_file_metadata.species[it].orbital_file = next_data_line(ifs, "NUMERICAL_ORBITAL body");
            }
            continue;
        }
        if (line == "NUMERICAL_DESCRIPTOR")
        {
            metadata.stru_file_metadata.descriptor_file = next_data_line(ifs, "NUMERICAL_DESCRIPTOR body");
            continue;
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
        MdStruFileSpecies species;
        iss >> species.pseudo_file >> species.pseudo_type;
        metadata.stru_file_metadata.species.push_back(species);
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

std::vector<LocalAtom> read_owned_atoms(std::ifstream& ifs,
                                         StruMetadata& metadata,
                                         const ModuleBase::Matrix3& primitive_latvec,
                                         const ModuleBase::Matrix3& primitive_gt,
                                         const std::vector<int>& cell_replica,
                                         double cutoff,
                                         double skin,
                                         std::int64_t& nat,
                                         const ModuleBase::CommunicationDomain& comm_domain)
{
    int rank = 0;
#ifdef __MPI
    DomainDecomposition decomposition;
    decomposition.init(comm_domain.communicator(), metadata.latvec, metadata.lat0, cutoff, skin);
    rank = comm_domain.rank();
#endif

    int begin[3] = {0, 0, 0};
    int end[3] = {cell_replica[0], cell_replica[1], cell_replica[2]};
#ifdef __MPI
    const std::array<int, 3>& dims = decomposition.dims();
    const std::array<int, 3>& coords = decomposition.coords();
    for (int idim = 0; idim < 3; ++idim)
    {
        begin[idim] = std::max(0, static_cast<int>(std::floor(
            static_cast<double>(coords[idim]) * cell_replica[idim] / dims[idim])) - 1);
        end[idim] = std::min(cell_replica[idim], static_cast<int>(std::ceil(
            static_cast<double>(coords[idim] + 1) * cell_replica[idim] / dims[idim])) + 1);
    }
#endif

    expect_keyword(ifs, "ATOMIC_POSITIONS");
    const std::string coord_type = next_data_line(ifs, "ATOMIC_POSITIONS type");
    const bool is_cartesian = coord_type == "Cartesian";
    const bool is_direct = coord_type == "Direct";
    if (!is_cartesian && !is_direct)
    {
        throw std::runtime_error("Only Direct and Cartesian ATOMIC_POSITIONS are supported for MD.");
    }

    std::vector<LocalAtom> owned_atoms;
    nat = 0;
    for (std::size_t it = 0; it < metadata.labels.size(); ++it)
    {
        const std::string label = next_data_line(ifs, "atom label");
        if (label != metadata.labels[it])
        {
            throw std::runtime_error("ATOMIC_POSITIONS label order does not match ATOMIC_SPECIES.");
        }
        std::istringstream magnetism(next_data_line(ifs, "magnetism"));
        magnetism >> metadata.stru_file_metadata.species[it].start_mag;
        const std::int64_t nat_type = parse_int64(next_data_line(ifs, "atom count"), "atom count");

        for (std::int64_t ia = 0; ia < nat_type; ++ia)
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
            if (is_cartesian)
            {
                frac = wrap_fractional(ModuleBase::Vector3<double>(c1, c2, c3) * primitive_gt);
            }
            else
            {
                frac = wrap_fractional(ModuleBase::Vector3<double>(c1, c2, c3));
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
                    if (!iss) throw std::runtime_error("Invalid move flag record in STRU.");
                    mbl.set(parse_int(mx, "move flag x"), parse_int(my, "move flag y"), parse_int(mz, "move flag z"));
                }
                else if (token == "v" || token == "vel" || token == "velocity")
                {
                    std::string vx;
                    std::string vy;
                    std::string vz;
                    iss >> vx >> vy >> vz;
                    if (!iss) throw std::runtime_error("Invalid velocity record in STRU.");
                    vel.set(parse_double(vx, "velocity x"),
                            parse_double(vy, "velocity y"),
                            parse_double(vz, "velocity z"));
                }
            }

            for (int ix = begin[0]; ix < end[0]; ++ix)
            {
                for (int iy = begin[1]; iy < end[1]; ++iy)
                {
                    for (int iz = begin[2]; iz < end[2]; ++iz)
                    {
                        ModuleBase::Vector3<double> final_frac(
                            (ix + frac.x) / cell_replica[0],
                            (iy + frac.y) / cell_replica[1],
                            (iz + frac.z) / cell_replica[2]);
                        int owner = 0;
#ifdef __MPI
                        owner = decomposition.owner_rank_from_frac(final_frac);
#endif
                        if (owner == rank)
                        {
                            owned_atoms.push_back(LocalAtom(final_frac * metadata.latvec,
                                                             final_frac,
                                                             vel,
                                                             ModuleBase::Vector3<double>(0.0, 0.0, 0.0),
                                                             mbl,
                                                             metadata.masses[it] / ModuleBase::AU_to_MASS,
                                                             static_cast<int>(it),
                                                             ((static_cast<std::int64_t>(ix) * cell_replica[1] + iy)
                                                                  * cell_replica[2] + iz) * nat_type + ia,
                                                             owner));
                        }
                    }
                }
            }
        }
        const std::int64_t replicated_atom_count = nat_type * cell_replica[0] * cell_replica[1] * cell_replica[2];
        metadata.type_atom_counts.push_back(replicated_atom_count);
        nat += replicated_atom_count;
    }
    return owned_atoms;
}
} // namespace

MDCell MDCellReader::read_stru(const std::string& stru_file,
                                          const std::vector<int>& cell_replica,
                                          double cutoff,
                                          double skin,
                                          const ModuleBase::CommunicationDomain& comm_domain)
{
    if (cutoff <= 0.0)
    {
        throw std::runtime_error("MDCell requires a positive cutoff.");
    }

    std::ifstream ifs(stru_file.c_str(), std::ios::in);
    if (!ifs)
    {
        throw std::runtime_error("Failed to open STRU file: " + stru_file);
    }

    if (cell_replica.size() != 3 || cell_replica[0] <= 0 || cell_replica[1] <= 0 || cell_replica[2] <= 0)
    {
        throw std::runtime_error("cell_replica requires three positive integers.");
    }
    StruMetadata metadata = parse_stru_metadata(ifs);
    const ModuleBase::Matrix3 primitive_latvec = metadata.latvec;
    const ModuleBase::Matrix3 primitive_gt = metadata.gt;
    metadata.latvec.e11 *= cell_replica[0]; metadata.latvec.e12 *= cell_replica[0]; metadata.latvec.e13 *= cell_replica[0];
    metadata.latvec.e21 *= cell_replica[1]; metadata.latvec.e22 *= cell_replica[1]; metadata.latvec.e23 *= cell_replica[1];
    metadata.latvec.e31 *= cell_replica[2]; metadata.latvec.e32 *= cell_replica[2]; metadata.latvec.e33 *= cell_replica[2];
    metadata.gt = metadata.latvec.Inverse();
    metadata.omega = std::abs(metadata.latvec.Det()) * metadata.lat0 * metadata.lat0 * metadata.lat0;
    std::int64_t nat = 0;
    const std::vector<LocalAtom> owned_atoms = read_owned_atoms(ifs, metadata, primitive_latvec, primitive_gt,
                                                                  cell_replica, cutoff, skin, nat, comm_domain);
    MDCell mdcell;
    mdcell.initialize_from_owned_atoms(metadata.latvec,
                                       metadata.gt,
                                       metadata.lat0,
                                       metadata.omega,
                                       nat,
                                       owned_atoms,
                                       metadata.labels,
                                       metadata.masses,
                                       metadata.type_atom_counts,
                                       cutoff,
                                       skin,
                                       comm_domain);
    mdcell.mutable_stru_file_metadata() = metadata.stru_file_metadata;
    return mdcell;
}
