#include "source_cell/distributed_mdcell_reader.h"

#include "source_base/constants.h"
#include "source_base/vector3.h"
#include "source_cell/md_cell.h"

#ifdef __MPI
#include "source_cell/module_neighlist/domain_decomposition.h"
#endif

#include <cctype>
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
    MdStruMetadata stru_metadata;
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
            for (std::size_t it = 0; it < metadata.stru_metadata.species.size(); ++it)
            {
                metadata.stru_metadata.species[it].orbital_file = next_data_line(ifs, "NUMERICAL_ORBITAL body");
            }
            continue;
        }
        if (line == "NUMERICAL_DESCRIPTOR")
        {
            metadata.stru_metadata.descriptor_file = next_data_line(ifs, "NUMERICAL_DESCRIPTOR body");
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
        MdStruSpecies species;
        species.label = label;
        species.mass = metadata.masses.back();
        iss >> species.pseudo_file >> species.pseudo_type;
        metadata.stru_metadata.species.push_back(species);
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
                                         double cutoff_bohr,
                                         double skin_bohr,
                                         int& nat)
{
    int rank = 0;
#ifdef __MPI
    DomainDecomposition decomposition;
    decomposition.init(MPI_COMM_WORLD, metadata.latvec, metadata.lat0, cutoff_bohr, skin_bohr);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

    expect_keyword(ifs, "ATOMIC_POSITIONS");
    const std::string coord_type = next_data_line(ifs, "ATOMIC_POSITIONS type");
    const bool is_cartesian = coord_type == "Cartesian";
    const bool is_direct = coord_type == "Direct";
    if (!is_cartesian && !is_direct)
    {
        throw std::runtime_error("Only Direct and Cartesian ATOMIC_POSITIONS are supported for LJ MD.");
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
        magnetism >> metadata.stru_metadata.species[it].start_mag;
        const int nat_type = parse_int(next_data_line(ifs, "atom count"), "atom count");
        metadata.stru_metadata.species[it].atom_count = nat_type;

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

            int owner = 0;
#ifdef __MPI
            owner = decomposition.owner_rank_from_frac(frac);
#endif
            if (owner == rank)
            {
                owned_atoms.push_back(LocalAtom(cart,
                                                 frac,
                                                 vel,
                                                 ModuleBase::Vector3<double>(0.0, 0.0, 0.0),
                                                 mbl,
                                                 metadata.masses[it] / ModuleBase::AU_to_MASS,
                                                 static_cast<int>(it),
                                                 ia,
                                                 owner,
                                                 false));
            }
            ++nat;
        }
    }
    return owned_atoms;
}
} // namespace

MdCell DistributedMdCellReader::read_lj_stru(const std::string& stru_file,
                                             double cutoff_bohr,
                                             double skin_bohr)
{
    if (cutoff_bohr <= 0.0)
    {
        throw std::runtime_error("MdCell requires a positive LJ cutoff from Parameter.");
    }

    std::ifstream ifs(stru_file.c_str(), std::ios::in);
    if (!ifs)
    {
        throw std::runtime_error("Failed to open STRU file: " + stru_file);
    }

    StruMetadata metadata = parse_stru_metadata(ifs);
    int nat = 0;
    const std::vector<LocalAtom> owned_atoms = read_owned_atoms(ifs, metadata, cutoff_bohr, skin_bohr, nat);
    MdCell mdcell(metadata.latvec,
                  metadata.gt,
                  metadata.lat0,
                  metadata.omega,
                  nat,
                  owned_atoms,
                  metadata.labels,
                  metadata.masses,
                  cutoff_bohr,
                  skin_bohr);
    mdcell.set_stru_metadata(metadata.stru_metadata);
    return mdcell;
}
