#include "source_cell/mdcell_cif_io.h"

#include "source_base/parallel_cell.h"
#include "source_cell/cif_io.h"
#include "source_cell/mdcell.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace
{
double vector_length(const ModuleBase::Vector3<double>& vector)
{
    const double squared_length = vector.x * vector.x + vector.y * vector.y + vector.z * vector.z;
    return std::sqrt(squared_length);
}

double vector_angle(const ModuleBase::Vector3<double>& left, const ModuleBase::Vector3<double>& right)
{
    const double cosine = (left.x * right.x + left.y * right.y + left.z * right.z)
                          / (vector_length(left) * vector_length(right));
    const double upper_bounded_cosine = std::min(1.0, cosine);
    const double bounded_cosine = std::max(-1.0, upper_bounded_cosine);
    const double negative_one = -1.0;
    const double radians_to_degrees = 180.0 / std::acos(negative_one);
    return std::acos(bounded_cosine) * radians_to_degrees;
}
}

namespace mdcell
{
void print_cif_file(const MDCell& cell, const std::string& filename, const std::string& frame_comment)
{
    ModuleBase::CommunicationDomain domain;
#ifdef __MPI
    domain.initialize(cell.communicator());
#endif
    if (cell.nat() <= 0 || cell.nat() > std::numeric_limits<int>::max() / 3)
    {
        throw std::runtime_error("MDCell CIF atom count is outside the supported range.");
    }
    const int nat = static_cast<int>(cell.nat());
    const int coordinate_count = 3 * nat;
    std::vector<int> owners(nat, 0);
    std::vector<double> coordinates(coordinate_count, 0.0);
    std::vector<int> offsets(cell.type_atom_counts().size() + 1, 0);
    int invalid_identity = 0;
    for (std::size_t it = 0; it < cell.type_atom_counts().size(); ++it)
    {
        if (cell.type_atom_counts()[it] < 0 || cell.type_atom_counts()[it] > nat - offsets[it])
        {
            invalid_identity = 1;
            break;
        }
        offsets[it + 1] = offsets[it] + static_cast<int>(cell.type_atom_counts()[it]);
    }
    if (offsets.back() != nat || cell.type_atom_counts().size() != cell.type_labels().size())
    {
        invalid_identity = 1;
    }
    if (domain.max(invalid_identity) != 0)
    {
        throw std::runtime_error("MDCell CIF type counts do not match the atom count.");
    }
    for (const LocalAtom& atom : cell.owned_atoms())
    {
        if (atom.type < 0 || static_cast<std::size_t>(atom.type) >= cell.type_atom_counts().size()
            || atom.type_index < 0 || atom.type_index >= cell.type_atom_counts()[atom.type])
        {
            invalid_identity = 1;
            continue;
        }
        const int atom_index = offsets[atom.type] + static_cast<int>(atom.type_index);
        ++owners[atom_index];
        coordinates[3 * atom_index] += atom.frac.x;
        coordinates[3 * atom_index + 1] += atom.frac.y;
        coordinates[3 * atom_index + 2] += atom.frac.z;
    }
    if (domain.max(invalid_identity) != 0)
    {
        throw std::runtime_error("Invalid MDCell CIF atom identity.");
    }
    domain.sum(owners.data(), nat);
    domain.sum(coordinates.data(), coordinate_count);
    int invalid_ownership = 0;
    for (const int owner_count : owners)
    {
        if (owner_count != 1)
        {
            invalid_ownership = 1;
        }
    }
    if (domain.max(invalid_ownership) != 0)
    {
        throw std::runtime_error("Duplicate or missing MDCell CIF atom index.");
    }
    if (domain.rank() != 0)
    {
        return;
    }
    const ModuleBase::Matrix3& lattice = cell.latvec();
    const ModuleBase::Vector3<double> a(lattice.e11, lattice.e12, lattice.e13);
    const ModuleBase::Vector3<double> b(lattice.e21, lattice.e22, lattice.e23);
    const ModuleBase::Vector3<double> c(lattice.e31, lattice.e32, lattice.e33);
    const double bohr_to_angstrom = 0.52917721067;
    const double length_scale = cell.lat0() * bohr_to_angstrom;
    std::vector<double> angles(6);
    angles[0] = vector_length(a) * length_scale;
    angles[1] = vector_length(b) * length_scale;
    angles[2] = vector_length(c) * length_scale;
    angles[3] = vector_angle(b, c);
    angles[4] = vector_angle(a, c);
    angles[5] = vector_angle(a, b);
    std::vector<std::string> labels(nat);
    for (std::size_t it = 0; it < cell.type_labels().size(); ++it)
    {
        for (int iat = offsets[it]; iat < offsets[it + 1]; ++iat)
        {
            labels[iat] = cell.type_labels()[it];
        }
    }
    const std::string data_tag = "data_?";
    const ModuleBase::CommunicationDomain world = ModuleBase::world_comm_domain();
    const int root_rank = world.rank();
    ModuleIO::CifParser::write(filename, angles, labels, coordinates, frame_comment, data_tag, root_rank);
}
}
