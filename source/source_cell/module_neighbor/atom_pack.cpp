#include "atom_pack.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

namespace ModuleNeighbor
{
namespace
{
struct ExpandLayers
{
    int x_plus;
    int x_minus;
    int y_plus;
    int y_minus;
    int z_plus;
    int z_minus;
};

ExpandLayers compute_expand_layers(const UnitCell& ucell, const double radius_lat0, const bool pbc)
{
    if (radius_lat0 < 0.0)
    {
        throw std::invalid_argument("AtomPack radius must be non-negative.");
    }
    if (!pbc)
    {
        // Non-periodic search should not manufacture image atoms. The loop bounds
        // below use [-minus, plus), so {plus=1, minus=0} keeps only shift (0,0,0).
        return ExpandLayers{1, 0, 1, 0, 1, 0};
    }

    // Use the same projection-based layer estimate as Grid::Check_Expand_Condition().
    // The three cross products measure the cell-face areas used to convert a
    // Cartesian cutoff into lattice-vector image counts for skewed cells.
    const double a23_1 = ucell.latvec.e22 * ucell.latvec.e33 - ucell.latvec.e23 * ucell.latvec.e32;
    const double a23_2 = ucell.latvec.e21 * ucell.latvec.e33 - ucell.latvec.e23 * ucell.latvec.e31;
    const double a23_3 = ucell.latvec.e21 * ucell.latvec.e32 - ucell.latvec.e22 * ucell.latvec.e31;
    const double a23_norm = std::sqrt(a23_1 * a23_1 + a23_2 * a23_2 + a23_3 * a23_3);
    const int extend_x = static_cast<int>(
        std::ceil(a23_norm * radius_lat0 / ucell.omega * ucell.lat0 * ucell.lat0 * ucell.lat0));

    const double a31_1 = ucell.latvec.e32 * ucell.latvec.e13 - ucell.latvec.e33 * ucell.latvec.e12;
    const double a31_2 = ucell.latvec.e31 * ucell.latvec.e13 - ucell.latvec.e33 * ucell.latvec.e11;
    const double a31_3 = ucell.latvec.e31 * ucell.latvec.e12 - ucell.latvec.e32 * ucell.latvec.e11;
    const double a31_norm = std::sqrt(a31_1 * a31_1 + a31_2 * a31_2 + a31_3 * a31_3);
    const int extend_y = static_cast<int>(
        std::ceil(a31_norm * radius_lat0 / ucell.omega * ucell.lat0 * ucell.lat0 * ucell.lat0));

    const double a12_1 = ucell.latvec.e12 * ucell.latvec.e23 - ucell.latvec.e13 * ucell.latvec.e22;
    const double a12_2 = ucell.latvec.e11 * ucell.latvec.e23 - ucell.latvec.e13 * ucell.latvec.e21;
    const double a12_3 = ucell.latvec.e11 * ucell.latvec.e22 - ucell.latvec.e12 * ucell.latvec.e21;
    const double a12_norm = std::sqrt(a12_1 * a12_1 + a12_2 * a12_2 + a12_3 * a12_3);
    const int extend_z = static_cast<int>(
        std::ceil(a12_norm * radius_lat0 / ucell.omega * ucell.lat0 * ucell.lat0 * ucell.lat0));

    // Grid expands images with loops over [-minus, plus). Keeping plus one layer
    // larger preserves the current asymmetric loop convention and test comparability.
    return ExpandLayers{extend_x + 1, extend_x, extend_y + 1, extend_y, extend_z + 1, extend_z};
}

int count_images(const ExpandLayers& layers)
{
    return (layers.x_plus + layers.x_minus) * (layers.y_plus + layers.y_minus)
           * (layers.z_plus + layers.z_minus);
}

int clamp_index(const int index, const int size)
{
    if (index < 0)
    {
        return 0;
    }
    if (index >= size)
    {
        return size - 1;
    }
    return index;
}
} // namespace

void AtomPack::clear()
{
    x.clear();
    y.clear();
    z.clear();
    type.clear();
    natom.clear();
    cell_x.clear();
    cell_y.clear();
    cell_z.clear();
    global_index.clear();
    is_ghost.clear();
}

int AtomPack::size() const
{
    return static_cast<int>(x.size());
}

bool AtomPack::empty() const
{
    return x.empty();
}

void AtomPack::reserve(const int count)
{
    x.reserve(count);
    y.reserve(count);
    z.reserve(count);
    type.reserve(count);
    natom.reserve(count);
    cell_x.reserve(count);
    cell_y.reserve(count);
    cell_z.reserve(count);
    global_index.reserve(count);
    is_ghost.reserve(count);
}

void AtomPack::append_atom(const double x_in,
                           const double y_in,
                           const double z_in,
                           const int type_in,
                           const int natom_in,
                           const int cell_x_in,
                           const int cell_y_in,
                           const int cell_z_in,
                           const int global_index_in,
                           const bool is_ghost_in)
{
    x.push_back(x_in);
    y.push_back(y_in);
    z.push_back(z_in);
    type.push_back(type_in);
    natom.push_back(natom_in);
    cell_x.push_back(cell_x_in);
    cell_y.push_back(cell_y_in);
    cell_z.push_back(cell_z_in);
    global_index.push_back(global_index_in);
    is_ghost.push_back(is_ghost_in);
}

void AtomPack::append_mpi_record(const MpiAtomRecord& record, const int type_in, const int natom_in)
{
    append_atom(record.x,
                record.y,
                record.z,
                type_in,
                natom_in,
                record.pbc_shift[0],
                record.pbc_shift[1],
                record.pbc_shift[2],
                record.global_index,
                record.is_ghost);
}

void GridStorage::clear()
{
    atoms_in_box.clear();
    box_offset.clear();
    box_count.clear();
    box_nx = 0;
    box_ny = 0;
    box_nz = 0;
    x_min = 0.0;
    y_min = 0.0;
    z_min = 0.0;
    x_max = 0.0;
    y_max = 0.0;
    z_max = 0.0;
    box_edge_length = 0.0;
}

int GridStorage::box_size() const
{
    return box_nx * box_ny * box_nz;
}

int GridStorage::get_box_id(const int bx, const int by, const int bz) const
{
    if (box_nx <= 0 || box_ny <= 0 || box_nz <= 0)
    {
        throw std::runtime_error("GridStorage has not been initialized.");
    }
    const int bx_safe = clamp_index(bx, box_nx);
    const int by_safe = clamp_index(by, box_ny);
    const int bz_safe = clamp_index(bz, box_nz);
    return (bx_safe * box_ny + by_safe) * box_nz + bz_safe;
}

int GridStorage::get_box_id_from_coord(const double x, const double y, const double z) const
{
    if (box_edge_length <= 0.0)
    {
        throw std::runtime_error("GridStorage box edge length must be positive.");
    }
    const int bx = static_cast<int>(std::floor((x - x_min) / box_edge_length));
    const int by = static_cast<int>(std::floor((y - y_min) / box_edge_length));
    const int bz = static_cast<int>(std::floor((z - z_min) / box_edge_length));
    return get_box_id(bx, by, bz);
}

AtomPack build_atom_pack_from_unitcell(const UnitCell& ucell, const double radius_lat0, const bool pbc)
{
    const ExpandLayers layers = compute_expand_layers(ucell, radius_lat0, pbc);

    // UnitCell::tau is already stored in lattice-constant units after check_dtau().
    // A periodic image is therefore obtained by adding integer multiples of latvec.
    ModuleBase::Vector3<double> vec1(ucell.latvec.e11, ucell.latvec.e12, ucell.latvec.e13);
    ModuleBase::Vector3<double> vec2(ucell.latvec.e21, ucell.latvec.e22, ucell.latvec.e23);
    ModuleBase::Vector3<double> vec3(ucell.latvec.e31, ucell.latvec.e32, ucell.latvec.e33);

    AtomPack pack;
    pack.reserve(ucell.nat * count_images(layers));

    int global_index_base = 0;
    for (int it = 0; it < ucell.ntype; ++it)
    {
        for (int ia = 0; ia < ucell.atoms[it].na; ++ia)
        {
            // Periodic images are search candidates only. They keep the original
            // (type, natom, global_index) identity so neighbor results can be mapped
            // back to the physical atom after distance checks.
            for (int ix = -layers.x_minus; ix < layers.x_plus; ++ix)
            {
                for (int iy = -layers.y_minus; iy < layers.y_plus; ++iy)
                {
                    for (int iz = -layers.z_minus; iz < layers.z_plus; ++iz)
                    {
                        const double x = ucell.atoms[it].tau[ia].x + vec1[0] * ix + vec2[0] * iy + vec3[0] * iz;
                        const double y = ucell.atoms[it].tau[ia].y + vec1[1] * ix + vec2[1] * iy + vec3[1] * iz;
                        const double z = ucell.atoms[it].tau[ia].z + vec1[2] * ix + vec2[2] * iy + vec3[2] * iz;
                        pack.append_atom(x, y, z, it, ia, ix, iy, iz, global_index_base + ia, false);
                    }
                }
            }
        }
        global_index_base += ucell.atoms[it].na;
    }

    return pack;
}

GridStorage build_grid_storage_from_atom_pack(const AtomPack& pack, const double box_edge_length)
{
    if (box_edge_length <= 0.0)
    {
        throw std::invalid_argument("GridStorage box edge length must be positive.");
    }
    if (pack.empty())
    {
        throw std::invalid_argument("GridStorage cannot be built from an empty AtomPack.");
    }

    GridStorage storage;
    storage.box_edge_length = box_edge_length;
    storage.x_min = storage.x_max = pack.x[0];
    storage.y_min = storage.y_max = pack.y[0];
    storage.z_min = storage.z_max = pack.z[0];

    for (int i = 1; i < pack.size(); ++i)
    {
        storage.x_min = std::min(storage.x_min, pack.x[i]);
        storage.x_max = std::max(storage.x_max, pack.x[i]);
        storage.y_min = std::min(storage.y_min, pack.y[i]);
        storage.y_max = std::max(storage.y_max, pack.y[i]);
        storage.z_min = std::min(storage.z_min, pack.z[i]);
        storage.z_max = std::max(storage.z_max, pack.z[i]);
    }

    storage.box_nx = std::max(1, static_cast<int>(std::floor((storage.x_max - storage.x_min) / box_edge_length)) + 1);
    storage.box_ny = std::max(1, static_cast<int>(std::floor((storage.y_max - storage.y_min) / box_edge_length)) + 1);
    storage.box_nz = std::max(1, static_cast<int>(std::floor((storage.z_max - storage.z_min) / box_edge_length)) + 1);

    const int nbox = storage.box_size();
    storage.box_count.assign(nbox, 0);
    // First pass: count atoms per box. The count array is later reused as the
    // stable box length table, so it must not be modified during insertion.
    for (int i = 0; i < pack.size(); ++i)
    {
        ++storage.box_count[storage.get_box_id_from_coord(pack.x[i], pack.y[i], pack.z[i])];
    }

    // Prefix sum converts counts into disjoint ranges in atoms_in_box.
    // box_offset[0] stays 0; box_offset[b] is the sum of all previous counts.
    storage.box_offset.assign(nbox, 0);
    std::partial_sum(storage.box_count.begin(), storage.box_count.end() - 1, storage.box_offset.begin() + 1);

    storage.atoms_in_box.assign(pack.size(), -1);
    std::vector<int> next_offset = storage.box_offset;
    // Second pass: write AtomPack indices into their box ranges. next_offset is
    // a scratch cursor so box_offset remains the public range table.
    for (int i = 0; i < pack.size(); ++i)
    {
        const int box_id = storage.get_box_id_from_coord(pack.x[i], pack.y[i], pack.z[i]);
        storage.atoms_in_box[next_offset[box_id]++] = i;
    }

    return storage;
}

} // namespace ModuleNeighbor
