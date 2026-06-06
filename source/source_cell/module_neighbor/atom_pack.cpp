#include "atom_pack.h"

#include <algorithm>
#include <cmath>
#include <map>
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

void box_id_to_indices(const GridStorage& storage, const int box_id, int& bx, int& by, int& bz)
{
    const int yz_size = storage.box_ny * storage.box_nz;
    bx = box_id / yz_size;
    const int yz_id = box_id % yz_size;
    by = yz_id / storage.box_nz;
    bz = yz_id % storage.box_nz;
}

bool is_center_atom(const AtomPack& pack, const int index)
{
    return !pack.is_ghost[index] && pack.cell_x[index] == 0 && pack.cell_y[index] == 0 && pack.cell_z[index] == 0;
}

using AtomImageKey = std::tuple<int, int, int, int, int>;

AtomImageKey make_atom_image_key(const AtomPack& pack, const int index)
{
    return std::make_tuple(pack.type[index], pack.natom[index], pack.cell_x[index], pack.cell_y[index], pack.cell_z[index]);
}

std::map<AtomImageKey, int> build_atom_image_index(const AtomPack& pack)
{
    std::map<AtomImageKey, int> index_by_image;
    for (int i = 0; i < pack.size(); ++i)
    {
        if (!pack.is_ghost[i])
        {
            index_by_image[make_atom_image_key(pack, i)] = i;
        }
    }
    return index_by_image;
}

int find_atom_image_index(const std::map<AtomImageKey, int>& index_by_image,
                          const int type,
                          const int natom,
                          const int cell_x,
                          const int cell_y,
                          const int cell_z)
{
    const auto iter = index_by_image.find(std::make_tuple(type, natom, cell_x, cell_y, cell_z));
    return iter == index_by_image.end() ? -1 : iter->second;
}

void validate_neighbor_search_input(const GridStorage& storage, const double radius)
{
    if (radius < 0.0)
    {
        throw std::invalid_argument("Neighbor search radius must be non-negative.");
    }
    if (storage.box_edge_length <= 0.0 || storage.box_size() <= 0)
    {
        throw std::runtime_error("GridStorage has not been initialized.");
    }
}

NeighborPair make_neighbor_pair(const AtomPack& pack, const int center, const int candidate)
{
    // Keep explicit assignment for C++11 CI builds: NeighborPair has default
    // member initializers and member functions, so brace-list aggregate
    // initialization is not accepted by all toolchain variants.
    NeighborPair pair;
    pair.center_type = pack.type[center];
    pair.center_natom = pack.natom[center];
    pair.neighbor_type = pack.type[candidate];
    pair.neighbor_natom = pack.natom[candidate];
    pair.cell_x = pack.cell_x[candidate];
    pair.cell_y = pack.cell_y[candidate];
    pair.cell_z = pack.cell_z[candidate];
    pair.center_index = center;
    pair.neighbor_index = candidate;
    return pair;
}

NeighborPair make_restored_reverse_pair(const AtomPack& pack,
                                        const NeighborPair& pair,
                                        const std::map<AtomImageKey, int>& index_by_image)
{
    // Half-domain search visits only one direction. The restored pair uses the
    // neighbor's origin-cell image as the new center and the opposite image shift
    // for the original center atom.
    const int reverse_center = find_atom_image_index(index_by_image, pair.neighbor_type, pair.neighbor_natom, 0, 0, 0);
    const int reverse_neighbor = find_atom_image_index(index_by_image,
                                                       pair.center_type,
                                                       pair.center_natom,
                                                       -pair.cell_x,
                                                       -pair.cell_y,
                                                       -pair.cell_z);

    NeighborPair reverse;
    reverse.center_type = pair.neighbor_type;
    reverse.center_natom = pair.neighbor_natom;
    reverse.neighbor_type = pair.center_type;
    reverse.neighbor_natom = pair.center_natom;
    reverse.cell_x = -pair.cell_x;
    reverse.cell_y = -pair.cell_y;
    reverse.cell_z = -pair.cell_z;
    reverse.center_index = reverse_center;
    reverse.neighbor_index = reverse_neighbor;
    return reverse;
}

bool is_within_radius(const AtomPack& pack, const int center, const int candidate, const double radius2)
{
    const double dx = pack.x[center] - pack.x[candidate];
    const double dy = pack.y[center] - pack.y[candidate];
    const double dz = pack.z[center] - pack.z[candidate];
    const double dr = dx * dx + dy * dy + dz * dz;
    return dr != 0.0 && dr <= radius2;
}

bool is_half_domain_offset(const int dbx, const int dby, const int dbz)
{
    return dbx > 0 || (dbx == 0 && dby > 0) || (dbx == 0 && dby == 0 && dbz >= 0);
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

std::tuple<int, int, int, int, int, int, int> NeighborPair::key() const
{
    return std::make_tuple(center_type, center_natom, neighbor_type, neighbor_natom, cell_x, cell_y, cell_z);
}

bool NeighborPair::operator<(const NeighborPair& rhs) const
{
    return key() < rhs.key();
}

bool NeighborPair::operator==(const NeighborPair& rhs) const
{
    return key() == rhs.key();
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

std::vector<NeighborPair> build_neighbor_pairs_27(const AtomPack& pack,
                                                  const GridStorage& storage,
                                                  const double radius)
{
    validate_neighbor_search_input(storage, radius);

    std::vector<NeighborPair> pairs;
    const double radius2 = radius * radius;
    const int search_layer = std::max(1, static_cast<int>(std::ceil(radius / storage.box_edge_length)));

    for (int center = 0; center < pack.size(); ++center)
    {
        if (!is_center_atom(pack, center))
        {
            continue;
        }

        int center_bx = 0;
        int center_by = 0;
        int center_bz = 0;
        box_id_to_indices(storage,
                          storage.get_box_id_from_coord(pack.x[center], pack.y[center], pack.z[center]),
                          center_bx,
                          center_by,
                          center_bz);

        const int bx_begin = std::max(0, center_bx - search_layer);
        const int bx_end = std::min(storage.box_nx - 1, center_bx + search_layer);
        const int by_begin = std::max(0, center_by - search_layer);
        const int by_end = std::min(storage.box_ny - 1, center_by + search_layer);
        const int bz_begin = std::max(0, center_bz - search_layer);
        const int bz_end = std::min(storage.box_nz - 1, center_bz + search_layer);

        for (int bx = bx_begin; bx <= bx_end; ++bx)
        {
            for (int by = by_begin; by <= by_end; ++by)
            {
                for (int bz = bz_begin; bz <= bz_end; ++bz)
                {
                    const int box_id = storage.get_box_id(bx, by, bz);
                    const int begin = storage.box_offset[box_id];
                    const int end = begin + storage.box_count[box_id];
                    for (int offset = begin; offset < end; ++offset)
                    {
                        const int candidate = storage.atoms_in_box[offset];
                        if (is_within_radius(pack, center, candidate, radius2))
                        {
                            pairs.push_back(make_neighbor_pair(pack, center, candidate));
                        }
                    }
                }
            }
        }
    }

    std::sort(pairs.begin(), pairs.end());
    return pairs;
}

std::vector<NeighborPair> build_neighbor_pairs_14(const AtomPack& pack,
                                                  const GridStorage& storage,
                                                  const double radius)
{
    validate_neighbor_search_input(storage, radius);

    std::vector<NeighborPair> pairs;
    const double radius2 = radius * radius;
    const int search_layer = std::max(1, static_cast<int>(std::ceil(radius / storage.box_edge_length)));
    const std::map<AtomImageKey, int> index_by_image = build_atom_image_index(pack);

    for (int center = 0; center < pack.size(); ++center)
    {
        if (!is_center_atom(pack, center))
        {
            continue;
        }

        int center_bx = 0;
        int center_by = 0;
        int center_bz = 0;
        box_id_to_indices(storage,
                          storage.get_box_id_from_coord(pack.x[center], pack.y[center], pack.z[center]),
                          center_bx,
                          center_by,
                          center_bz);

        for (int dbx = -search_layer; dbx <= search_layer; ++dbx)
        {
            for (int dby = -search_layer; dby <= search_layer; ++dby)
            {
                for (int dbz = -search_layer; dbz <= search_layer; ++dbz)
                {
                    if (!is_half_domain_offset(dbx, dby, dbz))
                    {
                        continue;
                    }

                    const int bx = center_bx + dbx;
                    const int by = center_by + dby;
                    const int bz = center_bz + dbz;
                    if (bx < 0 || bx >= storage.box_nx || by < 0 || by >= storage.box_ny || bz < 0
                        || bz >= storage.box_nz)
                    {
                        continue;
                    }

                    const int box_id = storage.get_box_id(bx, by, bz);
                    const int begin = storage.box_offset[box_id];
                    const int end = begin + storage.box_count[box_id];
                    for (int offset = begin; offset < end; ++offset)
                    {
                        const int candidate = storage.atoms_in_box[offset];
                        if (!is_within_radius(pack, center, candidate, radius2))
                        {
                            continue;
                        }

                        const NeighborPair forward = make_neighbor_pair(pack, center, candidate);
                        const NeighborPair reverse = make_restored_reverse_pair(pack, forward, index_by_image);
                        if (forward.key() == reverse.key())
                        {
                            continue;
                        }

                        pairs.push_back(forward);
                        pairs.push_back(reverse);
                    }
                }
            }
        }
    }

    std::sort(pairs.begin(), pairs.end());
    pairs.erase(std::unique(pairs.begin(), pairs.end()), pairs.end());
    return pairs;
}

} // namespace ModuleNeighbor
