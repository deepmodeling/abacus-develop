#include "atom_pack.h"

#include <algorithm>
#include <cmath>
#include <limits>
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
    if (!std::isfinite(radius_lat0) || radius_lat0 < 0.0)
    {
        throw std::invalid_argument("AtomPack radius must be finite and non-negative.");
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
    const long long count_x = static_cast<long long>(layers.x_plus) + layers.x_minus;
    const long long count_y = static_cast<long long>(layers.y_plus) + layers.y_minus;
    const long long count_z = static_cast<long long>(layers.z_plus) + layers.z_minus;
    const long long count = count_x * count_y * count_z;
    if (count <= 0 || count > std::numeric_limits<int>::max())
    {
        throw std::overflow_error("AtomPack periodic image count exceeds the supported range.");
    }
    return static_cast<int>(count);
}

int checked_product(const int lhs, const int rhs, const char* message)
{
    if (lhs < 0 || rhs < 0
        || (lhs != 0 && rhs > std::numeric_limits<int>::max() / lhs))
    {
        throw std::overflow_error(message);
    }
    return lhs * rhs;
}

int checked_grid_dimension(const double minimum,
                           const double maximum,
                           const double edge_length)
{
    const double dimension = std::floor((maximum - minimum) / edge_length) + 1.0;
    if (!std::isfinite(dimension) || dimension < 1.0
        || dimension > static_cast<double>(std::numeric_limits<int>::max()))
    {
        throw std::overflow_error("GridStorage dimension exceeds the supported range.");
    }
    return static_cast<int>(dimension);
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

void validate_neighbor_search_input(const GridStorage& storage, const double radius)
{
    if (!std::isfinite(radius) || radius < 0.0)
    {
        throw std::invalid_argument("Neighbor search radius must be finite and non-negative.");
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
    if (count < 0)
    {
        throw std::invalid_argument("AtomPack reserve count cannot be negative.");
    }
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
    if (x.size() >= static_cast<std::size_t>(std::numeric_limits<int>::max()))
    {
        throw std::overflow_error("AtomPack size exceeds the supported integer index range.");
    }
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
    const int atom_type = type_in >= 0 ? type_in : record.type;
    const int atom_natom = natom_in >= 0 ? natom_in : record.natom;
    append_atom(record.x,
                record.y,
                record.z,
                atom_type,
                atom_natom,
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
    return checked_product(checked_product(box_nx,
                                           box_ny,
                                           "GridStorage box count exceeds the supported range."),
                           box_nz,
                           "GridStorage box count exceeds the supported range.");
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
    if (!std::isfinite(box_edge_length) || box_edge_length <= 0.0)
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

void PagedNeighborList::clear()
{
    page_data.clear();
    page_used.clear();
    page_next.clear();
    center_first_page.clear();
    center_last_page.clear();
    center_count.clear();
    type_offset.clear();
}

bool PagedNeighborList::empty() const
{
    return total_neighbors() == 0;
}

int PagedNeighborList::page_count() const
{
    return static_cast<int>(page_used.size());
}

int PagedNeighborList::center_size() const
{
    return static_cast<int>(center_count.size());
}

int PagedNeighborList::total_neighbors() const
{
    long long count = 0;
    for (const int center_neighbors: center_count)
    {
        count += center_neighbors;
        if (count > std::numeric_limits<int>::max())
        {
            throw std::overflow_error("PagedNeighborList neighbor count exceeds the supported range.");
        }
    }
    return static_cast<int>(count);
}

int PagedNeighborList::used_slots() const
{
    long long count = 0;
    for (const int used: page_used)
    {
        count += used;
        if (count > std::numeric_limits<int>::max())
        {
            throw std::overflow_error("PagedNeighborList used slot count exceeds the supported range.");
        }
    }
    return static_cast<int>(count);
}

int PagedNeighborList::capacity_slots() const
{
    return checked_product(page_count(),
                           PAGE_SIZE,
                           "PagedNeighborList capacity exceeds the supported range.");
}

double PagedNeighborList::utilization() const
{
    return capacity_slots() == 0
               ? 1.0
               : static_cast<double>(used_slots()) / static_cast<double>(capacity_slots());
}

long long PagedNeighborList::memory_usage_bytes() const
{
    const long long integer_count = static_cast<long long>(page_data.capacity())
                                    + static_cast<long long>(page_used.capacity())
                                    + static_cast<long long>(page_next.capacity())
                                    + static_cast<long long>(center_first_page.capacity())
                                    + static_cast<long long>(center_last_page.capacity())
                                    + static_cast<long long>(center_count.capacity())
                                    + static_cast<long long>(type_offset.capacity());
    return integer_count * static_cast<long long>(sizeof(int));
}

int PagedNeighborList::allocate_page()
{
    if (page_count() == std::numeric_limits<int>::max()
        || page_data.size() > page_data.max_size() - PAGE_SIZE)
    {
        throw std::overflow_error("PagedNeighborList page count exceeds the supported range.");
    }
    const int page = page_count();
    page_data.resize(page_data.size() + PAGE_SIZE, -1);
    page_used.push_back(0);
    page_next.push_back(-1);
    return page;
}

void PagedNeighborList::append_pair_index(const int center_id, const int pair_index)
{
    if (center_id < 0 || center_id >= center_size())
    {
        throw std::out_of_range("PagedNeighborList center id is out of range.");
    }

    int page = center_last_page[center_id];
    if (page < 0 || page_used[page] == PAGE_SIZE)
    {
        const int next_page = allocate_page();
        if (page < 0)
        {
            center_first_page[center_id] = next_page;
        }
        else
        {
            page_next[page] = next_page;
        }
        center_last_page[center_id] = next_page;
        page = next_page;
    }

    page_data[page * PAGE_SIZE + page_used[page]] = pair_index;
    ++page_used[page];
    ++center_count[center_id];
}

void PagedNeighborList::build(const std::vector<NeighborPair>& pairs,
                              const std::vector<int>& atom_counts)
{
    clear();
    type_offset.assign(atom_counts.size() + 1, 0);
    for (int type = 0; type < static_cast<int>(atom_counts.size()); ++type)
    {
        if (atom_counts[type] < 0)
        {
            throw std::invalid_argument("PagedNeighborList atom count cannot be negative.");
        }
        if (atom_counts[type] > std::numeric_limits<int>::max() - type_offset[type])
        {
            throw std::overflow_error("PagedNeighborList center count exceeds the supported range.");
        }
        type_offset[type + 1] = type_offset[type] + atom_counts[type];
    }

    const int centers = type_offset.back();
    center_first_page.assign(centers, -1);
    center_last_page.assign(centers, -1);
    center_count.assign(centers, 0);
    // The exact page count depends on per-center fragmentation. This upper
    // estimate prevents page_data from repeatedly reallocating while pages are
    // appended, without permanently materializing empty pages.
    if (pairs.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
    {
        throw std::overflow_error("PagedNeighborList pair count exceeds the supported range.");
    }
    const std::size_t estimated_pages = pairs.size() / PAGE_SIZE
                                        + static_cast<std::size_t>(centers);
    if (estimated_pages > page_used.max_size()
        || estimated_pages > page_data.max_size() / PAGE_SIZE)
    {
        throw std::overflow_error("PagedNeighborList reserve size exceeds the supported range.");
    }
    page_data.reserve(estimated_pages * PAGE_SIZE);
    page_used.reserve(estimated_pages);
    page_next.reserve(estimated_pages);
    for (int pair_index = 0; pair_index < static_cast<int>(pairs.size()); ++pair_index)
    {
        const NeighborPair& pair = pairs[pair_index];
        const int center_id = get_center_id(pair.center_type, pair.center_natom);
        append_pair_index(center_id, pair_index);
    }
}

int PagedNeighborList::get_center_id(const int type, const int natom) const
{
    if (type < 0 || type + 1 >= static_cast<int>(type_offset.size()))
    {
        throw std::out_of_range("PagedNeighborList atom type is out of range.");
    }
    const int count_for_type = type_offset[type + 1] - type_offset[type];
    if (natom < 0 || natom >= count_for_type)
    {
        throw std::out_of_range("PagedNeighborList atom index is out of range.");
    }
    return type_offset[type] + natom;
}

int PagedNeighborList::count(const int type, const int natom) const
{
    return center_count[get_center_id(type, natom)];
}

void PagedNeighborList::for_each_pair_index(const int type,
                                            const int natom,
                                            const std::function<void(int)>& visitor) const
{
    const int center_id = get_center_id(type, natom);
    int page = center_first_page[center_id];
    int visited_pages = 0;
    int visited_pairs = 0;
    while (page >= 0)
    {
        if (page >= page_count() || ++visited_pages > page_count())
        {
            throw std::runtime_error("PagedNeighborList page chain is invalid.");
        }
        if (page_used[page] < 0 || page_used[page] > PAGE_SIZE)
        {
            throw std::runtime_error("PagedNeighborList page usage is invalid.");
        }
        const int page_begin = page * PAGE_SIZE;
        if (page_begin < 0
            || page_begin + page_used[page] > static_cast<int>(page_data.size()))
        {
            throw std::runtime_error("PagedNeighborList page storage is truncated.");
        }
        for (int slot = 0; slot < page_used[page]; ++slot)
        {
            const int pair_index = page_data[page_begin + slot];
            if (pair_index < 0)
            {
                throw std::runtime_error("PagedNeighborList contains an invalid pair index.");
            }
            visitor(pair_index);
            ++visited_pairs;
        }
        page = page_next[page];
    }
    if (visited_pairs != center_count[center_id])
    {
        throw std::runtime_error("PagedNeighborList center count does not match its page chain.");
    }
}

std::vector<std::array<int, 3>> build_periodic_image_shifts(const UnitCell& ucell,
                                                            const double radius_lat0,
                                                            const bool pbc)
{
    const ExpandLayers layers = compute_expand_layers(ucell, radius_lat0, pbc);
    std::vector<std::array<int, 3>> shifts;
    shifts.reserve(count_images(layers));
    for (int ix = -layers.x_minus; ix < layers.x_plus; ++ix)
    {
        for (int iy = -layers.y_minus; iy < layers.y_plus; ++iy)
        {
            for (int iz = -layers.z_minus; iz < layers.z_plus; ++iz)
            {
                shifts.push_back(std::array<int, 3>{{ix, iy, iz}});
            }
        }
    }
    return shifts;
}

AtomPack build_atom_pack_from_unitcell(const UnitCell& ucell, const double radius_lat0, const bool pbc)
{
    const std::vector<std::array<int, 3>> shifts
        = build_periodic_image_shifts(ucell, radius_lat0, pbc);

    // UnitCell::tau is already stored in lattice-constant units after check_dtau().
    // A periodic image is therefore obtained by adding integer multiples of latvec.
    ModuleBase::Vector3<double> vec1(ucell.latvec.e11, ucell.latvec.e12, ucell.latvec.e13);
    ModuleBase::Vector3<double> vec2(ucell.latvec.e21, ucell.latvec.e22, ucell.latvec.e23);
    ModuleBase::Vector3<double> vec3(ucell.latvec.e31, ucell.latvec.e32, ucell.latvec.e33);

    AtomPack pack;
    if (ucell.nat < 0
        || shifts.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())
        || (ucell.nat != 0
            && shifts.size()
                   > static_cast<std::size_t>(std::numeric_limits<int>::max() / ucell.nat)))
    {
        throw std::overflow_error("AtomPack expanded atom count exceeds the supported range.");
    }
    pack.reserve(ucell.nat * static_cast<int>(shifts.size()));

    int global_index_base = 0;
    for (int it = 0; it < ucell.ntype; ++it)
    {
        for (int ia = 0; ia < ucell.atoms[it].na; ++ia)
        {
            // Periodic images are search candidates only. They keep the original
            // (type, natom, global_index) identity so neighbor results can be mapped
            // back to the physical atom after distance checks.
            for (const std::array<int, 3>& shift: shifts)
            {
                const int ix = shift[0];
                const int iy = shift[1];
                const int iz = shift[2];
                const double x = ucell.atoms[it].tau[ia].x + vec1[0] * ix + vec2[0] * iy + vec3[0] * iz;
                const double y = ucell.atoms[it].tau[ia].y + vec1[1] * ix + vec2[1] * iy + vec3[1] * iz;
                const double z = ucell.atoms[it].tau[ia].z + vec1[2] * ix + vec2[2] * iy + vec3[2] * iz;
                pack.append_atom(x, y, z, it, ia, ix, iy, iz, global_index_base + ia, false);
            }
        }
        global_index_base += ucell.atoms[it].na;
    }

    return pack;
}

GridStorage build_grid_storage_from_atom_pack(const AtomPack& pack, const double box_edge_length)
{
    if (!std::isfinite(box_edge_length) || box_edge_length <= 0.0)
    {
        throw std::invalid_argument("GridStorage box edge length must be finite and positive.");
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

    storage.box_nx = checked_grid_dimension(storage.x_min, storage.x_max, box_edge_length);
    storage.box_ny = checked_grid_dimension(storage.y_min, storage.y_max, box_edge_length);
    storage.box_nz = checked_grid_dimension(storage.z_min, storage.z_max, box_edge_length);

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

    for (int base_box_id = 0; base_box_id < storage.box_size(); ++base_box_id)
    {
        int center_bx = 0;
        int center_by = 0;
        int center_bz = 0;
        box_id_to_indices(storage, base_box_id, center_bx, center_by, center_bz);

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

                    const int adj_box_id = storage.get_box_id(bx, by, bz);
                    const int begin = storage.box_offset[base_box_id];
                    const int end = begin + storage.box_count[base_box_id];
                    const int adj_begin = storage.box_offset[adj_box_id];
                    const int adj_end = adj_begin + storage.box_count[adj_box_id];
                    for (int offset = begin; offset < end; ++offset)
                    {
                        const int atom_a = storage.atoms_in_box[offset];
                        const int candidate_begin = adj_box_id == base_box_id ? offset + 1 : adj_begin;
                        for (int adj_offset = candidate_begin; adj_offset < adj_end; ++adj_offset)
                        {
                            const int atom_b = storage.atoms_in_box[adj_offset];
                            if (!is_within_radius(pack, atom_a, atom_b, radius2))
                            {
                                continue;
                            }

                            if (is_center_atom(pack, atom_a))
                            {
                                pairs.push_back(make_neighbor_pair(pack, atom_a, atom_b));
                            }
                            if (is_center_atom(pack, atom_b))
                            {
                                pairs.push_back(make_neighbor_pair(pack, atom_b, atom_a));
                            }
                        }
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
