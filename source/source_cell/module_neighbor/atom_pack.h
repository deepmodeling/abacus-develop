#ifndef MODULE_NEIGHBOR_ATOM_PACK_H
#define MODULE_NEIGHBOR_ATOM_PACK_H

#include "mpi_domain.h"
#include "source_cell/unitcell.h"

#include <tuple>
#include <vector>

namespace ModuleNeighbor
{

struct AtomPack
{
    // A compact SoA container for neighbor-search input atoms.
    // All vectors are kept at the same length and index i refers to one atom record.
    // This layout is intended to replace pointer-heavy temporary records in later
    // cell-list and MPI paths.
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;
    std::vector<int> type;
    std::vector<int> natom;

    // Periodic image shift relative to the origin cell. The origin atom uses
    // (0, 0, 0); ghost atoms reuse the shift carried by MpiAtomRecord.
    std::vector<int> cell_x;
    std::vector<int> cell_y;
    std::vector<int> cell_z;

    // global_index identifies the original atom before PBC image expansion.
    // is_ghost separates atoms imported from neighboring MPI domains from local atoms.
    std::vector<int> global_index;
    std::vector<bool> is_ghost;

    void clear();
    int size() const;
    bool empty() const;
    void reserve(const int count);
    void append_atom(const double x_in,
                     const double y_in,
                     const double z_in,
                     const int type_in,
                     const int natom_in,
                     const int cell_x_in,
                     const int cell_y_in,
                     const int cell_z_in,
                     const int global_index_in,
                     const bool is_ghost_in);
    // Convert one MpiAtomRecord into the same layout used by local atoms, so local
    // and ghost atoms can be traversed by one later neighbor-search kernel.
    void append_mpi_record(const MpiAtomRecord& record, const int type_in = -1, const int natom_in = -1);
};

struct GridStorage
{
    // Flat box storage for AtomPack indices. Atoms in box b are stored in
    // atoms_in_box[box_offset[b], box_offset[b] + box_count[b]).
    std::vector<int> atoms_in_box;
    std::vector<int> box_offset;
    std::vector<int> box_count;

    int box_nx = 0;
    int box_ny = 0;
    int box_nz = 0;
    double x_min = 0.0;
    double y_min = 0.0;
    double z_min = 0.0;
    double x_max = 0.0;
    double y_max = 0.0;
    double z_max = 0.0;
    double box_edge_length = 0.0;

    void clear();
    int box_size() const;
    int get_box_id(const int bx, const int by, const int bz) const;
    int get_box_id_from_coord(const double x, const double y, const double z) const;
};

struct NeighborPair
{
    // Complete directed pair key used for correctness comparisons and stable
    // sorting. cell_* is the periodic image shift of the neighbor atom.
    int center_type = 0;
    int center_natom = 0;
    int neighbor_type = 0;
    int neighbor_natom = 0;
    int cell_x = 0;
    int cell_y = 0;
    int cell_z = 0;

    // Internal AtomPack indices used to recover coordinates when converting the
    // pair list back to AdjacentAtomInfo. They are deliberately ignored by
    // operator< and operator== so tests compare the physical neighbor relation.
    int center_index = -1;
    int neighbor_index = -1;

    std::tuple<int, int, int, int, int, int, int> key() const;
    bool operator<(const NeighborPair& rhs) const;
    bool operator==(const NeighborPair& rhs) const;
};

// Build a flat atom pack from UnitCell. When pbc is true, the image layers follow
// Grid::Check_Expand_Condition() so this helper remains comparable with Grid.
AtomPack build_atom_pack_from_unitcell(const UnitCell& ucell, const double radius_lat0, const bool pbc);

GridStorage build_grid_storage_from_atom_pack(const AtomPack& pack, const double box_edge_length);

std::vector<NeighborPair> build_neighbor_pairs_27(const AtomPack& pack,
                                                  const GridStorage& storage,
                                                  const double radius);

// Build the same directed neighbor-pair result as build_neighbor_pairs_27(), but
// traverse only one half of the box-neighbor domain and restore the opposite
// direction explicitly. This keeps the result directly comparable with the
// current 27-direction baseline.
std::vector<NeighborPair> build_neighbor_pairs_14(const AtomPack& pack,
                                                  const GridStorage& storage,
                                                  const double radius);

} // namespace ModuleNeighbor

#endif
