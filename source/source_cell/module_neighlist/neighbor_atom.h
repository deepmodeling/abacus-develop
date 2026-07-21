#ifndef NEIGHBOR_ATOM_H
#define NEIGHBOR_ATOM_H

#include "source_cell/module_neighlist/neighbor_types.h"

#include <vector>

/**
 * @brief Represents an atom with neighbor search related properties.
 *
 * This class stores atom position, type, and index information used during
 * neighbor search operations. It distinguishes between atoms inside the
 * local MPI domain and ghost atoms from neighboring domains.
 */
class NeighborAtom
{
public:
    /// X coordinate of the atom in Cartesian coordinates
    double position_x;

    /// Y coordinate of the atom in Cartesian coordinates
    double position_y;

    /// Z coordinate of the atom in Cartesian coordinates
    double position_z;

    /// Atom type index
    int atom_type;

    /// Index of the atom within its type
    int atom_index;

    /// Rank-local atom ID used by the neighbor list.
    ModuleNeighList::LocalAtomIndex atom_id;

    /// Global atom ID in the primary cell. Rank-local images share this ID.
    ModuleNeighList::GlobalAtomId global_id;

    /// MPI rank that owns the primary atom.
    int owner_rank;

    /**
     * @brief Construct a NeighborAtom.
     *
     * @param x X coordinate.
     * @param y Y coordinate.
     * @param z Z coordinate.
     * @param type Atom type index.
     * @param index Index within the atom type.
     * @param id Unique atom ID.
     */
    NeighborAtom(double x,
                 double y,
                 double z,
                 int type,
                 int index,
                 ModuleNeighList::LocalAtomIndex id)
        : position_x(x), position_y(y), position_z(z),
          atom_type(type), atom_index(index), atom_id(id),
          global_id(id), owner_rank(0) {}

    NeighborAtom(double x,
                 double y,
                 double z,
                 int type,
                 int index,
                 ModuleNeighList::LocalAtomIndex id,
                 ModuleNeighList::GlobalAtomId global_id_in,
                 int owner_rank_in)
        : position_x(x),
          position_y(y),
          position_z(z),
          atom_type(type),
          atom_index(index),
          atom_id(id),
          global_id(global_id_in),
          owner_rank(owner_rank_in)
    {
    }
};

#endif // NEIGHBOR_ATOM_H
