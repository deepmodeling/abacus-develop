#ifndef LOCAL_ATOM_H
#define LOCAL_ATOM_H

#include "source_cell/module_neighlist/neighbor_types.h"
#include "source_base/vector3.h"

/**
 * @brief Atom record owned by a distributed neighbor-search rank.
 *
 * cart is in lattice-coordinate units, matching UnitCell::tau and the existing
 * NeighborSearch implementation. frac is wrapped into [0, 1) for owned atoms.
 * Ghost atoms may have shifted cartesian coordinates while retaining the
 * original wrapped frac coordinate for ownership metadata.
 */
struct LocalAtom
{
    ModuleBase::Vector3<double> cart;
    ModuleBase::Vector3<double> frac;
    int type;
    int type_index;
    ModuleNeighList::GlobalAtomId global_id;
    int owner_rank;
    bool is_ghost;

    LocalAtom()
        : cart(0.0, 0.0, 0.0),
          frac(0.0, 0.0, 0.0),
          type(0),
          type_index(0),
          global_id(-1),
          owner_rank(0),
          is_ghost(false)
    {
    }

    LocalAtom(const ModuleBase::Vector3<double>& cart_in,
              const ModuleBase::Vector3<double>& frac_in,
              int type_in,
              int type_index_in,
              ModuleNeighList::GlobalAtomId global_id_in,
              int owner_rank_in,
              bool is_ghost_in)
        : cart(cart_in),
          frac(frac_in),
          type(type_in),
          type_index(type_index_in),
          global_id(global_id_in),
          owner_rank(owner_rank_in),
          is_ghost(is_ghost_in)
    {
    }
};

#endif // LOCAL_ATOM_H
