#ifndef LOCAL_ATOM_H
#define LOCAL_ATOM_H

#include "source_cell/module_neighlist/neighbor_types.h"
#include "source_base/vector3.h"

#include <cstdint>

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
    ModuleBase::Vector3<double> vel;
    ModuleBase::Vector3<double> force;
    ModuleBase::Vector3<int> mbl;
    double mass;
    int type;
    std::int64_t type_index;
    int owner_rank;

    LocalAtom()
        : cart(0.0, 0.0, 0.0),
          frac(0.0, 0.0, 0.0),
          vel(0.0, 0.0, 0.0),
          force(0.0, 0.0, 0.0),
          mbl(1, 1, 1),
          mass(1.0),
          type(0),
          type_index(0),
          owner_rank(0)
    {
    }

    LocalAtom(const ModuleBase::Vector3<double>& cart_in,
              const ModuleBase::Vector3<double>& frac_in,
              const ModuleBase::Vector3<double>& vel_in,
              const ModuleBase::Vector3<double>& force_in,
              const ModuleBase::Vector3<int>& mbl_in,
              double mass_in,
              int type_in,
              std::int64_t type_index_in,
              int owner_rank_in)
        : cart(cart_in),
          frac(frac_in),
          vel(vel_in),
          force(force_in),
          mbl(mbl_in),
          mass(mass_in),
          type(type_in),
          type_index(type_index_in),
          owner_rank(owner_rank_in)
    {
    }
};

#endif // LOCAL_ATOM_H
