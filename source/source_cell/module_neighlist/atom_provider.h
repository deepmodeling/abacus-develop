#ifndef ATOM_PROVIDER_H
#define ATOM_PROVIDER_H

#include "source_base/vector3.h"
#include "source_base/matrix3.h"

/**
 * @brief Interface for providing atom and lattice information.
 *
 * This abstract interface defines the minimum set of methods needed by
 * the neighbor search module to access atom positions and lattice parameters.
 * Any class implementing this interface can be used with NeighborSearch.
 *
 * @see UnitCell
 * @see UnitCellLite
 */
class AtomProvider
{
public:
    /**
     * @brief Default destructor.
     */
    virtual ~AtomProvider() = default;

    /**
     * @brief Get the lattice constant.
     * @return Lattice constant in Bohr.
     */
    virtual double get_lat0() const = 0;

    /**
     * @brief Get the volume of the unit cell.
     * @return Unit cell volume in Bohr^3.
     */
    virtual double get_omega() const = 0;

    /**
     * @brief Get the lattice vectors.
     * @return Const reference to the 3x3 lattice vector matrix.
     */
    virtual const ModuleBase::Matrix3& get_latvec() const = 0;

    /**
     * @brief Get the total number of atoms.
     * @return Total atom count.
     */
    virtual int get_natom() const = 0;

    /**
     * @brief Get the number of atoms of a specific type.
     * @param i Type index.
     * @return Number of atoms of type i.
     */
    virtual int get_na(int i) const = 0;

    /**
     * @brief Get the number of atom types.
     * @return Number of atom types.
     */
    virtual int get_ntype() const = 0;

    /**
     * @brief Get the Cartesian coordinates of a specific atom.
     *
     * Returns the position of the j-th atom of type i.
     *
     * @param i Type index.
     * @param j Atom index within type i.
     * @return Cartesian position vector.
     */
    virtual ModuleBase::Vector3<double> get_tau(int i, int j) const = 0;
};

#endif // ATOM_PROVIDER_H