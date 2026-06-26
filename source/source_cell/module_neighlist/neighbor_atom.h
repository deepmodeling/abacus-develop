#ifndef NEIGHBOR_ATOM_H
#define NEIGHBOR_ATOM_H

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

    /// Unique atom ID across all domains and periodic images
    int atom_id;

    /// Whether this atom is inside the local MPI domain
    bool is_inside;

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
    NeighborAtom(double x, double y, double z, int type, int index, int id)
        : position_x(x), position_y(y), position_z(z),
          atom_type(type), atom_index(index), atom_id(id), is_inside(false) {}
};

/**
 * @brief Input structure for neighbor search initialization.
 *
 * Contains atom data and spatial bounds computed from input atoms,
 * used to initialize the binning grid.
 */
class InputAtoms
{
public:
    /// List of input atoms
    std::vector<NeighborAtom> InputAtom;

    /// Minimum X coordinate of the atom bounding box
    double x_low;

    /// Maximum X coordinate of the atom bounding box
    double x_high;

    /// Minimum Y coordinate of the atom bounding box
    double y_low;

    /// Maximum Y coordinate of the atom bounding box
    double y_high;

    /// Minimum Z coordinate of the atom bounding box
    double z_low;

    /// Maximum Z coordinate of the atom bounding box
    double z_high;

    /// Total number of atoms
    int n_atoms;

    /**
     * @brief Default constructor.
     *
     * Initializes bounds to zero and atom count to zero.
     */
    InputAtoms()
        : x_low(0), x_high(0), y_low(0), y_high(0), z_low(0), z_high(0), n_atoms(0) {}
};

#endif // NEIGHBOR_ATOM_H