#ifndef BIN_MANAGER_H
#define BIN_MANAGER_H

#include <vector>
#include "source_cell/module_neighlist/neighbor_atom.h"
#include "source_cell/module_neighlist/neighbor_list.h"

/**
 * @brief A single bin in the 3D binning grid for neighbor search.
 *
 * Each bin stores atoms that fall within its spatial region,
 * along with its position indices in the 3D grid.
 */
class Bin
{
public:
    /**
     * @brief Default constructor.
     */
    Bin() = default;

    /**
     * @brief Default destructor.
     */
    ~Bin() = default;

    // ========== Getter methods ==========

    /**
     * @brief Get the X index of this bin in the grid.
     * @return X index.
     */
    int get_id_x() const;

    /**
     * @brief Get the Y index of this bin in the grid.
     * @return Y index.
     */
    int get_id_y() const;

    /**
     * @brief Get the Z index of this bin in the grid.
     * @return Z index.
     */
    int get_id_z() const;

    /**
     * @brief Get the atoms stored in this bin.
     * @return Const reference to the atom vector.
     */
    const std::vector<NeighborAtom>& get_atoms() const;

    // ========== Setter methods (internal use) ==========

    /**
     * @brief Set the grid indices for this bin.
     * @param ix X index.
     * @param iy Y index.
     * @param iz Z index.
     */
    void set_id(int ix, int iy, int iz);

    /**
     * @brief Clear all atoms from this bin.
     */
    void clear_atoms();

    /**
     * @brief Add an atom to this bin.
     * @param atom The atom to add.
     */
    void add_atom(const NeighborAtom& atom);

private:
    /// X index in the 3D bin grid
    int id_x_ = 0;

    /// Y index in the 3D bin grid
    int id_y_ = 0;

    /// Z index in the 3D bin grid
    int id_z_ = 0;

    /// Atoms contained in this bin
    std::vector<NeighborAtom> atoms_;
};

/**
 * @brief Manager for binning atoms to accelerate neighbor search.
 *
 * This class implements a spatial binning strategy where atoms are
 * sorted into a 3D grid of bins. Neighbor search then only needs to
 * check adjacent bins, significantly reducing the search complexity.
 *
 * The workflow is:
 * 1. Call init_bins() to set up the bin grid based on atom positions
 * 2. Call do_binning() to assign atoms to bins
 * 3. Call build_atom_neighbors() to construct the neighbor list
 * 4. Call clear() to reset for the next search
 */
class BinManager
{
public:
    /**
     * @brief Initialize the bin grid based on atom positions and search radius.
     *
     * Computes the spatial bounds from atom positions and divides the
     * region into bins of size approximately equal to the search radius.
     *
     * @param sr Search radius in lattice units.
     * @param inside_atoms Atoms inside the local MPI domain.
     * @param ghost_atoms Ghost atoms from neighboring domains.
     */
    void init_bins(
        double sr,
        const std::vector<NeighborAtom>& inside_atoms,
        const std::vector<NeighborAtom>& ghost_atoms
    );

    /**
     * @brief Assign atoms to their corresponding bins.
     *
     * Must be called after init_bins(). Each atom is placed into the
     * bin that contains its spatial position.
     *
     * @param inside_atoms Atoms inside the local MPI domain.
     * @param ghost_atoms Ghost atoms from neighboring domains.
     */
    void do_binning(
        const std::vector<NeighborAtom>& inside_atoms,
        const std::vector<NeighborAtom>& ghost_atoms
    );

    /**
     * @brief Build neighbor list by searching adjacent bins.
     *
     * For each atom, searches its containing bin and all adjacent bins
     * to find neighbors within the search radius.
     *
     * @param neighbor_list Output neighbor list to populate.
     * @param atoms Atoms for which to build neighbors.
     */
    void build_atom_neighbors(
        NeighborList& neighbor_list,
        std::vector<NeighborAtom>& atoms
    );

    /**
     * @brief Clear all bins and reset internal state.
     */
    void clear();

    // ========== Getter methods ==========

    /**
     * @brief Get the number of bins in X direction.
     * @return Number of bins in X.
     */
    int get_nbinx() const;

    /**
     * @brief Get the number of bins in Y direction.
     * @return Number of bins in Y.
     */
    int get_nbiny() const;

    /**
     * @brief Get the number of bins in Z direction.
     * @return Number of bins in Z.
     */
    int get_nbinz() const;

    /**
     * @brief Get the total number of bins.
     * @return Total bin count (nbinx * nbiny * nbinz).
     */
    int get_total_bins() const;

    /**
     * @brief Get the number of atoms in a specific bin.
     * @param bin_index Index of the bin in the flat array.
     * @return Number of atoms in that bin.
     */
    int get_bin_atom_count(int bin_index) const;

private:
    /// Search radius in lattice units
    double sradius_ = 0.0;

    /// Minimum coordinates of the binned region
    double x_min_ = 0.0;
    double y_min_ = 0.0;
    double z_min_ = 0.0;

    /// Maximum coordinates of the binned region
    double x_max_ = 0.0;
    double y_max_ = 0.0;
    double z_max_ = 0.0;

    /// Size of each bin in each direction
    double bin_sizex_ = 0.0;
    double bin_sizey_ = 0.0;
    double bin_sizez_ = 0.0;

    /// Number of bins in each direction
    int nbinx_ = 1;
    int nbiny_ = 1;
    int nbinz_ = 1;

    /// All bins in the 3D grid (stored as flat array)
    std::vector<Bin> bins_;

    /**
     * @brief Compute the flat index for a bin from its 3D coordinates.
     * @param ix X index of the bin.
     * @param iy Y index of the bin.
     * @param iz Z index of the bin.
     * @return Flat index in the bins_ array.
     */
    int bin_index(int ix, int iy, int iz) const;
};

#endif // BIN_MANAGER_H