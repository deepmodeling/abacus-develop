#ifndef NEIGHBOR_SEARCH_H
#define NEIGHBOR_SEARCH_H

#include "source_cell/module_neighlist/neighbor_atom.h"
#include "source_cell/module_neighlist/bin_manager.h"
#include "source_cell/module_neighlist/neighbor_list.h"
#include "source_cell/module_neighlist/atom_provider.h"

/**
 * @brief Neighbor search algorithm for building atom neighbor lists.
 *
 * This class implements a neighbor search algorithm that finds all atoms
 * within a given cutoff radius. It uses a binning strategy for efficiency
 * and supports MPI parallelization by decomposing the simulation domain.
 *
 * The workflow is:
 * 1. Call init() to initialize with unit cell and search radius
 * 2. Call build_neighbors() to construct the neighbor list
 * 3. Access results via get_neighbor_list()
 */
class NeighborSearch
{
public:
    /**
     * @brief Default constructor.
     */
    NeighborSearch() = default;

    /**
     * @brief Default destructor.
     */
    ~NeighborSearch() = default;

    // ========== Main public interface ==========

    /**
     * @brief Initialize the neighbor search with unit cell and search radius.
     *
     * This method sets up the domain decomposition, identifies inside and
     * ghost atoms, and prepares internal data structures.
     *
     * @param ucell Unit cell providing atom positions and lattice info.
     * @param sr Search radius (cutoff distance) in Bohr.
     * @param mpi_rank MPI rank of this process.
     */
    void init(const AtomProvider& ucell, double sr, int mpi_rank);

    /**
     * @brief Build the neighbor list for all inside atoms.
     *
     * Must be called after init(). Uses binning to efficiently find
     * all neighbors within the search radius.
     */
    void build_neighbors();

    /**
     * @brief Get the constructed neighbor list.
     * @return Reference to the NeighborList object.
     */
    NeighborList& get_neighbor_list();

    /**
     * @brief Get the constructed neighbor list (const version).
     * @return Const reference to the NeighborList object.
     */
    const NeighborList& get_neighbor_list() const;

    // ========== Utility methods (public for testing) ==========

    /**
     * @brief Calculate squared distance from a point to the local domain box.
     *
     * Used to determine if an atom is within the search radius of the
     * local MPI domain.
     *
     * @param position_x X coordinate of the point.
     * @param position_y Y coordinate of the point.
     * @param position_z Z coordinate of the point.
     * @param x_low Lower bound of the global domain in X.
     * @param y_low Lower bound of the global domain in Y.
     * @param z_low Lower bound of the global domain in Z.
     * @return Squared distance to the domain box.
     */
    double distance(double position_x,
                    double position_y,
                    double position_z,
                    double x_low,
                    double y_low,
                    double z_low);

    /**
     * @brief Decompose MPI size into a 3D grid.
     *
     * Finds a balanced decomposition of mpi_size into nx * ny * nz.
     *
     * @param mpi_size Total number of MPI processes.
     * @param nx Output: number of divisions in X.
     * @param ny Output: number of divisions in Y.
     * @param nz Output: number of divisions in Z.
     */
    void decompose(int mpi_size, int& nx, int& ny, int& nz);

    // ========== Getter methods ==========

    /**
     * @brief Get the search radius.
     * @return Search radius in lattice units.
     */
    double get_search_radius() const;

    /**
     * @brief Get the X position of this MPI domain.
     * @return Domain index in X.
     */
    int get_x() const;

    /**
     * @brief Get the Y position of this MPI domain.
     * @return Domain index in Y.
     */
    int get_y() const;

    /**
     * @brief Get the Z position of this MPI domain.
     * @return Domain index in Z.
     */
    int get_z() const;

    /**
     * @brief Get the width of this MPI domain in X.
     * @return Domain width in X.
     */
    double get_wide_x() const;

    /**
     * @brief Get the width of this MPI domain in Y.
     * @return Domain width in Y.
     */
    double get_wide_y() const;

    /**
     * @brief Get the width of this MPI domain in Z.
     * @return Domain width in Z.
     */
    double get_wide_z() const;

    /**
     * @brief Get the number of expansion layers in +X direction.
     * @return Number of layers.
     */
    int get_glayerX() const;

    /**
     * @brief Get the number of expansion layers in +Y direction.
     * @return Number of layers.
     */
    int get_glayerY() const;

    /**
     * @brief Get the number of expansion layers in +Z direction.
     * @return Number of layers.
     */
    int get_glayerZ() const;

    /**
     * @brief Get the number of expansion layers in -X direction.
     * @return Number of layers.
     */
    int get_glayerX_minus() const;

    /**
     * @brief Get the number of expansion layers in -Y direction.
     * @return Number of layers.
     */
    int get_glayerY_minus() const;

    /**
     * @brief Get the number of expansion layers in -Z direction.
     * @return Number of layers.
     */
    int get_glayerZ_minus() const;

    /**
     * @brief Get all atoms (including periodic images).
     * @return Const reference to the vector of all atoms.
     */
    const std::vector<NeighborAtom>& get_all_atoms() const;

    /**
     * @brief Get atoms inside the local MPI domain.
     * @return Const reference to the vector of inside atoms.
     */
    const std::vector<NeighborAtom>& get_inside_atoms() const;

    /**
     * @brief Get ghost atoms (neighbors of inside atoms).
     * @return Const reference to the vector of ghost atoms.
     */
    const std::vector<NeighborAtom>& get_ghost_atoms() const;

    // ========== Setter methods ==========

    /**
     * @brief Set the search radius.
     * @param sr Search radius in lattice units.
     */
    void set_search_radius(double sr);

    /**
     * @brief Set the position of this MPI domain.
     * @param x Domain index in X.
     * @param y Domain index in Y.
     * @param z Domain index in Z.
     */
    void set_position(int x, int y, int z);

    /**
     * @brief Set the width of this MPI domain.
     * @param wx Domain width in X.
     * @param wy Domain width in Y.
     * @param wz Domain width in Z.
     */
    void set_width(double wx, double wy, double wz);

private:
    // ========== Internal methods ==========

    /**
     * @brief Convert unit cell atoms to InputAtoms format.
     * @param ucell Unit cell providing atom info.
     * @return InputAtoms structure for processing.
     */
    InputAtoms ucell_to_input_atoms(const AtomProvider& ucell);

    /**
     * @brief Check and compute expansion layer counts.
     *
     * Determines how many periodic images are needed to cover
     * the search radius in each lattice direction.
     *
     * @param ucell Unit cell providing lattice vectors.
     */
    void check_expand_condition(const AtomProvider& ucell);

    /**
     * @brief Set member variables by generating periodic images.
     *
     * Populates all_atoms_ with atoms from the unit cell and
     * all required periodic images.
     *
     * @param ucell Unit cell providing atom positions.
     */
    void set_member_variables(const AtomProvider& ucell);

    /**
     * @brief Compute the norm of the cross product of two 3D vectors.
     *
     * @param a1, a2, a3 Components of the first vector.
     * @param b1, b2, b3 Components of the second vector.
     * @return Norm of the cross product.
     */
    static double cross_product_norm(double a1, double a2, double a3,
                                     double b1, double b2, double b3);

    // ========== Data members ==========

    /// Search radius in lattice units
    double search_radius_ = 0.0;

    /// Position of this MPI domain in the 3D grid
    int x_ = 0;
    int y_ = 0;
    int z_ = 0;

    /// Width of this MPI domain
    double wide_x_ = 0.0;
    double wide_y_ = 0.0;
    double wide_z_ = 0.0;

    /// Number of expansion layers in positive directions
    int glayerX_ = 0;
    int glayerY_ = 0;
    int glayerZ_ = 0;

    /// Number of expansion layers in negative directions
    int glayerX_minus_ = 0;
    int glayerY_minus_ = 0;
    int glayerZ_minus_ = 0;

    /// All atoms including periodic images
    std::vector<NeighborAtom> all_atoms_;

    /// Atoms inside the local MPI domain
    std::vector<NeighborAtom> inside_atoms_;

    /// Ghost atoms (neighbors from other domains or images)
    std::vector<NeighborAtom> ghost_atoms_;

    /// The constructed neighbor list
    NeighborList neighbor_list_;

    /// Bin manager for efficient neighbor search
    BinManager bin_manager_;
    // ========== Compile-time constants ==========

    /// Tolerance for coordinate comparisons in lattice units
    static constexpr double coord_tolerance = 1e-8;

    /// Offset added to expansion layers in positive directions
    static constexpr int positive_layer_offset = 1;

    /// Reserve factor for neighbor list capacity estimation
    static constexpr int neighbor_reserve_factor = 2;
};

#endif // NEIGHBOR_SEARCH_H