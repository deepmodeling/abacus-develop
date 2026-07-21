#ifndef NEIGHBOR_SEARCH_H
#define NEIGHBOR_SEARCH_H

#include "source_cell/module_neighlist/neighbor_atom.h"
#include "source_cell/module_neighlist/bin_manager.h"
#include "source_cell/module_neighlist/neighbor_list.h"
#include "source_cell/module_neighlist/atom_provider.h"
#include "source_cell/module_neighlist/local_atom.h"

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
     */
    void init(const AtomProvider& ucell, double sr);

    /**
     * @brief Initialize from rank-local owned atoms and exchanged ghost atoms.
     *
     * This distributed entry point does not inspect a global UnitCell. The
     * caller is responsible for domain ownership and ghost exchange.
     *
     * @param owned_atoms Atoms owned by this rank and used as list centers.
     * @param ghost_atoms Cutoff halo atoms received from neighboring ranks.
     * @param sr Search radius (cutoff distance) in Bohr.
     * @param lat0 Lattice constant in Bohr.
     */
    void init_distributed(const std::vector<LocalAtom>& owned_atoms,
                          const std::vector<LocalAtom>& ghost_atoms,
                          double sr,
                          double lat0);

    /**
     * @brief Build the neighbor list for all inside atoms.
     *
     * Must be called after init(). Uses binning to efficiently find
     * all neighbors within the search radius.
     */
    void build_neighbors();


    // ========== Getter methods ==========
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

    /**
     * @brief Get the search radius.
     * @return Search radius in lattice units.
     */
    double get_search_radius() const;

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

private:
    // ========== Internal methods ==========

    double cross_product_norm(double a1, double a2, double a3,
                                          double b1, double b2, double b3);

    /**
     * @brief Check and compute expansion layer counts.
     *
     * Determines how many periodic images are needed to cover
     * the search radius in each lattice direction.
     *
     * @param ucell Unit cell providing lattice vectors.
     */
    void check_expand_condition(const AtomProvider& ucell, int& glayerX_minus, int& glayerX, int& glayerY_minus, int& glayerY, int& glayerZ_minus, int& glayerZ);

    /**
     * @brief Set member variables by generating periodic images.
     *
     * Populates all_atoms_ with atoms from the unit cell and
     * all required periodic images.
     *
     * @param ucell Unit cell providing atom positions.
     */
    void set_member_variables(const AtomProvider& ucell, int glayerX_minus, int glayerX, int glayerY_minus, int glayerY, int glayerZ_minus, int glayerZ);

    // ========== Data members ==========

    /// Search radius in lattice units
    double search_radius_ = 0.0;

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

    /// Offset added to expansion layers in positive directions
    static constexpr int positive_layer_offset = 1;

    /// Reserve factor for neighbor list capacity estimation
    static constexpr int neighbor_reserve_factor = 2;
};

#endif // NEIGHBOR_SEARCH_H
