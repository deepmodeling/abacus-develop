/**
 * @file qlist.h
 * @brief QList class for managing q-points.
 * @author Mohan Chen (added on 2026-05-18)
 * @note The q-point mesh generation and star (IBZ) reduction share the
 *       spin-free functionality of ModuleCell::ReciprocalGrid with
 *       K_Vectors; the irreducible-representation analysis is added on top.
 */
#ifndef QLIST_H
#define QLIST_H

#include "source_base/vector3.h"
#include "module_symmetry/little_group.h"
#include "module_symmetry/symmetry.h"
#include "unitcell.h"
#include "reciprocal_grid.h"
#include <vector>

namespace ModuleCell {

/**
 * @brief QList class for managing q-points.
 *
 * Inherits the spin-free reciprocal-grid functionality (mesh generation,
 * coordinate conversion, weights, star reduction primitive) from
 * ReciprocalGrid and adds the q-specific star reduction (always including
 * the time-reversal partner -q, no spin expansion) together with the
 * little-group irreducible-representation data (placeholder in the current
 * version).
 */
class QList : public ModuleCell::ReciprocalGrid {
public:
    /**
     * @brief Default constructor.
     */
    QList();

    /**
     * @brief Destructor.
     */
    ~QList();

    /**
     * @brief Generate the Monkhorst-Pack q-point mesh and reduce it by star.
     *
     * @param ucell unit cell
     * @param symm symmetry object
     * @param mp_grid Monkhorst-Pack grid
     * @param use_irreps whether to use irreps
     */
    void generate_mesh(UnitCell& ucell, ModuleSymmetry::Symmetry& symm,
                       const std::vector<int>& mp_grid, bool use_irreps);

    /**
     * @brief Read q-points from file.
     *
     * @param filename filename
     * @param ucell unit cell
     */
    void read_from_file(const std::string& filename, UnitCell& ucell);

    /**
     * @brief Get the number of q-points.
     * @return number of q-points
     */
    int get_nq() const { return this->nkstot; }

    /**
     * @brief Get q-point at given index.
     * @param idx q-point index
     * @return q-point vector (direct coordinates)
     */
    ModuleBase::Vector3<double> get_q(int idx) const { return this->kvec_d[idx]; }

    /**
     * @brief Get the number of irreps at given q-point.
     * @param idx q-point index
     * @return number of irreps
     */
    int get_nirr(int idx) const { return nirr_[idx]; }

    /**
     * @brief Get irrep modes at given q-point and irrep index.
     * @param q_idx q-point index
     * @param irrep_idx irrep index
     * @return irrep modes
     */
    std::vector<int> get_irrep_modes(int q_idx, int irrep_idx) const;

    /**
     * @brief Reduce the q-points by star (time-reversal included).
     *
     * Implements the pure-virtual hook of ReciprocalGrid: builds the
     * reciprocal-space point-group operations via build_star_ops, always
     * doubles them by -q, and folds the mesh with ReciprocalGrid::reduce_ibz.
     *
     * @param ucell unit cell
     * @param symm symmetry object
     * @param use_symm whether symmetry reduction is enabled
     * @param skpt output string (unused for q-points)
     * @param match set to false if the reciprocal lattice is not compatible
     *              with the real-space lattice
     */
    void reduce_by_symmetry(const UnitCell& ucell,
                            const ModuleSymmetry::Symmetry& symm,
                            bool use_symm,
                            std::string& skpt,
                            bool& match) override;

private:
    std::vector<int> nirr_; ///< number of irreps for each q-point
    std::vector<std::vector<std::vector<int>>> irrep_modes_; ///< irrep modes
    ModuleSymmetry::LittleGroup little_group_; ///< little group of the current q-point

    /**
     * @brief Get irreps for each q-point.
     *
     * Currently fills a fully-symmetric placeholder (one A1 irrep per
     * q-point); the LittleGroup decomposition is implemented in Phase 3.
     *
     * @param ucell unit cell
     * @param symm symmetry object
     */
    void get_irreps(const UnitCell& ucell, const ModuleSymmetry::Symmetry& symm);
};

} // namespace ModuleCell

#endif // QLIST_H