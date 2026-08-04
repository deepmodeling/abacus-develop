/**
 * @file qlist.h
 * @brief QList class for managing q-points.
 * @author Mohan Chen (added on 2026-05-18)
 * @note This code is currently in the design phase and has not been
 *       put into production yet. It may change in the future.
 *       Please use this code with caution. Only developers who know
 *       what they are doing should use this code.
 */
#ifndef QLIST_H
#define QLIST_H

#include "source_base/vector3.h"
#include "module_symmetry/symmetry.h"
#include "unitcell.h"
#include <vector>

namespace ModuleCell {

/**
 * @brief QList class for managing q-points.
 */
class QList {
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
     * @brief Generate q-point mesh.
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
    int get_nq() const { return nq_; }
    
    /**
     * @brief Get q-point at given index.
     * @param idx q-point index
     * @return q-point vector
     */
    ModuleBase::Vector3<double> get_q(int idx) const { return qvec_[idx]; }
    
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

private:
    int nq_ = 0; ///< number of q-points
    std::vector<ModuleBase::Vector3<double>> qvec_; ///< q-point vectors
    std::vector<int> nirr_; ///< number of irreps for each q-point
    std::vector<std::vector<std::vector<int>>> irrep_modes_; ///< irrep modes
    
    /**
     * @brief Reduce q-points using symmetry.
     *
     * @param ucell unit cell
     * @param symm symmetry object
     */
    void reduce(UnitCell& ucell, ModuleSymmetry::Symmetry& symm);
    
    /**
     * @brief Get irreps for each q-point.
     *
     * @param ucell unit cell
     * @param symm symmetry object
     */
    void get_irreps(UnitCell& ucell, ModuleSymmetry::Symmetry& symm);
};

} // namespace ModuleCell

#endif // QLIST_H