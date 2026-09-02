/**
 * @file little_group.h
 * @brief Little group (group of the wavevector) of a q-point.
 * @note Added 2026-08-14 for the DFPT q-point irreducible-representation
 *       analysis. The little group of a q-point is the subgroup of the
 *       crystal space group whose rotations R satisfy R q ≡ q modulo a
 *       reciprocal-lattice vector; its irreducible representations classify
 *       the atomic displacement patterns at that q-point (phonon irreps).
 *       This version exposes the interface and the little-group operation
 *       list; the full irrep table / projection-operator decomposition is
 *       added in a later iteration.
 */
#ifndef LITTLE_GROUP_H
#define LITTLE_GROUP_H

#include "source_base/vector3.h"
#include <vector>

namespace ModuleSymmetry
{

class Symmetry;

/**
 * @brief Little group (group of the wavevector) of a q-point.
 *
 * QList aggregates one LittleGroup per q-point. The operation list is built
 * from the reciprocal-space rotation matrices (Symmetry::kgmatrix) of the
 * crystal space group.
 */
class LittleGroup
{
  public:
    LittleGroup() = default;
    ~LittleGroup() = default;

    /**
     * @brief Set the q-point and determine its little group.
     *
     * The little group consists of the operations R with R q - q a vector of
     * integers (a reciprocal-lattice vector in direct coordinates), within the
     * symmetry tolerance.
     *
     * @param q    q-point in direct (fractional) coordinates
     * @param symm symmetry of the system (kgmatrix / nrotk / epsilon)
     */
    void set_q(const ModuleBase::Vector3<double>& q, const Symmetry& symm);

    /// @brief Number of irreducible representations of the little group.
    ///        Placeholder: returns 1 (the fully-symmetric A1) until the
    ///        full irrep table is implemented.
    int get_nirr() const { return nirr_; }

    /**
     * @brief Representative basis modes of an irrep.
     *        Placeholder: empty until the projection-operator decomposition
     *        is implemented.
     * @param irrep irrep index
     */
    std::vector<int> get_mode_basis(int irrep) const;

    /// @brief Indices (into Symmetry::kgmatrix / Symmetry::gtrans) of the
    ///        little-group operations.
    const std::vector<int>& get_little_group_ops() const { return little_group_ops_; }

    /// @brief The current q-point (direct coordinates).
    ModuleBase::Vector3<double> get_q() const { return q_; }

  private:
    ModuleBase::Vector3<double> q_;
    std::vector<int> little_group_ops_;
    int nirr_ = 1; ///< placeholder: fully-symmetric A1
};

} // namespace ModuleSymmetry

#endif // LITTLE_GROUP_H
