/**
 * @file little_group.cpp
 * @brief Implementation of ModuleSymmetry::LittleGroup.
 */
#include "little_group.h"

#include "symmetry.h"

#include <cmath>

namespace ModuleSymmetry
{

void LittleGroup::set_q(const ModuleBase::Vector3<double>& q, const Symmetry& symm)
{
    this->q_ = q;
    this->little_group_ops_.clear();

    // A rotation R belongs to the little group of q when R q - q is a
    // reciprocal-lattice vector. In direct coordinates of the reciprocal
    // lattice such a vector has integer components, so we test each component
    // of (R q - q) for integrality within the symmetry tolerance. The
    // row-vector convention (q * R) matches the one used by reduce_ibz.
    const double eps = symm.epsilon;
    for (int i = 0; i < symm.nrotk; ++i)
    {
        const ModuleBase::Vector3<double> rq = q * symm.kgmatrix[i];
        const double dx = rq.x - q.x;
        const double dy = rq.y - q.y;
        const double dz = rq.z - q.z;
        const double fx = dx - std::floor(dx + 0.5); // signed distance to nearest integer
        const double fy = dy - std::floor(dy + 0.5);
        const double fz = dz - std::floor(dz + 0.5);
        if (std::abs(fx) < eps && std::abs(fy) < eps && std::abs(fz) < eps)
        {
            this->little_group_ops_.push_back(i);
        }
    }

    // Placeholder: one fully-symmetric irrep (A1) per q-point. The real
    // little-group representation analysis (kgmatrix + gtrans phase factors)
    // is implemented in a later iteration.
    this->nirr_ = 1;
}

std::vector<int> LittleGroup::get_mode_basis(int irrep) const
{
    (void)irrep;
    // Placeholder: the projection-operator basis is not implemented yet.
    return std::vector<int>();
}

} // namespace ModuleSymmetry
