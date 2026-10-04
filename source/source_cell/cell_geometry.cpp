#include "cell_geometry.h"

#include "source_base/matrix3.h"

#include <cmath>

namespace
{
bool finite_vector(const ModuleBase::Vector3<double>& value)
{
    return std::isfinite(value.x) && std::isfinite(value.y) && std::isfinite(value.z);
}
} // namespace

namespace unitcell
{

bool make_orthogonal_cell(const ModuleBase::Matrix3& lattice,
                          const double lattice_scale,
                          const double relative_tolerance,
                          OrthogonalCell& cell,
                          std::string& error)
{
    error.clear();
    if (!std::isfinite(lattice_scale) || lattice_scale <= 0.0
        || !std::isfinite(relative_tolerance) || relative_tolerance <= 0.0)
    {
        error = "cell scale and relative tolerance must be positive and finite";
        return false;
    }

    const ModuleBase::Vector3<double> a(lattice.e11, lattice.e12, lattice.e13);
    const ModuleBase::Vector3<double> b(lattice.e21, lattice.e22, lattice.e23);
    const ModuleBase::Vector3<double> c(lattice.e31, lattice.e32, lattice.e33);
    const std::array<ModuleBase::Vector3<double>, 3> vectors = {{a, b, c}};
    OrthogonalCell candidate;
    for (int axis = 0; axis < 3; ++axis)
    {
        const ModuleBase::Vector3<double> physical_vector = vectors[axis] * lattice_scale;
        const double length = physical_vector.norm();
        if (!finite_vector(physical_vector) || !std::isfinite(length) || length <= 0.0)
        {
            error = "cell lattice vectors must have positive finite lengths";
            return false;
        }
        candidate.lengths[axis] = length;
        candidate.axes[axis] = physical_vector / length;
    }
    for (int left = 0; left < 3; ++left)
    {
        for (int right = left + 1; right < 3; ++right)
        {
            const double projection = candidate.axes[left] * candidate.axes[right];
            if (std::abs(projection) > relative_tolerance)
            {
                error = "cell lattice vectors must be mutually perpendicular";
                return false;
            }
        }
    }
    const ModuleBase::Vector3<double> vector_sum = a + b + c;
    const double half_scale = 0.5 * lattice_scale;
    candidate.origin = vector_sum * half_scale;
    if (!finite_vector(candidate.origin))
    {
        error = "cell center must be finite";
        return false;
    }
    cell = candidate;
    return true;
}

ModuleBase::Vector3<double> relative_position(const ModuleBase::Vector3<double>& position,
                                             const OrthogonalCell& cell)
{
    const ModuleBase::Vector3<double> displacement = position - cell.origin;
    ModuleBase::Vector3<double> relative;
    for (int axis = 0; axis < 3; ++axis)
    {
        const double projection = displacement * cell.axes[axis];
        const double length = cell.lengths[axis];
        const double fractional_image = projection / length + 0.5;
        const double image = std::floor(fractional_image);
        const double wrapped = projection - length * image;
        relative += cell.axes[axis] * wrapped;
    }
    return relative;
}

bool weighted_center(const std::vector<ModuleBase::Vector3<double>>& positions,
                      const std::vector<double>& weights,
                      const OrthogonalCell& cell,
                      ModuleBase::Vector3<double>& center,
                      std::string& error)
{
    error.clear();
    if (positions.empty() || positions.size() != weights.size())
    {
        error = "weighted center requires matching non-empty positions and weights";
        return false;
    }
    OrthogonalCell reference = cell;
    reference.origin = positions.front();
    double total_weight = 0.0;
    ModuleBase::Vector3<double> weighted_displacement;
    for (std::size_t index = 0; index < positions.size(); ++index)
    {
        if (!finite_vector(positions[index]) || !std::isfinite(weights[index]) || weights[index] <= 0.0)
        {
            error = "weighted center requires finite positions and positive finite weights";
            return false;
        }
        const ModuleBase::Vector3<double> displacement = relative_position(positions[index], reference);
        weighted_displacement += displacement * weights[index];
        total_weight += weights[index];
    }
    if (!std::isfinite(total_weight) || !finite_vector(weighted_displacement))
    {
        error = "weighted center accumulation must be finite";
        return false;
    }
    const ModuleBase::Vector3<double> mean_displacement = weighted_displacement / total_weight;
    const ModuleBase::Vector3<double> unwrapped_center = positions.front() + mean_displacement;
    const ModuleBase::Vector3<double> wrapped_center = relative_position(unwrapped_center, cell);
    const ModuleBase::Vector3<double> candidate = cell.origin + wrapped_center;
    if (!finite_vector(candidate))
    {
        error = "weighted center must be finite";
        return false;
    }
    center = candidate;
    return true;
}

} // namespace unitcell
