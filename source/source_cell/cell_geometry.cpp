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

bool make_slab_cell(const ModuleBase::Matrix3& lattice,
                    const double lattice_scale,
                    const int open_axis,
                    const double relative_tolerance,
                    SlabCell& cell,
                    std::string& error)
{
    error.clear();
    if (open_axis < 0 || open_axis > 2 || !std::isfinite(lattice_scale) || lattice_scale <= 0.0
        || !std::isfinite(relative_tolerance) || relative_tolerance <= 0.0)
    {
        error = "slab geometry requires open axis 0/1/2 and positive finite scale and tolerance";
        return false;
    }
    const ModuleBase::Vector3<double> a(lattice.e11, lattice.e12, lattice.e13);
    const ModuleBase::Vector3<double> b(lattice.e21, lattice.e22, lattice.e23);
    const ModuleBase::Vector3<double> c(lattice.e31, lattice.e32, lattice.e33);
    const std::array<ModuleBase::Vector3<double>, 3> vectors = {{a, b, c}};
    const int first_axis = (open_axis + 1) % 3;
    const int second_axis = (open_axis + 2) % 3;
    const ModuleBase::Vector3<double> open = vectors[open_axis] * lattice_scale;
    const ModuleBase::Vector3<double> first = vectors[first_axis] * lattice_scale;
    const ModuleBase::Vector3<double> second = vectors[second_axis] * lattice_scale;
    const double length = open.norm();
    const double first_length = first.norm();
    const double second_length = second.norm();
    const ModuleBase::Vector3<double> plane_normal = first ^ second;
    const double area = plane_normal.norm();
    if (!finite_vector(open) || !finite_vector(first) || !finite_vector(second)
        || !std::isfinite(length) || length <= 0.0 || !std::isfinite(area)
        || area <= relative_tolerance * first_length * second_length)
    {
        error = "slab geometry requires finite nonzero open and independent periodic lattice vectors";
        return false;
    }
    const double orientation = (plane_normal * open < 0.0) ? -1.0 : 1.0;
    const double normal_scale = orientation / area;
    const ModuleBase::Vector3<double> normal = plane_normal * normal_scale;
    const ModuleBase::Vector3<double> alignment = open ^ normal;
    if (alignment.norm() > relative_tolerance * length)
    {
        error = "slab open lattice vector must be perpendicular to the periodic plane";
        return false;
    }
    SlabCell candidate;
    candidate.normal = normal;
    candidate.length = length;
    candidate.area = area;
    candidate.origin = 0.5 * length;
    cell = candidate;
    return true;
}

double relative_coordinate(const ModuleBase::Vector3<double>& position, const SlabCell& cell)
{
    const double displacement = position * cell.normal - cell.origin;
    const double image_coordinate = displacement / cell.length + 0.5;
    const double image = std::floor(image_coordinate);
    return displacement - cell.length * image;
}

bool weighted_center(const std::vector<ModuleBase::Vector3<double>>& positions,
                      const std::vector<double>& weights,
                      const SlabCell& cell,
                      double& center,
                      std::string& error)
{
    error.clear();
    if (positions.empty() || positions.size() != weights.size())
    {
        error = "slab center requires matching non-empty positions and weights";
        return false;
    }
    SlabCell reference = cell;
    reference.origin = positions.front() * cell.normal;
    double total_weight = 0.0;
    double weighted_displacement = 0.0;
    for (std::size_t index = 0; index < positions.size(); ++index)
    {
        if (!finite_vector(positions[index]) || !std::isfinite(weights[index]) || weights[index] <= 0.0)
        {
            error = "slab center requires finite positions and positive finite weights";
            return false;
        }
        const double displacement = relative_coordinate(positions[index], reference);
        weighted_displacement += weights[index] * displacement;
        total_weight += weights[index];
    }
    const double unwrapped = reference.origin + weighted_displacement / total_weight;
    const double image_coordinate = unwrapped / cell.length;
    const double image = std::floor(image_coordinate);
    const double candidate = unwrapped - image * cell.length;
    if (!std::isfinite(total_weight) || !std::isfinite(weighted_displacement) || !std::isfinite(candidate))
    {
        error = "slab center accumulation must be finite";
        return false;
    }
    center = candidate;
    return true;
}

} // namespace unitcell
