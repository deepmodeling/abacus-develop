#include "cell_geometry.h"

#include "source_base/matrix3.h"

#include <cmath>

namespace unitcell
{

bool make_orthogonal_cell(const ModuleBase::Matrix3& lattice,
                          const double lattice_scale,
                          const double relative_tolerance,
                          OrthogonalCell& cell)
{
    const ModuleBase::Vector3<double> a(lattice.e11, lattice.e12, lattice.e13);
    const ModuleBase::Vector3<double> b(lattice.e21, lattice.e22, lattice.e23);
    const ModuleBase::Vector3<double> c(lattice.e31, lattice.e32, lattice.e33);
    const std::array<ModuleBase::Vector3<double>, 3> vectors = {{a, b, c}};
    for (int axis = 0; axis < 3; ++axis)
    {
        const ModuleBase::Vector3<double> physical_vector = vectors[axis] * lattice_scale;
        const double length = physical_vector.norm();
        cell.lengths[axis] = length;
        cell.axes[axis] = physical_vector / length;
    }
    const ModuleBase::Vector3<double> vector_sum = a + b + c;
    const double half_scale = 0.5 * lattice_scale;
    cell.origin = vector_sum * half_scale;
    for (int left = 0; left < 3; ++left)
    {
        for (int right = left + 1; right < 3; ++right)
        {
            const double projection = cell.axes[left] * cell.axes[right];
            if (std::abs(projection) > relative_tolerance)
            {
                return false;
            }
        }
    }
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

ModuleBase::Vector3<double> weighted_center(const std::vector<ModuleBase::Vector3<double>>& positions,
                                            const std::vector<double>& weights,
                                            const OrthogonalCell& cell)
{
    OrthogonalCell reference = cell;
    reference.origin = positions.front();
    double total_weight = 0.0;
    ModuleBase::Vector3<double> weighted_displacement;
    for (std::size_t index = 0; index < positions.size(); ++index)
    {
        const ModuleBase::Vector3<double> displacement = relative_position(positions[index], reference);
        weighted_displacement += displacement * weights[index];
        total_weight += weights[index];
    }
    const ModuleBase::Vector3<double> mean_displacement = weighted_displacement / total_weight;
    const ModuleBase::Vector3<double> unwrapped_center = positions.front() + mean_displacement;
    const ModuleBase::Vector3<double> wrapped_center = relative_position(unwrapped_center, cell);
    return cell.origin + wrapped_center;
}

bool make_slab_cell(const ModuleBase::Matrix3& lattice,
                    const double lattice_scale,
                    const int open_axis,
                    const double relative_tolerance,
                    SlabCell& cell)
{
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
    const ModuleBase::Vector3<double> plane_normal = first ^ second;
    const double area = plane_normal.norm();
    const double orientation = (plane_normal * open < 0.0) ? -1.0 : 1.0;
    const double normal_scale = orientation / area;
    const ModuleBase::Vector3<double> normal = plane_normal * normal_scale;
    cell.normal = normal;
    cell.length = length;
    cell.area = area;
    cell.origin = 0.5 * length;
    const ModuleBase::Vector3<double> alignment = open ^ normal;
    return alignment.norm() <= relative_tolerance * length;
}

double relative_coordinate(const ModuleBase::Vector3<double>& position, const SlabCell& cell)
{
    const double displacement = position * cell.normal - cell.origin;
    const double image_coordinate = displacement / cell.length + 0.5;
    const double image = std::floor(image_coordinate);
    return displacement - cell.length * image;
}

double weighted_center(const std::vector<ModuleBase::Vector3<double>>& positions,
                       const std::vector<double>& weights,
                       const SlabCell& cell)
{
    SlabCell reference = cell;
    reference.origin = positions.front() * cell.normal;
    double total_weight = 0.0;
    double weighted_displacement = 0.0;
    for (std::size_t index = 0; index < positions.size(); ++index)
    {
        const double displacement = relative_coordinate(positions[index], reference);
        weighted_displacement += weights[index] * displacement;
        total_weight += weights[index];
    }
    const double unwrapped = reference.origin + weighted_displacement / total_weight;
    const double image_coordinate = unwrapped / cell.length;
    const double image = std::floor(image_coordinate);
    return unwrapped - image * cell.length;
}

} // namespace unitcell
