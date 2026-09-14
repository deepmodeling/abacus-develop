#include "xc_ncgga_radial.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace ModuleXC
{

double NcggaRadialPoint::jacobian(const int row, const int column) const
{
    const double identity = (row == column) ? 1.0 : 0.0;
    return transverse_hessian * identity + (radial_hessian - transverse_hessian) * direction[row] * direction[column];
}

NcggaRadialPoint make_ncgga_radial_point(const std::array<double, 3>& magnetization, const double eta)
{
    if (!(eta > 0.0))
    {
        throw std::invalid_argument("noncollinear GGA radial eta must be positive");
    }

    NcggaRadialPoint point;
    const double magnitude = std::sqrt(magnetization[0] * magnetization[0] + magnetization[1] * magnetization[1]
                                       + magnetization[2] * magnetization[2]);
    if (magnitude == 0.0)
    {
        return point;
    }

    for (int component = 0; component < 3; ++component)
    {
        point.direction[component] = magnetization[component] / magnitude;
    }

    if (magnitude < eta)
    {
        const double x = magnitude / eta;
        const double x2 = x * x;
        const double x3 = x2 * x;
        point.value = eta * x3 * (3.0 * x2 - 8.0 * x + 6.0);
        point.transverse_hessian = x * (15.0 * x2 - 32.0 * x + 18.0) / eta;
        point.radial_hessian = x * (60.0 * x2 - 96.0 * x + 36.0) / eta;
    }
    else
    {
        point.value = magnitude;
        point.transverse_hessian = 1.0 / magnitude;
        point.radial_hessian = 0.0;
    }

    for (int component = 0; component < 3; ++component)
    {
        point.gradient[component] = point.transverse_hessian * magnetization[component];
    }
    return point;
}

double ncgga_lca_radial_eta()
{
    return 1.0e-3;
}

double NcggaSpinMapPoint::jacobian(const int spin, const int channel) const
{
    if (channel == 0)
    {
        if (saturated)
        {
            return spin == 0 ? density_sign : 0.0;
        }
        return 0.5 * density_sign;
    }
    if (saturated)
    {
        return 0.0;
    }
    const double spin_sign = spin == 0 ? 0.5 : -0.5;
    return spin_sign * radial.gradient[channel - 1];
}

NcggaSpinMapPoint make_ncgga_spin_map_point(const double total_density, const NcggaRadialPoint& radial)
{
    NcggaSpinMapPoint point;
    point.radial = radial;
    point.absolute_density = std::abs(total_density);
    point.clipped_magnitude = std::min(radial.value, point.absolute_density);
    point.spin_density[0] = 0.5 * (point.absolute_density + point.clipped_magnitude);
    point.spin_density[1] = 0.5 * (point.absolute_density - point.clipped_magnitude);
    point.density_sign = total_density > 0.0 ? 1.0 : total_density < 0.0 ? -1.0 : 0.0;
    point.saturated = !(radial.value < point.absolute_density);
    return point;
}

} // namespace ModuleXC
