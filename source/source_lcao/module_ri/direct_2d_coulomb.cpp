#include "direct_2d_coulomb.h"

#include "source_base/constants.h"
#include "source_base/module_external/lapack_connector.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace Direct2dCoulomb
{
namespace
{
double column_norm(const ModuleBase::Matrix3& matrix, const int column)
{
    if (column == 0)
    {
        return std::sqrt(matrix.e11 * matrix.e11 + matrix.e21 * matrix.e21 + matrix.e31 * matrix.e31);
    }
    if (column == 1)
    {
        return std::sqrt(matrix.e12 * matrix.e12 + matrix.e22 * matrix.e22 + matrix.e32 * matrix.e32);
    }
    return std::sqrt(matrix.e13 * matrix.e13 + matrix.e23 * matrix.e23 + matrix.e33 * matrix.e33);
}

std::string json_escape(const std::string& input)
{
    std::string output;
    output.reserve(input.size());
    for (const char c: input)
    {
        if (c == '\\' || c == '"')
        {
            output.push_back('\\');
        }
        output.push_back(c);
    }
    return output;
}

std::complex<double> minus_i_to_power(const int l)
{
    switch ((l % 4 + 4) % 4)
    {
    case 0:
        return {1.0, 0.0};
    case 1:
        return {0.0, -1.0};
    case 2:
        return {-1.0, 0.0};
    default:
        return {0.0, 1.0};
    }
}

double interpolate_radial(const RadialTable& radial, const double k_bohr)
{
    if (radial.l < 0 || radial.dk <= 0.0 || radial.values == nullptr || radial.size < 4)
    {
        throw std::invalid_argument("Direct G-sum Coulomb received an invalid auxiliary radial table.");
    }
    const double last_supported = radial.dk * static_cast<double>(radial.size - 4);
    if (k_bohr < -1.0e-14 || k_bohr > last_supported + 1.0e-12)
    {
        throw std::out_of_range("Direct G-sum cutoff exceeds the auxiliary reciprocal radial table.");
    }
    const double position = std::max(0.0, k_bohr) / radial.dk;
    const std::size_t iq = static_cast<std::size_t>(position);
    if (iq > radial.size - 4)
    {
        return 0.0;
    }
    const double x0 = position - static_cast<double>(iq);
    const double x1 = 1.0 - x0;
    const double x2 = 2.0 - x0;
    const double x3 = 3.0 - x0;
    return x1 * x2 * (radial.values[iq] * x3 + radial.values[iq + 3] * x0) / 6.0
           + x0 * x3 * (radial.values[iq + 1] * x2 - radial.values[iq + 2] * x1) / 2.0;
}

double associated_legendre(const int l, const int m, const double x)
{
    double pmm = 1.0;
    if (m > 0)
    {
        const double root = std::sqrt(std::max(0.0, 1.0 - x * x));
        double factor = 1.0;
        for (int i = 1; i <= m; ++i)
        {
            pmm *= -factor * root;
            factor += 2.0;
        }
    }
    if (l == m)
    {
        return pmm;
    }
    double pmmp1 = x * static_cast<double>(2 * m + 1) * pmm;
    if (l == m + 1)
    {
        return pmmp1;
    }
    double pll = 0.0;
    for (int ll = m + 2; ll <= l; ++ll)
    {
        pll = (static_cast<double>(2 * ll - 1) * x * pmmp1
               - static_cast<double>(ll + m - 1) * pmm)
              / static_cast<double>(ll - m);
        pmm = pmmp1;
        pmmp1 = pll;
    }
    return pll;
}

std::vector<double> real_spherical_harmonics(const int lmax,
                                              const ModuleBase::Vector3<double>& vector)
{
    if (lmax < 0 || lmax > 30)
    {
        throw std::invalid_argument("Direct G-sum Coulomb supports auxiliary angular momentum 0 <= l <= 30.");
    }
    const double radius = vector.norm();
    if (radius <= 1.0e-14)
    {
        throw std::invalid_argument("Real spherical harmonics require non-zero momentum.");
    }
    const double cos_theta = std::max(-1.0, std::min(1.0, vector.z / radius));
    const double phi = std::atan2(vector.y, vector.x);
    std::vector<double> result(static_cast<std::size_t>((lmax + 1) * (lmax + 1)), 0.0);
    for (int l = 0; l <= lmax; ++l)
    {
        const double norm0 = std::sqrt(static_cast<double>(2 * l + 1) / ModuleBase::FOUR_PI);
        result[static_cast<std::size_t>(l * l)] = norm0 * associated_legendre(l, 0, cos_theta);
        for (int m = 1; m <= l; ++m)
        {
            const double log_ratio = std::lgamma(static_cast<double>(l - m + 1))
                                     - std::lgamma(static_cast<double>(l + m + 1));
            const double norm = std::sqrt(2.0 * static_cast<double>(2 * l + 1) / ModuleBase::FOUR_PI
                                          * std::exp(log_ratio));
            const double radial = associated_legendre(l, m, cos_theta);
            result[static_cast<std::size_t>(l * l + 2 * m - 1)] = norm * radial * std::cos(m * phi);
            result[static_cast<std::size_t>(l * l + 2 * m)] = norm * radial * std::sin(m * phi);
        }
    }
    return result;
}
} // namespace

ModuleBase::Vector3<double> cartesian_to_direct(
    const ModuleBase::Matrix3& reciprocal,
    const ModuleBase::Vector3<double>& cartesian)
{
    if (std::abs(reciprocal.Det()) <= 1.0e-14)
    {
        throw std::invalid_argument("Direct G-sum Coulomb received a singular reciprocal lattice.");
    }
    return cartesian * reciprocal.Inverse();
}

std::vector<ReciprocalPoint> enumerate_reciprocal_points(
    const ModuleBase::Matrix3& reciprocal,
    const ModuleBase::Vector3<double>& q_direct,
    const double tpiba,
    const double ecut_ry)
{
    if (!std::isfinite(tpiba) || tpiba <= 0.0)
    {
        throw std::invalid_argument("Direct G-sum Coulomb requires positive tpiba.");
    }
    if (!std::isfinite(ecut_ry) || ecut_ry <= 0.0)
    {
        throw std::invalid_argument("Direct G-sum Coulomb requires a positive reciprocal cutoff.");
    }
    if (std::abs(reciprocal.Det()) <= std::numeric_limits<double>::epsilon())
    {
        throw std::invalid_argument("Direct G-sum Coulomb received a singular reciprocal lattice.");
    }

    const double cart_radius = std::sqrt(ecut_ry) / tpiba;
    const ModuleBase::Matrix3 inverse = reciprocal.Inverse();
    std::array<int, 3> lower{};
    std::array<int, 3> upper{};
    for (int axis = 0; axis < 3; ++axis)
    {
        const double direct_bound = cart_radius * column_norm(inverse, axis);
        lower[static_cast<std::size_t>(axis)]
            = static_cast<int>(std::ceil(-direct_bound - q_direct[axis] - 1.0e-12));
        upper[static_cast<std::size_t>(axis)]
            = static_cast<int>(std::floor(direct_bound - q_direct[axis] + 1.0e-12));
    }

    std::vector<ReciprocalPoint> result;
    for (int i = lower[0]; i <= upper[0]; ++i)
    {
        for (int j = lower[1]; j <= upper[1]; ++j)
        {
            for (int k = lower[2]; k <= upper[2]; ++k)
            {
                const ModuleBase::Vector3<double> qg_direct(q_direct.x + i,
                                                            q_direct.y + j,
                                                            q_direct.z + k);
                const ModuleBase::Vector3<double> qg_cart = qg_direct * reciprocal;
                const double energy_ry = qg_cart.norm2() * tpiba * tpiba;
                if (energy_ry <= ecut_ry + 1.0e-12)
                {
                    result.push_back({{i, j, k}, qg_cart, std::sqrt(std::max(0.0, energy_ry))});
                }
            }
        }
    }
    std::sort(result.begin(), result.end(), [](const ReciprocalPoint& lhs, const ReciprocalPoint& rhs) {
        const double lhs_norm = lhs.q_plus_g_cart.norm2();
        const double rhs_norm = rhs.q_plus_g_cart.norm2();
        if (std::abs(lhs_norm - rhs_norm) > 1.0e-14)
        {
            return lhs_norm < rhs_norm;
        }
        return lhs.g_direct < rhs.g_direct;
    });
    return result;
}

std::vector<QuadratureNode> gauss_legendre_interval(const int order,
                                                    const double lower_bohr,
                                                    const double upper_bohr)
{
    if (order <= 0 || !std::isfinite(lower_bohr) || !std::isfinite(upper_bohr)
        || !(lower_bohr < upper_bohr))
    {
        throw std::invalid_argument("Direct mixed-2D Coulomb received an invalid quadrature interval.");
    }
    const double midpoint = 0.5 * (lower_bohr + upper_bohr);
    const double half_width = 0.5 * (upper_bohr - lower_bohr);
    std::vector<QuadratureNode> result(static_cast<std::size_t>(order));
    const int half_order = (order + 1) / 2;
    for (int i = 0; i < half_order; ++i)
    {
        double root = std::cos(ModuleBase::PI * (static_cast<double>(i) + 0.75)
                               / (static_cast<double>(order) + 0.5));
        double derivative = 0.0;
        for (int iteration = 0; iteration < 100; ++iteration)
        {
            double pnm2 = 1.0;
            double pnm1 = root;
            for (int n = 2; n <= order; ++n)
            {
                const double pn = ((2.0 * n - 1.0) * root * pnm1 - (n - 1.0) * pnm2)
                                  / static_cast<double>(n);
                pnm2 = pnm1;
                pnm1 = pn;
            }
            const double pn = order == 1 ? root : pnm1;
            const double pnm1_final = order == 1 ? 1.0 : pnm2;
            derivative = static_cast<double>(order) * (root * pn - pnm1_final)
                         / (root * root - 1.0);
            const double update = pn / derivative;
            root -= update;
            if (std::abs(update) < 4.0 * std::numeric_limits<double>::epsilon())
            {
                break;
            }
        }
        const double weight = 2.0 / ((1.0 - root * root) * derivative * derivative);
        const int left = i;
        const int right = order - 1 - i;
        result[static_cast<std::size_t>(left)] = {midpoint - half_width * root, half_width * weight};
        result[static_cast<std::size_t>(right)] = {midpoint + half_width * root, half_width * weight};
    }
    return result;
}

std::vector<ReciprocalLine2D> enumerate_reciprocal_lines(
    const ModuleBase::Matrix3& reciprocal,
    const ModuleBase::Vector3<double>& q_direct,
    const double tpiba,
    const double ecut_ry)
{
    if (!std::isfinite(tpiba) || tpiba <= 0.0 || !std::isfinite(ecut_ry) || ecut_ry <= 0.0)
    {
        throw std::invalid_argument("Direct mixed-2D Coulomb requires positive tpiba and cutoff.");
    }
    if (std::abs(q_direct.z) > 1.0e-12)
    {
        throw std::invalid_argument("Direct mixed-2D Coulomb requires qz=0.");
    }
    if (std::abs(reciprocal.Det()) <= std::numeric_limits<double>::epsilon())
    {
        throw std::invalid_argument("Direct mixed-2D Coulomb received a singular reciprocal lattice.");
    }

    const double cart_radius = std::sqrt(ecut_ry) / tpiba;
    const ModuleBase::Matrix3 inverse = reciprocal.Inverse();
    std::array<int, 2> lower{};
    std::array<int, 2> upper{};
    for (int axis = 0; axis < 2; ++axis)
    {
        const double direct_bound = cart_radius * column_norm(inverse, axis);
        lower[static_cast<std::size_t>(axis)]
            = static_cast<int>(std::ceil(-direct_bound - q_direct[axis] - 1.0e-12));
        upper[static_cast<std::size_t>(axis)]
            = static_cast<int>(std::floor(direct_bound - q_direct[axis] + 1.0e-12));
    }

    std::vector<ReciprocalLine2D> result;
    for (int i = lower[0]; i <= upper[0]; ++i)
    {
        for (int j = lower[1]; j <= upper[1]; ++j)
        {
            const ModuleBase::Vector3<double> qg_direct(q_direct.x + i, q_direct.y + j, 0.0);
            ModuleBase::Vector3<double> qg_cart = qg_direct * reciprocal;
            if (std::abs(qg_cart.z) > 1.0e-10)
            {
                throw std::invalid_argument(
                    "Direct mixed-2D Coulomb currently requires the reciprocal slab plane to be xy.");
            }
            qg_cart.z = 0.0;
            const double k_parallel_sq = qg_cart.norm2() * tpiba * tpiba;
            if (k_parallel_sq < ecut_ry - 1.0e-14)
            {
                result.push_back({{i, j},
                                  qg_cart,
                                  std::sqrt(std::max(0.0, k_parallel_sq)),
                                  std::sqrt(std::max(0.0, ecut_ry - k_parallel_sq))});
            }
        }
    }
    std::sort(result.begin(), result.end(), [](const ReciprocalLine2D& lhs, const ReciprocalLine2D& rhs) {
        if (std::abs(lhs.k_parallel_bohr - rhs.k_parallel_bohr) > 1.0e-14)
        {
            return lhs.k_parallel_bohr < rhs.k_parallel_bohr;
        }
        return lhs.g_direct < rhs.g_direct;
    });
    return result;
}

std::vector<MixedFourierPoint> expand_mixed_fourier_points(
    const std::vector<ReciprocalLine2D>& lines,
    const double tpiba,
    const int kz_order)
{
    if (!std::isfinite(tpiba) || tpiba <= 0.0 || kz_order <= 0)
    {
        throw std::invalid_argument("Direct mixed-2D Coulomb received invalid kz quadrature settings.");
    }
    std::vector<MixedFourierPoint> result;
    result.reserve(lines.size() * static_cast<std::size_t>(kz_order));
    for (const auto& line: lines)
    {
        if (line.k_parallel_bohr <= 1.0e-14)
        {
            throw std::invalid_argument(
                "Direct mixed-2D regular quadrature requires separate treatment of the Gamma G_parallel=0 line.");
        }
        const double theta_max = std::atan(line.kz_max_bohr / line.k_parallel_bohr);
        const auto quadrature = gauss_legendre_interval(kz_order, -theta_max, theta_max);
        for (const auto& node: quadrature)
        {
            const double tangent = std::tan(node.coordinate_bohr);
            const double kz_bohr = line.k_parallel_bohr * tangent;
            const double jacobian = line.k_parallel_bohr * (1.0 + tangent * tangent);
            ModuleBase::Vector3<double> k_cart = line.q_plus_g_parallel_cart;
            k_cart.z = kz_bohr / tpiba;
            const double k_bohr = line.k_parallel_bohr * std::sqrt(1.0 + tangent * tangent);
            result.push_back({{{line.g_direct[0], line.g_direct[1], 0}, k_cart, k_bohr},
                              node.weight_bohr * jacobian});
        }
    }
    return result;
}

std::vector<InPlaneQuadraturePoint> gamma_plane_quadrature(
    const ModuleBase::Matrix3& reciprocal,
    const double tpiba,
    const std::array<int, 3>& nmp,
    const int plane_order,
    const double massidda_chi)
{
    if (!std::isfinite(tpiba) || tpiba <= 0.0 || nmp[0] <= 0 || nmp[1] <= 0
        || nmp[2] != 1 || plane_order <= 0 || plane_order % 2 != 0
        || !std::isfinite(massidda_chi) || massidda_chi <= 0.0)
    {
        throw std::invalid_argument(
            "Direct mixed-2D Gamma quadrature requires a positive even plane order, "
            "a 2D mesh, and positive Massidda chi.");
    }
    if (std::abs(reciprocal.Det()) <= std::numeric_limits<double>::epsilon())
    {
        throw std::invalid_argument("Direct mixed-2D Gamma quadrature received a singular lattice.");
    }

    const double half_u = 0.5 / static_cast<double>(nmp[0]);
    const double half_v = 0.5 / static_cast<double>(nmp[1]);
    const auto u_nodes = gauss_legendre_interval(plane_order, -half_u, half_u);
    const auto v_nodes = gauss_legendre_interval(plane_order, -half_v, half_v);
    std::vector<InPlaneQuadraturePoint> result;
    result.reserve(static_cast<std::size_t>(plane_order * plane_order));
    double inverse_q_average = 0.0;
    for (const auto& u: u_nodes)
    {
        for (const auto& v: v_nodes)
        {
            ModuleBase::Vector3<double> q_cart
                = ModuleBase::Vector3<double>(u.coordinate_bohr, v.coordinate_bohr, 0.0)
                  * reciprocal;
            if (std::abs(q_cart.z) > 1.0e-10)
            {
                throw std::invalid_argument(
                    "Direct mixed-2D Gamma quadrature requires the reciprocal slab plane to be xy.");
            }
            q_cart.z = 0.0;
            const double q_bohr = q_cart.norm() * tpiba;
            if (q_bohr <= 1.0e-14)
            {
                throw std::invalid_argument(
                    "Direct mixed-2D Gamma quadrature contains q=0; use an even plane order.");
            }
            const double weight = u.weight_bohr * v.weight_bohr
                                  * static_cast<double>(nmp[0] * nmp[1]);
            result.push_back({q_cart, q_bohr, weight});
            inverse_q_average += weight / q_bohr;
        }
    }

    // Uniformly rescale the mesh-cell quadrature so that its leading 1/q
    // average is exactly the same Massidda value used by the Ewald route.
    // The weights stay positive, hence the resulting Gamma contribution can
    // still be accumulated as a positive Gram matrix.
    const double scale = inverse_q_average / massidda_chi;
    if (!std::isfinite(scale) || scale <= 0.0)
    {
        throw std::runtime_error("Direct mixed-2D Gamma quadrature produced an invalid scale.");
    }
    for (auto& point: result)
    {
        point.q_parallel_cart = point.q_parallel_cart * scale;
        point.q_bohr *= scale;
    }
    return result;
}

std::vector<MixedFourierPoint> expand_gamma_zero_line(
    const std::vector<InPlaneQuadraturePoint>& plane_points,
    const double tpiba,
    const double ecut_ry,
    const int kz_order)
{
    if (plane_points.empty() || !std::isfinite(tpiba) || tpiba <= 0.0
        || !std::isfinite(ecut_ry) || ecut_ry <= 0.0 || kz_order <= 0)
    {
        throw std::invalid_argument("Direct mixed-2D Gamma expansion received invalid settings.");
    }
    std::vector<MixedFourierPoint> result;
    result.reserve(plane_points.size() * static_cast<std::size_t>(kz_order));
    for (const auto& plane: plane_points)
    {
        if (!std::isfinite(plane.q_bohr) || plane.q_bohr <= 1.0e-14
            || !std::isfinite(plane.weight) || plane.weight <= 0.0
            || plane.q_bohr * plane.q_bohr >= ecut_ry)
        {
            throw std::invalid_argument(
                "Direct mixed-2D Gamma plane point is invalid or outside the reciprocal cutoff.");
        }
        const double kz_max_bohr = std::sqrt(ecut_ry - plane.q_bohr * plane.q_bohr);
        const double theta_max = std::atan(kz_max_bohr / plane.q_bohr);
        const auto quadrature = gauss_legendre_interval(kz_order, -theta_max, theta_max);
        for (const auto& node: quadrature)
        {
            const double tangent = std::tan(node.coordinate_bohr);
            const double kz_bohr = plane.q_bohr * tangent;
            const double jacobian = plane.q_bohr * (1.0 + tangent * tangent);
            ModuleBase::Vector3<double> k_cart = plane.q_parallel_cart;
            k_cart.z = kz_bohr / tpiba;
            const double k_bohr = plane.q_bohr * std::sqrt(1.0 + tangent * tangent);
            result.push_back({{{0, 0, 0}, k_cart, k_bohr},
                              plane.weight * node.weight_bohr * jacobian});
        }
    }
    return result;
}

std::vector<std::complex<double>> fourier_coefficients(
    const std::vector<BasisFunction>& basis,
    const ReciprocalPoint& point,
    const double omega)
{
    if (!std::isfinite(omega) || omega <= 0.0)
    {
        throw std::invalid_argument("Direct G-sum Coulomb requires a positive cell volume.");
    }
    int lmax = 0;
    for (const BasisFunction& function: basis)
    {
        if (function.radial.l < 0 || function.m < 0 || function.m > 2 * function.radial.l)
        {
            throw std::invalid_argument("Direct G-sum Coulomb received an invalid auxiliary angular index.");
        }
        lmax = std::max(lmax, function.radial.l);
    }

    std::vector<double> harmonics(static_cast<std::size_t>((lmax + 1) * (lmax + 1)), 0.0);
    const bool zero_momentum = point.k_bohr <= 1.0e-14;
    if (zero_momentum)
    {
        harmonics[0] = ModuleBase::SQRT_INVERSE_FOUR_PI;
    }
    else
    {
        harmonics = real_spherical_harmonics(lmax, point.q_plus_g_cart);
    }

    const double normalization = std::sqrt(std::pow(ModuleBase::TWO_PI, 3) / omega);
    std::vector<std::complex<double>> coefficients(basis.size(), {0.0, 0.0});
    for (std::size_t i = 0; i < basis.size(); ++i)
    {
        const BasisFunction& function = basis[i];
        if (zero_momentum && function.radial.l != 0)
        {
            continue;
        }
        const double radial_value = interpolate_radial(function.radial, point.k_bohr);
        const std::size_t lm = static_cast<std::size_t>(function.radial.l * function.radial.l + function.m);
        const double phase_arg = -ModuleBase::TWO_PI * (point.q_plus_g_cart * function.tau);
        const std::complex<double> center_phase(std::cos(phase_arg), std::sin(phase_arg));
        coefficients[i] = normalization * minus_i_to_power(function.radial.l) * harmonics[lm]
                          * radial_value * center_phase;
    }
    return coefficients;
}

Result build_coulomb_matrix(const std::vector<BasisFunction>& basis,
                            const std::vector<ReciprocalPoint>& points,
                            const double omega,
                            const double alpha,
                            const double chi,
                            const bool is_gamma)
{
    const auto start = std::chrono::steady_clock::now();
    if (basis.empty())
    {
        throw std::invalid_argument("Direct G-sum Coulomb requires a non-empty auxiliary basis.");
    }
    if (!std::isfinite(alpha) || alpha < 0.0)
    {
        throw std::invalid_argument("Direct G-sum Coulomb requires a non-negative Fock alpha.");
    }
    if (!std::isfinite(chi) || chi < 0.0)
    {
        throw std::invalid_argument(
            "Direct G-sum Coulomb requires a non-negative finite singular correction.");
    }

    const std::size_t dimension = basis.size();
    std::vector<bool> index_seen(dimension, false);
    for (const BasisFunction& function: basis)
    {
        if (function.global_index >= dimension || index_seen[function.global_index])
        {
            throw std::invalid_argument("Direct G-sum auxiliary indices must be a contiguous permutation.");
        }
        index_seen[function.global_index] = true;
    }

    std::vector<const ReciprocalPoint*> regular_points;
    regular_points.reserve(points.size());
    bool zero_found = false;
    for (const ReciprocalPoint& point: points)
    {
        const bool gamma_integer_zero = is_gamma
                                        && point.g_direct == std::array<int, 3>{0, 0, 0};
        if (gamma_integer_zero || point.k_bohr <= 1.0e-14)
        {
            zero_found = true;
        }
        else
        {
            regular_points.push_back(&point);
        }
    }
    if (zero_found && !is_gamma)
    {
        throw std::invalid_argument("A non-Gamma direct G-sum set contains zero momentum.");
    }

    Result result;
    result.matrix.dimension = dimension;
    result.matrix.values.assign(dimension * dimension, {0.0, 0.0});
    result.regular_g_count = regular_points.size();

    if (!regular_points.empty() && alpha > 0.0)
    {
        const std::size_t ng = regular_points.size();
        if (dimension > static_cast<std::size_t>(std::numeric_limits<int>::max())
            || ng > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        {
            throw std::overflow_error("Direct G-sum matrix dimensions exceed the BLAS integer range.");
        }
        std::vector<std::complex<double>> gram_rows(dimension * ng, {0.0, 0.0});
        for (std::size_t ig = 0; ig < ng; ++ig)
        {
            const ReciprocalPoint& point = *regular_points[ig];
            const auto rho = fourier_coefficients(basis, point, omega);
            const double scale = std::sqrt(alpha * ModuleBase::FOUR_PI
                                           / (point.k_bohr * point.k_bohr));
            for (std::size_t ib = 0; ib < basis.size(); ++ib)
            {
                gram_rows[basis[ib].global_index * ng + ig] = scale * std::conj(rho[ib]);
            }
        }
        LapackConnector::herk('U',
                              'N',
                              static_cast<int>(dimension),
                              static_cast<int>(ng),
                              1.0,
                              gram_rows.data(),
                              static_cast<int>(ng),
                              0.0,
                              result.matrix.values.data(),
                              static_cast<int>(dimension));
        for (std::size_t row = 0; row < dimension; ++row)
        {
            result.matrix(row, row) = {result.matrix(row, row).real(), 0.0};
            for (std::size_t column = row + 1; column < dimension; ++column)
            {
                result.matrix(column, row) = std::conj(result.matrix(row, column));
            }
        }
    }

    if (is_gamma && chi != 0.0 && alpha != 0.0)
    {
        // Singular_Value returns chi before the Coulomb 4*pi prefactor.
        const ReciprocalPoint zero{{0, 0, 0}, ModuleBase::Vector3<double>(), 0.0};
        const auto rho0 = fourier_coefficients(basis, zero, omega);
        for (std::size_t mu = 0; mu < basis.size(); ++mu)
        {
            for (std::size_t nu = 0; nu < basis.size(); ++nu)
            {
                result.matrix(basis[mu].global_index, basis[nu].global_index)
                    += alpha * ModuleBase::FOUR_PI * chi * std::conj(rho0[mu]) * rho0[nu];
            }
        }
    }

    for (std::size_t row = 0; row < dimension; ++row)
    {
        for (std::size_t column = 0; column < dimension; ++column)
        {
            result.hermitian_residual
                = std::max(result.hermitian_residual,
                           std::abs(result.matrix(row, column) - std::conj(result.matrix(column, row))));
        }
    }
    result.elapsed_seconds
        = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    return result;
}

Result build_gram_matrix(const std::vector<BasisFunction>& basis,
                               const std::vector<MixedFourierPoint>& points,
                               const double omega,
                               const double lz_bohr,
                               const double alpha)
{
    const auto start = std::chrono::steady_clock::now();
    if (basis.empty() || points.empty())
    {
        throw std::invalid_argument("Direct mixed-2D Coulomb requires basis functions and quadrature points.");
    }
    if (!std::isfinite(omega) || omega <= 0.0 || !std::isfinite(lz_bohr) || lz_bohr <= 0.0
        || !std::isfinite(alpha) || alpha <= 0.0)
    {
        throw std::invalid_argument("Direct mixed-2D Coulomb received invalid normalization parameters.");
    }

    const std::size_t dimension = basis.size();
    std::vector<bool> index_seen(dimension, false);
    for (const BasisFunction& function: basis)
    {
        if (function.global_index >= dimension || index_seen[function.global_index])
        {
            throw std::invalid_argument("Direct mixed-2D auxiliary indices must be a contiguous permutation.");
        }
        index_seen[function.global_index] = true;
    }
    if (dimension > static_cast<std::size_t>(std::numeric_limits<int>::max())
        || points.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
    {
        throw std::overflow_error("Direct mixed-2D matrix dimensions exceed the BLAS integer range.");
    }

    Result result;
    result.matrix.dimension = dimension;
    result.matrix.values.assign(dimension * dimension, {0.0, 0.0});
    result.regular_g_count = points.size();
    const std::size_t npoints = points.size();
    std::vector<std::complex<double>> gram_rows(dimension * npoints, {0.0, 0.0});
    for (std::size_t ipoint = 0; ipoint < npoints; ++ipoint)
    {
        const auto& mixed_point = points[ipoint];
        if (mixed_point.point.k_bohr <= 1.0e-14 || mixed_point.kz_weight_bohr <= 0.0)
        {
            throw std::invalid_argument(
                "Direct mixed-2D regular quadrature cannot contain the Gamma singular point.");
        }
        const auto rho = fourier_coefficients(basis, mixed_point.point, omega);
        const double continuum_weight = mixed_point.kz_weight_bohr * lz_bohr / ModuleBase::TWO_PI;
        const double scale = std::sqrt(alpha * ModuleBase::FOUR_PI * continuum_weight
                                       / (mixed_point.point.k_bohr * mixed_point.point.k_bohr));
        for (std::size_t ib = 0; ib < basis.size(); ++ib)
        {
            gram_rows[basis[ib].global_index * npoints + ipoint] = scale * std::conj(rho[ib]);
        }
    }
    LapackConnector::herk('U',
                          'N',
                          static_cast<int>(dimension),
                          static_cast<int>(npoints),
                          1.0,
                          gram_rows.data(),
                          static_cast<int>(npoints),
                          0.0,
                          result.matrix.values.data(),
                          static_cast<int>(dimension));
    for (std::size_t row = 0; row < dimension; ++row)
    {
        result.matrix(row, row) = {result.matrix(row, row).real(), 0.0};
        for (std::size_t column = row + 1; column < dimension; ++column)
        {
            result.matrix(column, row) = std::conj(result.matrix(row, column));
        }
    }
    for (std::size_t row = 0; row < dimension; ++row)
    {
        for (std::size_t column = 0; column < dimension; ++column)
        {
            result.hermitian_residual
                = std::max(result.hermitian_residual,
                           std::abs(result.matrix(row, column) - std::conj(result.matrix(column, row))));
        }
    }
    result.elapsed_seconds
        = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    return result;
}

std::vector<std::size_t> atom_global_indices(
    const std::vector<int>& atom_types,
    const std::vector<std::vector<std::size_t>>& type_traversal_to_local)
{
    std::vector<std::size_t> result;
    std::size_t offset = 0;
    for (const int type: atom_types)
    {
        if (type < 0 || static_cast<std::size_t>(type) >= type_traversal_to_local.size())
        {
            throw std::invalid_argument("Direct G-sum Coulomb received an invalid atom type.");
        }
        const auto& permutation = type_traversal_to_local[static_cast<std::size_t>(type)];
        std::vector<bool> seen(permutation.size(), false);
        for (const std::size_t local: permutation)
        {
            if (local >= permutation.size() || seen[local])
            {
                throw std::invalid_argument(
                    "Direct G-sum Coulomb type-local auxiliary ordering is not bijective.");
            }
            seen[local] = true;
            result.push_back(offset + local);
        }
        offset += permutation.size();
    }
    return result;
}

void write_matrix_market(const DenseMatrix& matrix, const std::string& path)
{
    if (matrix.dimension == 0 || matrix.values.size() != matrix.dimension * matrix.dimension)
    {
        throw std::invalid_argument("Direct G-sum diagnostic received an invalid dense matrix.");
    }
    std::ofstream stream(path, std::ios::out | std::ios::trunc);
    if (!stream.good())
    {
        throw std::runtime_error("Failed to open direct G-sum MatrixMarket output: " + path);
    }
    stream << "%%MatrixMarket matrix coordinate complex general\n";
    stream << matrix.dimension << ' ' << matrix.dimension << ' '
           << matrix.dimension * matrix.dimension << '\n';
    stream << std::scientific << std::setprecision(17);
    for (std::size_t row = 0; row < matrix.dimension; ++row)
    {
        for (std::size_t column = 0; column < matrix.dimension; ++column)
        {
            const auto value = matrix(row, column);
            stream << row + 1 << ' ' << column + 1 << ' ' << value.real() << ' ' << value.imag() << '\n';
        }
    }
    if (!stream.good())
    {
        throw std::runtime_error("Failed while writing direct G-sum MatrixMarket output: " + path);
    }
}

void write_reader_v1(const DenseMatrix& matrix,
                     const std::vector<int>& atom_naux,
                     const int iq_one_based,
                     const int mpi_rank,
                     const std::string& prefix)
{
    constexpr std::int32_t marker = -20129433;
    constexpr std::int32_t complex_flag = 1;
    static_assert(sizeof(std::complex<double>) == 2 * sizeof(double),
                  "LibRPA v1 Coulomb output expects complex<double> as two doubles.");
    if (matrix.dimension == 0 || matrix.values.size() != matrix.dimension * matrix.dimension)
    {
        throw std::invalid_argument("Direct G-sum reader-v1 output received an invalid dense matrix.");
    }
    if (atom_naux.empty() || iq_one_based <= 0 || mpi_rank < 0 || prefix.empty())
    {
        throw std::invalid_argument("Direct G-sum reader-v1 output received invalid metadata.");
    }
    std::size_t naux = 0;
    for (const int count: atom_naux)
    {
        if (count <= 0)
        {
            throw std::invalid_argument("Direct G-sum reader-v1 output requires positive per-atom sizes.");
        }
        naux += static_cast<std::size_t>(count);
    }
    if (naux != matrix.dimension
        || naux > static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max())
        || atom_naux.size() > static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max()))
    {
        throw std::invalid_argument("Direct G-sum reader-v1 auxiliary dimensions are inconsistent.");
    }

    struct Block
    {
        std::int32_t pair_index = 0;
        std::int64_t offset = 0;
        std::size_t row_begin = 0;
        std::size_t column_begin = 0;
        std::size_t nrows = 0;
        std::size_t ncolumns = 0;
    };
    std::vector<Block> blocks;
    std::size_t row_begin = 0;
    for (std::size_t atom_i = 0; atom_i < atom_naux.size(); ++atom_i)
    {
        std::size_t column_begin = row_begin;
        for (std::size_t atom_j = atom_i; atom_j < atom_naux.size(); ++atom_j)
        {
            const std::size_t pair_index
                = atom_i * atom_naux.size() - atom_i * (atom_i - 1) / 2 + (atom_j - atom_i);
            if (pair_index > static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max()))
            {
                throw std::overflow_error("Direct G-sum reader-v1 atom-pair index exceeds int32 range.");
            }
            blocks.push_back({static_cast<std::int32_t>(pair_index),
                              0,
                              row_begin,
                              column_begin,
                              static_cast<std::size_t>(atom_naux[atom_i]),
                              static_cast<std::size_t>(atom_naux[atom_j])});
            column_begin += static_cast<std::size_t>(atom_naux[atom_j]);
        }
        row_begin += static_cast<std::size_t>(atom_naux[atom_i]);
    }
    if (blocks.size() > static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max()))
    {
        throw std::overflow_error("Direct G-sum reader-v1 block count exceeds int32 range.");
    }

    const auto checked_add = [](const std::int64_t lhs, const std::size_t values) {
        if (values > static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max())
                         / sizeof(std::complex<double>))
        {
            throw std::overflow_error("Direct G-sum reader-v1 payload exceeds int64 range.");
        }
        const std::int64_t bytes
            = static_cast<std::int64_t>(values * sizeof(std::complex<double>));
        if (lhs > std::numeric_limits<std::int64_t>::max() - bytes)
        {
            throw std::overflow_error("Direct G-sum reader-v1 file offset exceeds int64 range.");
        }
        return lhs + bytes;
    };
    std::int64_t offset = static_cast<std::int64_t>(6 * sizeof(std::int32_t))
                          + static_cast<std::int64_t>(atom_naux.size() * sizeof(std::int32_t))
                          + static_cast<std::int64_t>(blocks.size())
                                * static_cast<std::int64_t>(sizeof(std::int32_t) + sizeof(std::int64_t));
    for (auto& block: blocks)
    {
        block.offset = offset;
        offset = checked_add(offset, block.nrows * block.ncolumns);
    }

    std::ostringstream name;
    name << prefix << iq_one_based << "_rank" << mpi_rank << ".dat";
    const std::string path = name.str();
    std::ofstream stream(path, std::ios::out | std::ios::binary | std::ios::trunc);
    if (!stream.good())
    {
        throw std::runtime_error("Failed to open direct G-sum reader-v1 output: " + path);
    }
    const auto write = [&stream, &path](const void* data, const std::size_t bytes) {
        stream.write(reinterpret_cast<const char*>(data), static_cast<std::streamsize>(bytes));
        if (!stream.good())
        {
            throw std::runtime_error("Failed while writing direct G-sum reader-v1 output: " + path);
        }
    };
    const std::int32_t iq = iq_one_based;
    const std::int32_t naux_i32 = static_cast<std::int32_t>(naux);
    const std::int32_t natoms = static_cast<std::int32_t>(atom_naux.size());
    const std::int32_t nblocks = static_cast<std::int32_t>(blocks.size());
    for (const auto value: {marker, iq, naux_i32, complex_flag, natoms, nblocks})
    {
        write(&value, sizeof(value));
    }
    for (const int count: atom_naux)
    {
        const std::int32_t count_i32 = count;
        write(&count_i32, sizeof(count_i32));
    }
    for (const auto& block: blocks)
    {
        write(&block.pair_index, sizeof(block.pair_index));
        write(&block.offset, sizeof(block.offset));
    }
    for (const auto& block: blocks)
    {
        for (std::size_t row = 0; row < block.nrows; ++row)
        {
            const auto* values = matrix.values.data()
                                 + (block.row_begin + row) * matrix.dimension
                                 + block.column_begin;
            write(values, block.ncolumns * sizeof(std::complex<double>));
        }
    }
}

void write_metadata(const Result& result, const Metadata& metadata, const std::string& path)
{
    std::ofstream stream(path, std::ios::out | std::ios::trunc);
    if (!stream.good())
    {
        throw std::runtime_error("Failed to open direct G-sum metadata output: " + path);
    }
    stream << std::scientific << std::setprecision(17);
    stream << "{\n"
           << "  \"iq\": " << metadata.iq << ",\n"
           << "  \"q_direct\": [" << metadata.q_direct.x << ", " << metadata.q_direct.y << ", "
           << metadata.q_direct.z << "],\n"
           << "  \"q_cart\": [" << metadata.q_cart.x << ", " << metadata.q_cart.y << ", "
           << metadata.q_cart.z << "],\n"
           << "  \"ecut_ry\": " << metadata.ecut_ry << ",\n"
           << "  \"alpha\": " << metadata.alpha << ",\n"
           << "  \"chi\": " << metadata.chi << ",\n"
           << "  \"dimension\": " << result.matrix.dimension << ",\n"
           << "  \"regular_g_count\": " << result.regular_g_count << ",\n"
           << "  \"hermitian_residual\": " << result.hermitian_residual << ",\n"
           << "  \"elapsed_seconds\": " << result.elapsed_seconds << ",\n"
           << "  \"mpi_owner\": " << metadata.mpi_owner << ",\n"
           << "  \"mpi_size\": " << metadata.mpi_size << ",\n"
           << "  \"build_version\": \"" << json_escape(metadata.build_version) << "\"\n"
           << "}\n";
    if (!stream.good())
    {
        throw std::runtime_error("Failed while writing direct G-sum metadata output: " + path);
    }
}
std::string format_method_metadata(const MethodMetadata& metadata)
{
    if (metadata.method.empty() || !std::isfinite(metadata.ecut_ry)
        || metadata.ecut_ry <= 0.0 || metadata.kz_order <= 0
        || metadata.gamma_order <= 0 || metadata.gamma_order % 2 != 0
        || metadata.nq <= 0 || metadata.naux == 0 || metadata.source_revision.empty())
    {
        throw std::invalid_argument("invalid direct strict-2D Coulomb method metadata");
    }
    std::ostringstream output;
    output << "# ABACUS reader-v1 strict 2D Coulomb method\n"
           << "version = 1\n"
           << "method = " << metadata.method << '\n'
           << std::setprecision(17)
           << "ecut_ry = " << metadata.ecut_ry << '\n'
           << "kz_order = " << metadata.kz_order << '\n'
           << "gamma_order = " << metadata.gamma_order << '\n'
           << "nq = " << metadata.nq << '\n'
           << "naux = " << metadata.naux << '\n'
           << "source_revision = " << metadata.source_revision << '\n';
    return output.str();
}

std::string format_3d_method_metadata(const MethodMetadata3D& metadata)
{
    if (metadata.method.empty() || !std::isfinite(metadata.ecut_ry)
        || metadata.ecut_ry <= 0.0 || metadata.nq <= 0
        || metadata.naux == 0 || metadata.source_revision.empty())
    {
        throw std::invalid_argument("invalid direct 3D Coulomb method metadata");
    }
    std::ostringstream output;
    output << "# ABACUS reader-v1 3D Coulomb method\n"
           << "version = 1\n"
           << "method = " << metadata.method << '\n'
           << std::setprecision(17)
           << "ecut_ry = " << metadata.ecut_ry << '\n'
           << "nq = " << metadata.nq << '\n'
           << "naux = " << metadata.naux << '\n'
           << "source_revision = " << metadata.source_revision << '\n';
    return output.str();
}

} // namespace Direct2dCoulomb
