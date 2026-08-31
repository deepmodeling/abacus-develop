#include "../direct_2d_coulomb.h"

#include "source_base/constants.h"

#include <array>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
void require(const bool condition, const char* message)
{
    if (!condition)
    {
        throw std::runtime_error(message);
    }
}

void require_close(const double actual,
                   const double expected,
                   const double tolerance,
                   const char* message)
{
    if (std::abs(actual - expected) > tolerance)
    {
        throw std::runtime_error(message);
    }
}

void require_positive_quadratic_forms(const Direct2dCoulomb::DenseMatrix& matrix)
{
    const std::vector<std::vector<std::complex<double>>> probes = {
        {{1.0, 0.0}, {0.0, 0.0}},
        {{0.0, 0.0}, {1.0, 0.0}},
        {{1.0, 0.0}, {1.0, 0.0}},
        {{1.0, 0.0}, {0.0, 1.0}},
        {{0.25, -0.75}, {-0.4, 0.2}},
    };
    for (const auto& vector: probes)
    {
        std::complex<double> value = 0.0;
        for (std::size_t i = 0; i < matrix.dimension; ++i)
        {
            for (std::size_t j = 0; j < matrix.dimension; ++j)
            {
                value += std::conj(vector[i]) * matrix(i, j) * vector[j];
            }
        }
        require(std::abs(value.imag()) < 1.0e-11, "3D quadratic form is not real");
        require(value.real() > -1.0e-11, "3D quadratic form is negative");
    }
}
} // namespace

int main()
{
    try
    {
        ModuleBase::Matrix3 reciprocal;
        reciprocal.e11 = 1.0;
        reciprocal.e22 = 1.0;
        reciprocal.e33 = 1.0;
        constexpr double tpiba = 1.0;
        constexpr double ecut_ry = 1.0;
        const ModuleBase::Vector3<double> gamma(0.0, 0.0, 0.0);
        const auto gamma_points = Direct2dCoulomb::enumerate_reciprocal_points(
            reciprocal, gamma, tpiba, ecut_ry);
        require(gamma_points.size() == 7,
                "cubic 3D reciprocal enumeration must contain G=0 and six unit vectors");
        require(gamma_points.front().g_direct == std::array<int, 3>{0, 0, 0},
                "3D reciprocal enumeration does not place G=0 first");

        constexpr std::size_t radial_size = 301;
        constexpr double dk = 0.01;
        std::vector<double> radial_s(radial_size, 1.0);
        constexpr double omega = 64.0;
        const std::vector<Direct2dCoulomb::BasisFunction> one_s_basis = {
            {{0, dk, radial_s.data(), radial_s.size()}, 0, 0, {0.0, 0.0, 0.0}},
        };
        constexpr double chi = 2.5;
        const auto one_s = Direct2dCoulomb::build_coulomb_matrix(
            one_s_basis, gamma_points, omega, 1.0, chi, true);
        require(one_s.regular_g_count == 6, "3D Gamma matrix did not remove only G=0");
        require(one_s.hermitian_residual < 1.0e-13, "3D one-function matrix is not Hermitian");
        const double expected_regular = 6.0 * 8.0 * std::pow(ModuleBase::PI, 3) / omega;
        const double expected_gamma = 8.0 * std::pow(ModuleBase::PI, 3) * chi / omega;
        require_close(one_s.matrix(0, 0).real(),
                      expected_regular + expected_gamma,
                      2.0e-12,
                      "3D direct Coulomb normalization does not match the analytic s result");
        require_close(one_s.matrix(0, 0).imag(),
                      0.0,
                      1.0e-14,
                      "3D direct Coulomb diagonal is complex");

        std::vector<double> radial_p(radial_size);
        for (std::size_t i = 0; i < radial_size; ++i)
        {
            const double k = dk * static_cast<double>(i);
            radial_p[i] = k * std::exp(-0.4 * k * k);
        }
        const std::vector<Direct2dCoulomb::BasisFunction> two_function_basis = {
            {{0, dk, radial_s.data(), radial_s.size()}, 0, 0, {0.0, 0.0, 0.0}},
            {{1, dk, radial_p.data(), radial_p.size()}, 0, 1, {0.17, 0.11, 0.23}},
        };
        const auto two_function = Direct2dCoulomb::build_coulomb_matrix(
            two_function_basis, gamma_points, omega, 1.0, chi, true);
        require(two_function.hermitian_residual < 1.0e-12,
                "3D two-function matrix is not Hermitian");
        require_positive_quadratic_forms(two_function.matrix);

        const ModuleBase::Vector3<double> q_direct(0.25, 0.0, 0.0);
        const auto non_gamma_points = Direct2dCoulomb::enumerate_reciprocal_points(
            reciprocal, q_direct, tpiba, ecut_ry);
        for (const auto& point: non_gamma_points)
        {
            require(point.k_bohr > 1.0e-14,
                    "non-Gamma 3D reciprocal enumeration contains zero momentum");
        }
        const auto non_gamma = Direct2dCoulomb::build_coulomb_matrix(
            two_function_basis, non_gamma_points, omega, 1.0, 0.0, false);
        require(non_gamma.hermitian_residual < 1.0e-12,
                "non-Gamma 3D matrix is not Hermitian");
        require_positive_quadratic_forms(non_gamma.matrix);

        Direct2dCoulomb::MethodMetadata3D metadata;
        metadata.method = "direct_reciprocal";
        metadata.ecut_ry = 110.0;
        metadata.nq = 65;
        metadata.naux = 1108;
        metadata.source_revision = "test-revision";
        const std::string formatted = Direct2dCoulomb::format_3d_method_metadata(metadata);
        const std::string expected_metadata =
            "# ABACUS reader-v1 3D Coulomb method\n"
            "version = 1\n"
            "method = direct_reciprocal\n"
            "ecut_ry = 110\n"
            "nq = 65\n"
            "naux = 1108\n"
            "source_revision = test-revision\n";
        require(formatted == expected_metadata, "3D method metadata is not deterministic");

        std::cout << "direct 3D reciprocal Coulomb tests passed" << std::endl;
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << std::endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
