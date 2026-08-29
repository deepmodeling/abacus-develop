#include "../direct_2d_coulomb.h"

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

double relative_frobenius(const Direct2dCoulomb::DenseMatrix& lhs,
                          const Direct2dCoulomb::DenseMatrix& rhs)
{
    require(lhs.dimension == rhs.dimension, "matrix dimensions differ");
    double numerator = 0.0;
    double denominator = 0.0;
    for (std::size_t i = 0; i < lhs.values.size(); ++i)
    {
        numerator += std::norm(lhs.values[i] - rhs.values[i]);
        denominator += std::norm(rhs.values[i]);
    }
    return std::sqrt(numerator / denominator);
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
        require(std::abs(value.imag()) < 1.0e-11, "quadratic form is not real");
        require(value.real() > -1.0e-11, "quadratic form is negative");
    }
}
} // namespace

int main()
{
    try
    {
        const auto quadrature = Direct2dCoulomb::gauss_legendre_interval(32, -3.0, 3.0);
        double weight_sum = 0.0;
        for (const auto& node: quadrature)
        {
            require(node.weight_bohr > 0.0, "kz quadrature contains a non-positive weight");
            weight_sum += node.weight_bohr;
        }
        require(std::abs(weight_sum - 6.0) < 1.0e-13,
                "kz quadrature does not integrate a constant");

        ModuleBase::Matrix3 reciprocal;
        reciprocal.e11 = 1.0;
        reciprocal.e22 = 1.0;
        reciprocal.e33 = 0.1;
        const ModuleBase::Vector3<double> q_direct(0.125, 0.0, 0.0);
        constexpr double tpiba = 1.0;
        constexpr double ecut_ry = 36.0;
        const auto lines = Direct2dCoulomb::enumerate_reciprocal_lines(
            reciprocal, q_direct, tpiba, ecut_ry);
        require(!lines.empty(), "2D reciprocal-line enumeration is empty");
        for (const auto& line: lines)
        {
            require(std::abs(line.q_plus_g_parallel_cart.z) < 1.0e-15,
                    "2D reciprocal-line enumeration retained a Gz image");
            require(line.kz_max_bohr > 0.0, "2D reciprocal line has no kz interval");
        }

        constexpr std::size_t radial_size = 1201;
        constexpr double dk = 0.01;
        std::vector<double> radial_s(radial_size);
        std::vector<double> radial_p(radial_size);
        for (std::size_t i = 0; i < radial_size; ++i)
        {
            const double k = dk * static_cast<double>(i);
            radial_s[i] = std::exp(-0.3 * k * k);
            radial_p[i] = k * std::exp(-0.25 * k * k);
        }
        const std::vector<Direct2dCoulomb::BasisFunction> basis = {
            {{0, dk, radial_s.data(), radial_s.size()}, 0, 0, {0.0, 0.0, 0.0}},
            {{1, dk, radial_p.data(), radial_p.size()}, 0, 1, {0.17, 0.11, 0.23}},
        };

        constexpr double omega = 100.0;
        constexpr double lz_bohr = 10.0;
        const auto points64 = Direct2dCoulomb::expand_mixed_fourier_points(lines, tpiba, 64);
        const auto points128 = Direct2dCoulomb::expand_mixed_fourier_points(lines, tpiba, 128);
        const auto points256 = Direct2dCoulomb::expand_mixed_fourier_points(lines, tpiba, 256);
        const auto matrix64 = Direct2dCoulomb::build_gram_matrix(basis, points64, omega, lz_bohr, 1.0);
        const auto matrix128 = Direct2dCoulomb::build_gram_matrix(basis, points128, omega, lz_bohr, 1.0);
        const auto matrix256 = Direct2dCoulomb::build_gram_matrix(basis, points256, omega, lz_bohr, 1.0);
        require(matrix128.hermitian_residual < 1.0e-12, "direct Gram matrix is not Hermitian");
        require_positive_quadratic_forms(matrix128.matrix);
        require(relative_frobenius(matrix128.matrix, matrix256.matrix) < 2.0e-10,
                "128-point kz quadrature is not converged");
        require(relative_frobenius(matrix64.matrix, matrix256.matrix) < 2.0e-6,
                "64-point kz quadrature is not within the WS2-tested tolerance");

        const std::array<int, 3> nmp{12, 12, 1};
        constexpr double massidda_chi = 3.25;
        const auto gamma_plane = Direct2dCoulomb::gamma_plane_quadrature(
            reciprocal, tpiba, nmp, 8, massidda_chi);
        double gamma_weight_sum = 0.0;
        double gamma_inverse_q_average = 0.0;
        for (const auto& node: gamma_plane)
        {
            require(node.weight > 0.0, "Gamma-plane quadrature contains a non-positive weight");
            require(node.q_bohr > 0.0, "Gamma-plane quadrature contains q=0");
            gamma_weight_sum += node.weight;
            gamma_inverse_q_average += node.weight / node.q_bohr;
        }
        require(std::abs(gamma_weight_sum - 1.0) < 2.0e-14,
                "Gamma-plane quadrature is not normalized");
        require(std::abs(gamma_inverse_q_average - massidda_chi) < 2.0e-13,
                "Gamma-plane quadrature does not reproduce the Massidda average");

        Direct2dCoulomb::MethodMetadata metadata;
        metadata.method = "direct_mixed_fourier";
        metadata.ecut_ry = 110.0;
        metadata.kz_order = 64;
        metadata.gamma_order = 8;
        metadata.nq = 144;
        metadata.naux = 807;
        metadata.source_revision = "test-revision";
        const std::string formatted = Direct2dCoulomb::format_method_metadata(metadata);
        const std::string expected_metadata =
            "# ABACUS reader-v1 strict 2D Coulomb method\n"
            "version = 1\n"
            "method = direct_mixed_fourier\n"
            "ecut_ry = 110\n"
            "kz_order = 64\n"
            "gamma_order = 8\n"
            "nq = 144\n"
            "naux = 807\n"
            "source_revision = test-revision\n";
        require(formatted == expected_metadata, "method metadata is not deterministic");
        require(formatted.find("method = direct_mixed_fourier") != std::string::npos,
                "metadata omits method");
        require(formatted.find("ecut_ry = 110") != std::string::npos,
                "metadata omits cutoff");
        require(formatted.find("kz_order = 64") != std::string::npos,
                "metadata omits kz order");
        require(formatted.find("gamma_order = 8") != std::string::npos,
                "metadata omits Gamma order");
        require(formatted.find("nq = 144") != std::string::npos,
                "metadata omits q count");
        require(formatted.find("naux = 807") != std::string::npos,
                "metadata omits auxiliary dimension");

        std::cout << "direct strict-2D Coulomb tests passed" << std::endl;
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << std::endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
