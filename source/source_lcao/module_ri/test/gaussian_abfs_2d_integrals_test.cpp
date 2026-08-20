#include "../gaussian_abfs_2d_integrals.h"

#include <cmath>
#include <cstdlib>
#include <iostream>

namespace
{
void require(const bool condition, const char* message)
{
    if (!condition)
    {
        std::cerr << message << '\n';
        std::exit(EXIT_FAILURE);
    }
}
} // namespace

int main()
{
    const double coefficient = std::sqrt(3.0 / (4.0 * GaussianAbfs2D::pi));

    const double q_in_plane = 0.06029258;
    const double expected_in_plane
        = -coefficient * GaussianAbfs2D::pi * std::exp(q_in_plane * q_in_plane)
          * std::erfc(q_in_plane);
    require(std::abs(GaussianAbfs2D::coulomb_l1(1, q_in_plane, 0.0, 0.0, 1.0).real()
                     - expected_in_plane)
                < 1.0e-13,
            "in-plane x L=1 analytic Coulomb integral mismatch");
    require(std::abs(GaussianAbfs2D::coulomb_l1(2, 0.0, q_in_plane, 0.0, 1.0).real()
                     - expected_in_plane)
                < 1.0e-13,
            "in-plane y L=1 analytic Coulomb integral mismatch");

    const double q_out_of_plane = 0.10048763;
    const double tau = 0.7;
    const double step = 1.0e-5;
    const double finite_difference
        = (GaussianAbfs2D::i2(q_out_of_plane, tau + step, 1.0)
           - GaussianAbfs2D::i2(q_out_of_plane, tau - step, 1.0))
          / (2.0 * step);
    require(std::abs(GaussianAbfs2D::di2_dtau(q_out_of_plane, tau, 1.0)
                     - finite_difference)
                < 2.0e-9,
            "out-of-plane I2 derivative mismatch");
    require(std::abs(GaussianAbfs2D::coulomb_l1(0, q_out_of_plane, 0.0, tau, 1.0).imag()
                     + coefficient * finite_difference)
                < 2.0e-9,
            "out-of-plane L=1 analytic Coulomb integral mismatch");

    const double small_q_value
        = GaussianAbfs2D::coulomb_l1(1, 1.0e-8, 0.0, 0.0, 1.0).real();
    require(std::abs(small_q_value + coefficient * GaussianAbfs2D::pi) < 2.0e-8,
            "in-plane L=1 small-q limit is not finite");

    return EXIT_SUCCESS;
}
