#include "../gaussian_abfs.h"
#include "source_cell/unitcell.h"

#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <vector>

namespace
{
void require(const bool condition, const char* message)
{
    if (!condition)
    {
        std::cerr << "FAILED: " << message << '\n';
        std::exit(1);
    }
}
} // namespace

// Gamma-point retention of the L=1,M=0 normal-dipole channel in the 2D
// Ewald lattice sum. That channel couples an s auxiliary to a pz auxiliary
// through the Gaunt coefficient <Y_00 Y_1z|Y_1z>: the Vq element (s,pz) is a
// prefactor times lattice_sum[0][lm=1]. At q=0 the channel has the unique
// finite limit -pi*erf(tau_z/(2*sqrt(beta))) (direction-independent, unlike
// the dropped L=1 in-plane channels), so the excluded Gamma term must keep
// it instead of discarding it together with the singular L=0 channel.
int main()
{
    // Square in-plane lattice with |a| = 1 bohr: the reciprocal cutoff
    // sqrt(35*lambda) ~ 5.92 bohr^-1 stays below the first nonzero
    // reciprocal shell 2*pi/|a| ~ 6.28 bohr^-1, so the 2D lattice sum
    // contains only the q=0 (Gamma) sample and the retained channel can be
    // asserted against its analytic limit without any other contribution.
    UnitCell ucell;
    ucell.lat0 = 1.0;
    ucell.tpiba = ModuleBase::TWO_PI / ucell.lat0;
    ucell.a1 = ModuleBase::Vector3<double>(1.0, 0.0, 0.0);
    ucell.a2 = ModuleBase::Vector3<double>(0.0, 1.0, 0.0);
    ucell.a3 = ModuleBase::Vector3<double>(0.0, 0.0, 3.0);
    ucell.omega = 3.0;

    ModuleBase::Matrix3 G;
    G.e11 = 1.0;
    G.e12 = 0.0;
    G.e13 = 0.0;
    G.e21 = 0.0;
    G.e22 = 1.0;
    G.e23 = 0.0;
    G.e31 = 0.0;
    G.e32 = 0.0;
    G.e33 = 1.0 / 3.0;

    const double lambda = 1.0;
    const double beta = 1.0 / lambda;
    Gaussian_Abfs gaussian_abfs;
    const std::vector<ModuleBase::Vector3<double>> kvec_c
        = {ModuleBase::Vector3<double>(0.0, 0.0, 0.0)};
    gaussian_abfs.init(ucell, 1, kvec_c, G, lambda, 2);

    // lz_bohr = omega / |a1_bohr x a2_bohr| = 3.0 / 1.0.
    const double conversion = 3.0 / ModuleBase::TWO_PI;
    const double coefficient = std::sqrt(3.0 / (4.0 * ModuleBase::PI));

    for (const double tau_z_alat : {0.6, -0.6})
    {
        const ModuleBase::Vector3<double> tau(0.13, 0.29, tau_z_alat);
        const double tau_z_bohr = tau.z * ucell.lat0;
        const auto sum = gaussian_abfs.get_lattice_sum_2d(ucell.tpiba, 0, -2.0, beta, true, 1, tau);

        require(sum.size() == 4, "lmax=1 lattice sum must carry four channels");

        // The singular L=0 channel is dropped here; DPcal_Vq_2d restores it
        // analytically together with the chi2d correction.
        require(std::abs(sum[0]) == 0.0, "the singular L=0 Gamma channel must be dropped");

        // The in-plane L=1 channels are direction-ambiguous at q->0 and dropped.
        require(std::abs(sum[2]) == 0.0 && std::abs(sum[3]) == 0.0,
                "the direction-ambiguous in-plane L=1 Gamma channels must be dropped");

        // The L=1,M=0 (pz) channel retains the analytic normal-dipole limit.
        const std::complex<double> expected_dipole(
            0.0,
            conversion * coefficient * ModuleBase::PI * std::erf(tau_z_bohr / (2.0 * std::sqrt(beta))));
        require(std::abs(sum[1] - expected_dipole)
                    < 1.0e-13 * std::max(1.0, std::abs(expected_dipole)),
                "the Gamma s-pz normal-dipole limit must be retained with the analytic value");
    }

    // Co-planar centers: the dipole limit vanishes, so no spurious coupling.
    const ModuleBase::Vector3<double> tau_coplanar(0.13, 0.29, 0.0);
    const auto sum_coplanar
        = gaussian_abfs.get_lattice_sum_2d(ucell.tpiba, 0, -2.0, beta, true, 1, tau_coplanar);
    require(std::abs(sum_coplanar[1]) == 0.0,
            "the Gamma s-pz normal-dipole limit must vanish at tau_z = 0");

    return 0;
}
