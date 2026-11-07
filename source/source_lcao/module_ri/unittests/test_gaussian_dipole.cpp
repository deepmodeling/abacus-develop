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

    // Ewald_Vq::init builds the radial Gaussian tables through Gauss() before
    // Gaussian_Abfs::init runs, so the dimension must be settable beforehand.
    {
        const double strict_cutoff = Gaussian_Abfs::radial_cutoff(lambda, 8, 2);
        const double historical_cutoff = Gaussian_Abfs::radial_cutoff(lambda, 8, 3);
        require(strict_cutoff > historical_cutoff, "strict-2D cutoff must exceed the historical 3D cutoff");

        auto probe_orbital = [](const double dr) {
            const int nr = 101;
            std::vector<double> rab(nr, dr);
            std::vector<double> r(nr, 0.0);
            std::vector<double> radial(nr, 0.0);
            for (int ir = 0; ir != nr; ++ir)
            {
                r[ir] = ir * dr;
            }
            Numerical_Orbital_Lm orbital;
            orbital.set_orbital_info("N",
                                     0,
                                     8,
                                     0,
                                     nr,
                                     rab.data(),
                                     r.data(),
                                     Numerical_Orbital_Lm::Psi_Type::Psi,
                                     radial.data(),
                                     101,
                                     0.05,
                                     dr,
                                     false,
                                     true,
                                     false);
            return orbital;
        };

        Gaussian_Abfs strict_gaussian;
        strict_gaussian.set_ewald_dimension(2);

        // Even ceil(rcut/dr): the odd point-count bump puts the endpoint past
        // the selected cutoff.
        const double bumped_dr = strict_cutoff / 99.5;
        const Numerical_Orbital_Lm bumped_table = strict_gaussian.Gauss(probe_orbital(bumped_dr), lambda);
        require(bumped_table.getNr() % 2 == 1, "the Gaussian radial grid must keep an odd point count");
        require(std::abs(bumped_table.getRcut() - (bumped_table.getNr() - 1) * bumped_dr)
                    <= 1.0e-12 * bumped_dr,
                "the Gaussian rcut must be the grid endpoint (Nr-1)*dr");
        require(bumped_table.getRcut() >= strict_cutoff,
                "a dimension set before Gauss() must select the strict-2D radial extent");

        // Odd ceil(rcut/dr): the inherited grid stops one step below the
        // cutoff, so only the endpoint-plus-step bound holds.
        const double odd_dr = strict_cutoff / 98.5;
        const Numerical_Orbital_Lm odd_table = strict_gaussian.Gauss(probe_orbital(odd_dr), lambda);
        require(odd_table.getNr() % 2 == 1, "the odd ceil(rcut/dr) grid must keep its odd point count");
        require(odd_table.getRcut() < strict_cutoff,
                "an odd ceil(rcut/dr) grid ends below the selected cutoff");
        require(odd_table.getRcut() + odd_dr >= strict_cutoff,
                "the odd ceil(rcut/dr) grid must stay within one step of the cutoff");

        Gaussian_Abfs historical_gaussian;
        const Numerical_Orbital_Lm historical_table
            = historical_gaussian.Gauss(probe_orbital(bumped_dr), lambda);
        require(historical_table.getNr() < bumped_table.getNr(),
                "the default 3D Gauss() must keep the historical radial extent");
    }

    return 0;
}
