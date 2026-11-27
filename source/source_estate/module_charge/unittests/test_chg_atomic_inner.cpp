#include "gtest/gtest.h"

#include "source_base/math_integral.h"
#include "source_base/matrix3.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_cell/unitcell.h"
#include "source_cell/atom_spec.h"
#include "source_estate/module_charge/charge.h"
#include "source_estate/module_charge/chg_atomic_detail.h"

#include <cmath>
#include <complex>
#include <sstream>
#include <vector>

// charge.cpp references Magnetism; provide a lightweight stub.
Magnetism::Magnetism()
{
    this->tot_mag = 0.0;
    this->abs_mag = 0.0;
}
Magnetism::~Magnetism()
{
}

/************************************************
 *  unit test of module_charge/chg_atomic_inner.cpp
 ***********************************************/

/**
 * - Tested Functions:
 *   - compute_rhoatm: builds the atomic radial charge density for USPP
 *     (tvanp=true: direct copy of rho_at) and NCPP (tvanp=false: divide by
 *     4*pi*r^2, extrapolate rho[0], Simpson-integrate, scale to zv, remultiply
 *     4*pi*r^2).
 *   - normalize_and_check: FFTs rho_g3d to real space, sums electrons,
 *     renormalizes rho to the given nelec.
 *
 * For compute_rhoatm we use a Gaussian rho_at whose analytic integral is known.
 */

namespace
{

/// Build a uniform radial mesh r[ir] = dr * ir, rab = dr.
void fill_uniform_mesh(Atom_pseudo& ncpp, int mesh, double dr, double zv, bool tvanp)
{
    ncpp.mesh = mesh;
    ncpp.msh = mesh;
    ncpp.zv = zv;
    ncpp.tvanp = tvanp;
    ncpp.r.assign(mesh, 0.0);
    ncpp.rab.assign(mesh, dr);
    ncpp.rho_at.assign(mesh, 0.0);
    for (int ir = 0; ir < mesh; ++ir)
    {
        ncpp.r[ir] = dr * ir;
    }
}

/// rho_at(r) = 4 pi r^2 * Gaussian, so the number density is a pure Gaussian.
void fill_gaussian_rho_at(Atom_pseudo& ncpp, double alpha, double norm)
{
    for (int ir = 0; ir < ncpp.mesh; ++ir)
    {
        const double r = ncpp.r[ir];
        ncpp.rho_at[ir] = norm * ModuleBase::FOUR_PI * r * r * std::exp(-alpha * r * r);
    }
}

} // namespace

TEST(ChgAtomicInnerTest, ComputeRhoatmUsppCopiesRhoAt)
{
    Atom atom;
    fill_uniform_mesh(atom.ncpp, 8, 0.5, 8.0, true);
    for (int ir = 0; ir < 8; ++ir)
    {
        atom.ncpp.rho_at[ir] = static_cast<double>(ir + 1);
    }

    std::stringstream ofs;
    const std::vector<double> rhoatm = module_charge::detail::compute_rhoatm(atom, 8, ofs);

    ASSERT_EQ(rhoatm.size(), 8u);
    for (int ir = 0; ir < 8; ++ir)
    {
        EXPECT_EQ(rhoatm[ir], atom.ncpp.rho_at[ir]);
    }
}

TEST(ChgAtomicInnerTest, ComputeRhoatmNcppIntegratesAndScalesToZv)
{
    Atom atom;
    const int mesh = 101;
    const double dr = 0.1;
    const double zv = 4.0;
    fill_uniform_mesh(atom.ncpp, mesh, dr, zv, false);
    // rho_at = 4 pi r^2 * exp(-r^2), integral over [0,inf) = pi^{3/2}.
    fill_gaussian_rho_at(atom.ncpp, 1.0, 1.0);

    std::stringstream ofs;
    const std::vector<double> rhoatm = module_charge::detail::compute_rhoatm(atom, mesh, ofs);

    // for NCPP, rhoatm[ir] = rho_at[ir] (after scaling) for ir>0 because the
    // /4pir^2 and *4pir^2 cancel; the net effect is scale = zv / charge.
    double charge = 0.0;
    ModuleBase::Integral::Simpson_Integral(atom.ncpp.msh,
                                           atom.ncpp.rho_at.data(),
                                           atom.ncpp.rab.data(),
                                           charge);
    const double scale = zv / charge;
    for (int ir = 1; ir < mesh; ++ir)
    {
        EXPECT_NEAR(rhoatm[ir], atom.ncpp.rho_at[ir] * scale, 1e-8);
    }
}

TEST(ChgAtomicInnerTest, NormalizeAndCheckRenormalizesToNelec)
{
    ModulePW::PW_Basis rhopw;
    rhopw.initgrids(1.0, ModuleBase::Matrix3(1, 0, 0, 0, 1, 0, 0, 0, 1), 10.0);
    rhopw.initparameters(false, 10.0);
    rhopw.setuptransform();
    rhopw.collect_local_pw();

    const int spin_number_need = 1;
    const double omega = rhopw.omega;
    const double nelec = 5.0;
    ModuleBase::ComplexMatrix rho_g3d(spin_number_need, rhopw.npw);
    // put all weight at G=0 => uniform real-space density.
    rho_g3d(0, 0) = std::complex<double>(1.0, 0.0);

    std::vector<double> rho_in(rhopw.nrxx, 0.0);
    double* rho_ptrs[1] = {rho_in.data()};

    std::stringstream ofs;
    module_charge::detail::normalize_and_check(rho_ptrs, rho_g3d, &rhopw,
                                               spin_number_need, omega, ofs, nelec);

    double ne = 0.0;
    for (int ir = 0; ir < rhopw.nrxx; ++ir)
    {
        ne += rho_in[ir];
    }
    ne *= omega / static_cast<double>(rhopw.nxyz);
    EXPECT_NEAR(ne, nelec, 1e-6);
}
