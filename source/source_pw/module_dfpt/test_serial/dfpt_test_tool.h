#ifndef DFPT_TEST_TOOL_H
#define DFPT_TEST_TOOL_H

// Shared free-function utilities for the DFPT serial unit tests
// (module_pw test/test_tool.{h,cpp} pattern): every test file derives
// its own gtest fixture from ::testing::Test, owns the objects it needs
// and calls these functions -- there is no shared fixture base class.

#include <complex>

class UnitCell;
namespace ModuleBase
{
class Matrix3;
template <class T>
class Vector3;
}
namespace ModulePW
{
class PW_Basis;
class PW_Basis_K;
}
namespace ModuleCell
{
class QList;
}
namespace ModuleDFPT
{
class DFPT_PW_Data;
}

namespace dfpt_test
{

// shared test-world conventions (the fixture values every test agrees on)
constexpr double lat0 = 1.8897261254578281; // bohr
constexpr double ecutwfc = 2.5;             // Ry
// rho cutoff inflated to 9x ecutwfc so every Delta = G''-G' of the
// convolution lies inside the rho ball and nothing aliases
constexpr double rho_mult = 9.0;
constexpr double cubic_edge = 10.0; // cubic cell edge in lat0 units

// Build the cubic single-atom cell (latvec / G / GT / tpiba / omega and
// the iat2it / iat2ia wiring; the atom itself is Coulomb flavor until a
// make_*_atom call replaces it).
void setup_cell(UnitCell& ucell, ModuleBase::Matrix3& latvec);

// Reconfigure the cell's single atom as a bare Coulomb Z=4 atom at tau.
void make_coulomb_atom(UnitCell& ucell);

// Reconfigure the cell's single atom as a norm-conserving pseudo atom
// (analytic Bessel projectors, no nlcc).
void make_nc_atom(UnitCell& ucell);

// Reconfigure the cell as a two-atom Z=4/Z=2 crystal breaking all
// symmetry (masses 12/4), the second atom at the anti-blind-spot tau
// (0.25, 0.31, 0.17).
void make_two_atom_cell(UnitCell& ucell);

// Extend the two-atom cell with a nonlinear core correction on the
// second species: a Gaussian rho_atc of total charge 0.24 e at the
// non-origin tau (0.25, 0.31, 0.17).
void make_nlcc_two_atom_cell(UnitCell& ucell);

// Assemble the shared rho/wfc bases for a (k, q) pair and bind the q
// list and the data object: the rho basis at rho_mult * ecutwfc, the
// wfc basis on the rho FFT box with the k list, the q list holding the
// single q, and q_cart = q_frac * ucell.G. Convention: k defaults to -q
// so the k+q ball stays inside the ground-state G list.
void setup_bases(ModulePW::PW_Basis& pw_rho,
                 ModulePW::PW_Basis_K& pw_wfc,
                 ModuleCell::QList& qlist,
                 ModuleDFPT::DFPT_PW_Data& data,
                 const UnitCell& ucell,
                 const ModuleBase::Vector3<double>& k_d,
                 const ModuleBase::Vector3<double>& q_d,
                 int nbands,
                 ModuleBase::Vector3<double>& q_cart);

// Analytic Coulomb local potential (Ry) at g2_bohr (bohr^-2) for valence
// zv in cell volume omega, mirroring
// vl_pw.cpp::vloc_coulomb independently of DFPT_Pert::vloc_at_g.
double vloc_coulomb(double zv, double omega, double g2_bohr);

// Analytic dVloc/dtau_alpha coefficient at displacement vector w
// (1/lat0): dV/dtau = -i (w)_alpha tpiba Vloc exp(-i 2pi w.tau) with the
// GS structure-factor phase convention.
std::complex<double> analytic_dvloc(const UnitCell& ucell,
                                    const ModuleBase::Vector3<double>& tau,
                                    int dir,
                                    const ModuleBase::Vector3<double>& w);

// Integer FFT triple key (ix, iy, iz) on the cubic cell.
long long fkey(int ix, int iy, int iz);
// Key of g on the cubic cell (gcar * a is integral).
long long gkey(const ModuleBase::Vector3<double>& g, double a);

} // namespace dfpt_test

#endif // DFPT_TEST_TOOL_H
