// Shared free-function utilities for the DFPT serial unit tests,
// implementation (see dfpt_test_tool.h). Everything runs without __MPI:
// the plane-wave bases are built through the real serial
// initgrids/initparameters/setuptransform path on a shared FFT grid,
// exactly like the production setup_pwrho/setup_pwwfc sequence.

#include "dfpt_test_tool.h"

#include "source_base/constants.h"
#include "source_base/matrix3.h"
#include "source_base/vector3.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_basis/module_pw/pw_basis_k.h"
#include "source_cell/atom_pseudo.h"
#include "source_cell/atom_spec.h"
#include "source_cell/qlist.h"
#include "source_cell/unitcell.h"
#include "source_psi/psi.h"
#include "source_pw/module_dfpt/dfpt_pw_data.h"

#include <cmath>
#include <vector>

namespace dfpt_test
{

// nlcc Gaussian conventions of make_nlcc_two_atom_cell
constexpr double nlcc_sigma = 1.5;        // bohr, grid-resolvable (spacing ~0.63 bohr)
constexpr double nlcc_core_charge = 0.24; // e, integrated

void setup_cell(UnitCell& ucell, ModuleBase::Matrix3& latvec)
{
    latvec = ModuleBase::Matrix3(cubic_edge, 0.0, 0.0, 0.0, cubic_edge, 0.0, 0.0, 0.0, cubic_edge);
    ucell.ntype = 1;
    ucell.nat = 1;
    ucell.atoms = new Atom[1];
    ucell.atoms[0].na = 1;
    ucell.atoms[0].tau.resize(1);
    ucell.atoms[0].tau[0] = ModuleBase::Vector3<double>(1.1, 2.3, 0.7);
    ucell.latvec = latvec;
    ucell.GT = latvec.Inverse();
    ucell.G = ucell.GT.Transpose();
    ucell.lat0 = lat0;
    ucell.tpiba = ModuleBase::TWO_PI / lat0;
    ucell.tpiba2 = ucell.tpiba * ucell.tpiba;
    ucell.omega = cubic_edge * cubic_edge * cubic_edge * lat0 * lat0 * lat0;
    ucell.iat2it.resize(1);
    ucell.iat2ia.resize(1);
    ucell.iat2it[0] = 0;
    ucell.iat2ia[0] = 0;
    make_coulomb_atom(ucell);
}

void make_coulomb_atom(UnitCell& ucell)
{
    Atom& at = ucell.atoms[0];
    at.label = "C";
    at.coulomb_potential = true;
    at.ncpp.zv = 4.0;
    at.ncpp.tvanp = false;
    at.ncpp.has_so = false;
    at.ncpp.nbeta = 0;
    at.ncpp.nh = 0;
    at.ncpp.msh = 0;
    at.ncpp.kkbeta = 0;
    at.mass = 12.0;
}

void make_nc_atom(UnitCell& ucell)
{
    Atom& at = ucell.atoms[0];
    at.label = "Si";
    at.coulomb_potential = false;
    pseudo& p = at.ncpp;
    p.zv = 4.0;
    p.tvanp = false;
    p.has_so = false;
    p.nbeta = 2;
    p.lll = {0, 1};
    p.nh = 4;
    p.msh = 121;
    p.kkbeta = 121;
    p.r.resize(121);
    p.rab.resize(121);
    p.vloc_at.assign(121, 0.0);
    const double dx = 0.025;
    for (int i = 0; i < 121; ++i)
    {
        p.r[i] = i * dx;
        p.rab[i] = dx;
    }
    p.betar.create(2, 121);
    for (int i = 0; i < 121; ++i)
    {
        const double r = p.r[i];
        p.betar(0, i) = std::exp(-std::pow(r - 1.0, 2) / (2.0 * 0.3 * 0.3));
        p.betar(1, i) = std::exp(-std::pow(r - 1.2, 2) / (2.0 * 0.35 * 0.35));
    }
    p.dion.create(2, 2);
    p.dion(0, 0) = 0.8;
    p.dion(0, 1) = 0.15;
    p.dion(1, 0) = -0.25;
    p.dion(1, 1) = 1.1;
}

void make_two_atom_cell(UnitCell& ucell)
{
    ucell.ntype = 2;
    ucell.nat = 2;
    delete[] ucell.atoms;
    ucell.atoms = new Atom[2];
    ucell.atoms[0].na = 1;
    ucell.atoms[1].na = 1;
    ucell.atoms[0].tau.resize(1);
    ucell.atoms[1].tau.resize(1);
    ucell.atoms[0].tau[0] = ModuleBase::Vector3<double>(0.0, 0.0, 0.0);
    ucell.atoms[1].tau[0] = ModuleBase::Vector3<double>(0.25, 0.31, 0.17);
    for (int it = 0; it < 2; ++it)
    {
        Atom& at = ucell.atoms[it];
        at.label = (it == 0) ? "A" : "B";
        at.coulomb_potential = true;
        at.ncpp.zv = (it == 0) ? 4.0 : 2.0;
        at.ncpp.tvanp = false;
        at.ncpp.has_so = false;
        at.ncpp.nbeta = 0;
        at.ncpp.nh = 0;
        at.mass = (it == 0) ? 12.0 : 4.0;
    }
    ucell.iat2it.resize(2);
    ucell.iat2ia.resize(2);
    ucell.iat2it[0] = 0;
    ucell.iat2ia[0] = 0;
    ucell.iat2it[1] = 1;
    ucell.iat2ia[1] = 0;
}

void make_nlcc_two_atom_cell(UnitCell& ucell)
{
    // start from the symmetry-broken two-atom cell, then upgrade the
    // second species (the one at the non-origin tau) to an NC pseudo
    // with a nonlinear core correction
    make_two_atom_cell(ucell);
    Atom& at = ucell.atoms[1];
    at.label = "Bc";
    at.coulomb_potential = false;
    pseudo& p = at.ncpp;
    p.zv = 2.0;
    p.tvanp = false;
    p.has_so = false;
    p.nbeta = 0;
    p.lll = {};
    p.nh = 0;
    p.msh = 121;
    p.kkbeta = 121;
    p.r.resize(121);
    p.rab.resize(121);
    p.vloc_at.assign(121, 0.0);
    const double dx = 0.025;
    for (int i = 0; i < 121; ++i)
    {
        p.r[i] = i * dx;
        p.rab[i] = dx;
    }
    // Gaussian core charge normalized to the integrated nlcc_core_charge
    // (rho_atc carries the 4 pi r^2 measure of the radial grid)
    const double s2 = nlcc_sigma * nlcc_sigma;
    double norm = 0.0;
    for (int i = 0; i < 121; ++i)
    {
        const double r = p.r[i];
        norm += r * r * std::exp(-r * r / (2.0 * s2)) * p.rab[i];
    }
    const double amp = nlcc_core_charge / (norm * ModuleBase::FOUR_PI);
    p.rho_atc.resize(121);
    for (int i = 0; i < 121; ++i)
    {
        const double r = p.r[i];
        p.rho_atc[i] = amp * std::exp(-r * r / (2.0 * s2));
    }
    p.nlcc = true;
}

void setup_bases(ModulePW::PW_Basis& pw_rho,
                 ModulePW::PW_Basis_K& pw_wfc,
                 ModuleCell::QList& qlist,
                 ModuleDFPT::DFPT_PW_Data& data,
                 const UnitCell& ucell,
                 const ModuleBase::Vector3<double>& k_d,
                 const ModuleBase::Vector3<double>& q_d,
                 int nbands,
                 ModuleBase::Vector3<double>& q_cart)
{
    // The tool keeps its own five-call assembly chains instead of
    // calling the production pw::setup_pwrho / setup_pwwfc wrappers
    // (which also run without __MPI): the tests need rho_mult = 9x
    // ecutwfc for the DFPT convolution (every Delta = G''-G' must fit
    // in the rho ball) while production hardcodes 4x ecutwfc; they
    // re-assemble for arbitrary (k, q) pairs per case, which the
    // K_Vectors-driven wrapper is not shaped for; and the wrappers
    // read GlobalV::ofs_running and require a fabricated Input_para.
    const ModuleBase::Matrix3& latvec = ucell.latvec;
    pw_rho.initgrids(lat0, latvec, rho_mult * ecutwfc);
    pw_rho.initparameters(false, rho_mult * ecutwfc);
    pw_rho.fft_bundle.initfftmode(0);
    pw_rho.setuptransform();
    pw_rho.collect_local_pw();

    const ModuleBase::Vector3<double> klist[1] = {k_d};
    pw_wfc.initgrids(lat0, latvec, pw_rho.nx, pw_rho.ny, pw_rho.nz);
    pw_wfc.initparameters(false, ecutwfc, 1, klist);
    pw_wfc.fft_bundle.initfftmode(0);
    pw_wfc.setuptransform();
    pw_wfc.collect_local_pw();

    qlist.nkstot = 1;
    qlist.kvec_d.clear();
    qlist.kvec_d.push_back(q_d);
    q_cart = q_d * ucell.G;

    data.init(&qlist, 1, nbands, pw_wfc.npwk_max, pw_rho.nrxx, 1, 1, nullptr);
}

double vloc_coulomb(double zv, double omega, double g2_bohr)
{
    return -zv * ModuleBase::e2 * ModuleBase::FOUR_PI / omega / g2_bohr;
}

std::complex<double> analytic_dvloc(const UnitCell& ucell,
                                    const ModuleBase::Vector3<double>& tau,
                                    int dir,
                                    const ModuleBase::Vector3<double>& w)
{
    const double w2 = w * w;
    if (w2 < 1.0e-12)
    {
        return std::complex<double>(0.0, 0.0);
    }
    const double zv = ucell.atoms[0].ncpp.zv;
    const double arg = -ModuleBase::TWO_PI * (w * tau);
    return std::complex<double>(0.0, -1.0) * (ucell.tpiba * w[dir])
           * vloc_coulomb(zv, ucell.omega, w2 * ucell.tpiba2)
           * std::complex<double>(std::cos(arg), std::sin(arg));
}

long long fkey(int ix, int iy, int iz)
{
    return (static_cast<long long>(ix + 64) * 128 + (iy + 64)) * 128 + (iz + 64);
}

long long gkey(const ModuleBase::Vector3<double>& g, double a)
{
    return fkey(static_cast<int>(std::llround(g.x * a)),
                static_cast<int>(std::llround(g.y * a)),
                static_cast<int>(std::llround(g.z * a)));
}

} // namespace dfpt_test
