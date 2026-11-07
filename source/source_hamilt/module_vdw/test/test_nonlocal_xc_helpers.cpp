#include "test_nonlocal_xc.h"

namespace rvv10_test
{

double valence_density(const ModulePW::PW_Basis& pw, const int ix, const int iy, const int iz)
{
    const double ax = ModuleBase::TWO_PI * ix / pw.nx;
    const double ay = ModuleBase::TWO_PI * iy / pw.ny;
    const double az = ModuleBase::TWO_PI * iz / pw.nz;
    const double ayz = ay + az;
    return 0.004 + 0.006 * (1.0 + 0.25 * std::sin(ax) + 0.15 * std::cos(ayz));
}

double core_density(const ModulePW::PW_Basis& pw, const int ix, const int iy, const int iz)
{
    const double az = ModuleBase::TWO_PI * iz / pw.nz;
    return 0.002 * (1.0 + std::cos(az));
}

double kinetic_density(const ModulePW::PW_Basis& pw, const int ix, const int iy, const int iz)
{
    const double ay = ModuleBase::TWO_PI * iy / pw.ny;
    return 0.003 * (1.0 + 0.2 * std::cos(ay));
}

double perturbation(const ModulePW::PW_Basis& pw, const int ix, const int iy, const int iz)
{
    const double ax = ModuleBase::TWO_PI * ix / pw.nx;
    const double ay = ModuleBase::TWO_PI * iy / pw.ny;
    const double axy = ax + ay;
    const double ay2 = 2.0 * ay;
    return std::sin(axy) + 0.3 * std::cos(ay2);
}

std::vector<double> local_grid(const ModulePW::PW_Basis& pw, GridFunction f)
{
    std::vector<double> values(pw.nrxx);
    for (int ix = 0; ix < pw.nx; ++ix)
    {
        for (int iy = 0; iy < pw.ny; ++iy)
        {
            for (int iz = 0; iz < pw.nplane; ++iz)
            {
                values[(ix * pw.ny + iy) * pw.nplane + iz] = f(pw, ix, iy, pw.startz_current + iz);
            }
        }
    }
    return values;
}

ModuleBase::Matrix3 oblique_lattice(const double z_sign)
{
    const double c3 = 2.0 * z_sign;
    return ModuleBase::Matrix3(1.0, 0.0, 0.0, -0.5, std::sqrt(0.75), 0.0, 0.0, 0.0, c3);
}

ModuleBase::Matrix3 strained_lattice(const ModuleBase::Matrix3& lattice,
                                     const int a,
                                     const int b,
                                     const double h)
{
    double m[3][3] = {{lattice.e11, lattice.e12, lattice.e13},
                      {lattice.e21, lattice.e22, lattice.e23},
                      {lattice.e31, lattice.e32, lattice.e33}};
    for (int i = 0; i < 3; ++i)
    {
        m[i][a] += h * m[i][b];
    }
    return ModuleBase::Matrix3(m[0][0], m[0][1], m[0][2], m[1][0], m[1][1], m[1][2], m[2][0], m[2][1], m[2][2]);
}

void setup_grid(ModulePW::PW_Basis& basis,
                const ModuleBase::Matrix3& lattice,
                const double cutoff,
                const int nx,
                const int ny,
                const int nz)
{
    basis.initgrids(10.0, lattice, nx, ny, nz);
    basis.initparameters(false, cutoff, 1, true);
    basis.setuptransform();
    basis.collect_local_pw();
}

void setup(ModulePW::PW_Basis& basis, const ModuleBase::Matrix3& lattice, const double cutoff)
{
    setup_grid(basis, lattice, cutoff, 5, 7, 9);
}

XCFunctionalParameters rvv10_parameters(const bool rvv10)
{
    XCFunctionalParameters parameters = XC_Functional::get_runtime_parameters();
    parameters.rvv10 = rvv10;
    parameters.rvv10_b = rvv10_b_rpw86;
    parameters.rvv10_c = rvv10_c_rpw86;
    return parameters;
}

VxcResult rvv10_result(ModulePW::PW_Basis& pw,
                       const std::vector<std::vector<double>>& channels,
                       const std::vector<double>& core)
{
    const ChargeOnGrid charge(pw, channels, core);
    const XCFunctionalParameters parameters = rvv10_parameters(true);
    double etxc = 0.0;
    double vtxc = 0.0;
    ModuleBase::matrix v(charge.chr.nspin, pw.nrxx);
    vdw::add_rvv10(parameters, charge.chr, etxc, vtxc, v);
    return VxcResult(etxc, vtxc, v);
}

ModuleBase::matrix rvv10_stress(ModulePW::PW_Basis& pw,
                                const std::vector<std::vector<double>>& channels,
                                const std::vector<double>& core)
{
    const ChargeOnGrid charge(pw, channels, core);
    const XCFunctionalParameters parameters = rvv10_parameters(true);
    ModuleBase::matrix sigma(3, 3);
    vdw::add_rvv10_stress(parameters, charge.chr, sigma);
    return sigma;
}

void expect_same_result(const VxcResult& actual,
                        const VxcResult& expected,
                        const double scalar_tolerance,
                        const double point_tolerance)
{
    EXPECT_NEAR(std::get<0>(actual), std::get<0>(expected), scalar_tolerance);
    EXPECT_NEAR(std::get<1>(actual), std::get<1>(expected), scalar_tolerance);
    const ModuleBase::matrix& v = std::get<2>(actual);
    const ModuleBase::matrix& w = std::get<2>(expected);
    ASSERT_EQ(v.nr, w.nr);
    ASSERT_EQ(v.nc, w.nc);
    for (int is = 0; is < v.nr; ++is)
    {
        for (int ir = 0; ir < v.nc; ++ir)
        {
            EXPECT_NEAR(v(is, ir), w(is, ir), point_tolerance);
        }
    }
}

void expect_rvv10_added(const VxcResult& total,
                        VxcResult semilocal,
                        const Charge& chr,
                        const double point_tolerance)
{
    const XCFunctionalParameters parameters = rvv10_parameters(true);
    double& etxc = std::get<0>(semilocal);
    double& vtxc = std::get<1>(semilocal);
    ModuleBase::matrix& v = std::get<2>(semilocal);
    vdw::add_rvv10(parameters, chr, etxc, vtxc, v);
    expect_same_result(total, semilocal, rounding_tolerance, point_tolerance);
}

void expect_vxc_adds_rvv10(const char* functional,
                           const Charge& chr,
                           const UnitCell& ucell,
                           const bool domag,
                           const int gga_grad,
                           const double point_tolerance)
{
    XC_Functional::ScopedType scope(functional);
    const int nrxx = chr.nrxx;
    const int nspin = chr.nspin;
    const bool domag_z = false;
    const double hybrid_alpha = 0.0;
    const double hse_omega = 0.0;
    VxcResult total;
    VxcResult semilocal;
    for (int with_rvv10 = 0; with_rvv10 < 2; ++with_rvv10)
    {
        const bool rvv10 = with_rvv10 == 1;
        RuntimeParameters parameters(rvv10);
        VxcResult& result = rvv10 ? total : semilocal;
        result = XC_Functional::v_xc(nrxx, &chr, &ucell, nspin, domag, domag_z, gga_grad, hybrid_alpha, hse_omega);
    }
    expect_rvv10_added(total, semilocal, chr, point_tolerance);
}

void expect_distributed_matches_complete(ModulePW::PW_Basis& distributed, ModulePW::PW_Basis& complete)
{
    const std::vector<double> valence = local_grid(distributed, valence_density);
    const std::vector<double> core = local_grid(distributed, core_density);
    const std::vector<double> valence_complete = local_grid(complete, valence_density);
    const std::vector<double> core_complete = local_grid(complete, core_density);
    const VxcResult result = rvv10_result(distributed, {valence}, core);
    const VxcResult reference = rvv10_result(complete, {valence_complete}, core_complete);
    EXPECT_NEAR(std::get<0>(result), std::get<0>(reference), rounding_tolerance);
    EXPECT_NEAR(std::get<1>(result), std::get<1>(reference), rounding_tolerance);
    for (int ixy = 0; ixy < distributed.nxy; ++ixy)
    {
        for (int iz = 0; iz < distributed.nplane; ++iz)
        {
            const int local = ixy * distributed.nplane + iz;
            const int global = ixy * distributed.nz + distributed.startz_current + iz;
            EXPECT_NEAR(std::get<2>(result)(0, local), std::get<2>(reference)(0, global), rounding_tolerance);
        }
    }
    const ModuleBase::matrix stress = rvv10_stress(distributed, {valence}, core);
    const ModuleBase::matrix stress_complete = rvv10_stress(complete, {valence_complete}, core_complete);
    EXPECT_GT(std::abs(stress_complete(0, 0)), 1.0e-10);
    for (int a = 0; a < 3; ++a)
    {
        for (int b = 0; b < 3; ++b)
        {
            EXPECT_NEAR(stress(a, b), stress_complete(a, b), rounding_tolerance);
        }
    }
}

ChargeOnGrid::ChargeOnGrid(ModulePW::PW_Basis& pw,
                           const std::vector<std::vector<double>>& channels,
                           const std::vector<double>& core)
{
    channels_ = channels;
    core_ = core;
    const std::size_t ncomponents = channels.size();
    const std::vector<std::complex<double>> reciprocal(pw.npw);
    rhog_.assign(ncomponents, reciprocal);
    core_g_ = reciprocal;
    for (std::size_t is = 0; is < ncomponents; ++is)
    {
        rows_.push_back(channels_[is].data());
        rhog_rows_.push_back(rhog_[is].data());
    }
    chr.rhopw = &pw;
    chr.nrxx = pw.nrxx;
    chr.nspin = static_cast<int>(ncomponents);
    chr.rho = rows_.data();
    chr.rhog = rhog_rows_.data();
    chr.rho_core = core_.data();
    chr.rhog_core = core_g_.data();
}

RuntimeParameters::RuntimeParameters(const bool rvv10)
{
    saved_ = XC_Functional::get_runtime_parameters();
    const XCFunctionalParameters parameters = rvv10_parameters(rvv10);
    XC_Functional::set_runtime_parameters(parameters);
}

RuntimeParameters::~RuntimeParameters()
{
    XC_Functional::set_runtime_parameters(saved_);
}

void NonlocalXCTest::SetUp()
{
    // The oblique cell catches row/column mistakes in the cell handed to libvdwxc.
    const ModuleBase::Matrix3 lattice = oblique_lattice(1.0);
    pw.initmpi(test_size, test_rank, MPI_COMM_WORLD);
    serial.initmpi(1, 0, MPI_COMM_SELF);
    setup(pw, lattice, basis_cutoff);
    setup(serial, lattice, basis_cutoff);
}

} // namespace rvv10_test
