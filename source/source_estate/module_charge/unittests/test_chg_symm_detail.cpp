#include "gtest/gtest.h"

#include "source_base/matrix3.h"
#include "source_base/vector3.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_cell/magnetism.h"
#include "source_cell/module_symmetry/symmetry.h"
#include "source_cell/module_symmetry/symm_rot_spin.h"
#include "source_estate/module_charge/chg_symm_detail.h"

#include <array>
#include <complex>
#include <vector>

// unitcell.cpp (pulled in via the cell_info object library) references
// Magnetism; provide a lightweight stub, mirroring test_chg_symm.cpp.
Magnetism::Magnetism()
{
    this->tot_mag = 0.0;
    this->abs_mag = 0.0;
}
Magnetism::~Magnetism()
{
}

/************************************************
 *  unit test of module_charge/chg_symm_detail.cpp
 ***********************************************/

/**
 * - Tested Functions:
 *   - psymmg: symmetrize one reciprocal-space density component
 *   - psymmg_soc: symmetrize three coupled spin components for nspin=4
 *
 * Both are driven with a manually built D_4 point group on a cubic lattice
 * (a=1), mirroring source_cell/module_symmetry/test/symm_rho_soc_test.cpp.
 * The PW_Basis is serial (single plane-wave per FFT point) so the
 * non-MPI path in psymmg/psymmg_soc is exercised.
 *
 * Checks:
 *   - Idempotence: symmetrizing twice gives the same result as once.
 */

namespace
{

// the 8 proper rotations of D_4 (column-vector convention r' = R r)
const std::array<std::array<std::array<int, 3>, 3>, 8> kRcol = {{
    {{{ 1, 0, 0}, { 0, 1, 0}, { 0, 0, 1}}},   // E
    {{{ 0,-1, 0}, { 1, 0, 0}, { 0, 0, 1}}},   // C4z
    {{{-1, 0, 0}, { 0,-1, 0}, { 0, 0, 1}}},   // C2z
    {{{ 0, 1, 0}, {-1, 0, 0}, { 0, 0, 1}}},   // C4z^3
    {{{ 1, 0, 0}, { 0,-1, 0}, { 0, 0,-1}}},   // C2x
    {{{-1, 0, 0}, { 0, 1, 0}, { 0, 0,-1}}},   // C2y
    {{{ 0, 1, 0}, { 1, 0, 0}, { 0, 0,-1}}},   // C2[110]
    {{{ 0,-1, 0}, {-1, 0, 0}, { 0, 0,-1}}},   // C2[1-10]
}};

ModuleBase::Matrix3 gmatc_of(int g)
{
    const auto& r = kRcol[g];
    return ModuleBase::Matrix3(r[0][0], r[1][0], r[2][0],
                               r[0][1], r[1][1], r[2][1],
                               r[0][2], r[1][2], r[2][2]);
}

/// Build a D_4 symmetry group on a cubic lattice (a=1).
void build_d4_group(ModuleSymmetry::Symmetry& symm)
{
    symm.epsilon = 1e-6;
    symm.nrot = 8;
    symm.nrotk = 8;
    symm.nrotk_anti = 0;
    symm.ncell = 1;
    symm.ptrans = {ModuleBase::Vector3<double>(0.0, 0.0, 0.0)};
    ModuleSymmetry::Symmetry::pricell_loop = false;
    for (int g = 0; g < 8; ++g)
    {
        const ModuleBase::Matrix3 gc = gmatc_of(g);
        symm.gmatrix[g] = gc;
        symm.kgmatrix[g] = gc;
        symm.gtrans[g] = ModuleBase::Vector3<double>(0.0, 0.0, 0.0);
    }
}

} // namespace

class ChgSymmDetailTest : public ::testing::Test
{
  protected:
    ModulePW::PW_Basis pw_basis;
    ModuleSymmetry::Symmetry symm;

    void SetUp() override
    {
        pw_basis.initgrids(1.0, ModuleBase::Matrix3(1, 0, 0, 0, 1, 0, 0, 0, 1), 10.0);
        pw_basis.initparameters(false, 10.0);
        pw_basis.setuptransform();
        pw_basis.collect_local_pw();
        build_d4_group(symm);
    }
};

TEST_F(ChgSymmDetailTest, PsymmgIdempotence)
{
    std::vector<std::complex<double>> rhog(pw_basis.npw);
    for (int ig = 0; ig < pw_basis.npw; ++ig)
    {
        rhog[ig] = std::complex<double>(0.3 * ig - 1.0, 0.7 * ((ig * 13) % 5) - 1.5);
    }
    std::vector<std::complex<double>> once = rhog;
    std::vector<std::complex<double>> twice = rhog;

    module_charge::detail::psymmg(once.data(), &pw_basis, symm);
    module_charge::detail::psymmg(twice.data(), &pw_basis, symm);
    module_charge::detail::psymmg(twice.data(), &pw_basis, symm);

    for (int ig = 0; ig < pw_basis.npw; ++ig)
    {
        EXPECT_NEAR(once[ig].real(), twice[ig].real(), 1e-8);
        EXPECT_NEAR(once[ig].imag(), twice[ig].imag(), 1e-8);
    }
}

TEST_F(ChgSymmDetailTest, PsymmgSocIdempotence)
{
    std::vector<std::complex<double>> x(pw_basis.npw);
    std::vector<std::complex<double>> y(pw_basis.npw);
    std::vector<std::complex<double>> z(pw_basis.npw);
    for (int ig = 0; ig < pw_basis.npw; ++ig)
    {
        x[ig] = std::complex<double>(0.3 * ig - 1.0, 0.7 * ((ig * 13) % 5) - 1.5);
        y[ig] = std::complex<double>(-0.5 * ((ig * 7) % 4) + 0.9, 0.2 * ig - 2.0);
        z[ig] = std::complex<double>(0.11 * ((ig * 3) % 6), -0.4 * ((ig * 5) % 7) + 1.0);
    }
    std::vector<std::complex<double>> x1 = x, y1 = y, z1 = z;

    module_charge::detail::psymmg_soc(x.data(), y.data(), z.data(), &pw_basis, symm);
    module_charge::detail::psymmg_soc(x1.data(), y1.data(), z1.data(), &pw_basis, symm);
    module_charge::detail::psymmg_soc(x1.data(), y1.data(), z1.data(), &pw_basis, symm);

    for (int ig = 0; ig < pw_basis.npw; ++ig)
    {
        EXPECT_NEAR(x[ig].real(), x1[ig].real(), 1e-8);
        EXPECT_NEAR(x[ig].imag(), x1[ig].imag(), 1e-8);
        EXPECT_NEAR(y[ig].real(), y1[ig].real(), 1e-8);
        EXPECT_NEAR(y[ig].imag(), y1[ig].imag(), 1e-8);
        EXPECT_NEAR(z[ig].real(), z1[ig].real(), 1e-8);
        EXPECT_NEAR(z[ig].imag(), z1[ig].imag(), 1e-8);
    }
}
