#include "source_cell/unitcell.h"

// Minimal mocks to satisfy the linker without pulling in heavy dependencies.
Magnetism::Magnetism() {}
Magnetism::~Magnetism() {}
SepPot::SepPot() {}
SepPot::~SepPot() {}
Sep_Cell::Sep_Cell() noexcept {}
Sep_Cell::~Sep_Cell() noexcept {}
UnitCell::UnitCell() {}
UnitCell::~UnitCell() {}

#include "../dftu_nao_fs_k.h"

#include "source_lcao/force_stress_arrays.h"
#include "source_pw/module_pwdft/dftu_base.h"

#include "gtest/gtest.h"

#include <string>
#include <vector>

/***********************************************************************
 * Unit tests for dftu_nao_fs_k.cpp.
 *
 * The bulk of the file (force_stress, cal_force_k, cal_stress_k) lives in
 * an anonymous namespace or needs Grid_Driver / BLACS / TwoCenterIntegrator
 * and is explicitly documented as hard to unit-test in the source. The
 * only publicly visible surface is the DftuFsEnv environment struct, whose
 * accessors bundle references to the shared dependencies. We verify that
 * the constructor stores each reference unchanged so downstream kernels
 * see the same objects the caller passed in.
 *
 * Note: we allocate the dependency objects on the heap and never destroy
 * them, so the linker never sees references to their (heavy) constructors
 * or destructors. This keeps the test target free of BLACS / neighbor /
 * two-center linkage while still exercising the reference semantics of
 * DftuFsEnv.
 ***********************************************************************/

namespace DFTU_LCAO
{
namespace
{

class DftuFsKTest : public ::testing::Test
{
  protected:
    void SetUp() override
    {
        // Intentionally leaked: their destructors pull in heavy symbols
        // (BLACS teardown, neighbor-list cleanup, ...). The test process
        // is short-lived, so the leak is harmless.
        dftu = static_cast<Plus_U_Base*>(::operator new(sizeof(Plus_U_Base)));
        fsr = static_cast<ForceStressArrays*>(::operator new(sizeof(ForceStressArrays)));
        gd = static_cast<Grid_Driver*>(::operator new(sizeof(Grid_Driver)));
        pv = static_cast<Parallel_Orbitals*>(::operator new(sizeof(Parallel_Orbitals)));
    }

    Plus_U_Base* dftu = nullptr;
    ForceStressArrays* fsr = nullptr;
    Grid_Driver* gd = nullptr;
    Parallel_Orbitals* pv = nullptr;
    UnitCell ucell;
};

TEST_F(DftuFsKTest, EnvStoresReferencesUnchanged)
{
    std::vector<double> orb_cutoff = {5.0, 6.0};
    const std::string solver = "scalapack";

    DftuFsEnv env(*dftu, ucell, *gd, *pv, *fsr, orb_cutoff, solver);

    EXPECT_EQ(&env.dftu(), dftu);
    EXPECT_EQ(&env.ucell(), &ucell);
    EXPECT_EQ(&env.gd(), gd);
    EXPECT_EQ(&env.pv(), pv);
    EXPECT_EQ(&env.fsr(), fsr);
    ASSERT_EQ(env.orb_cutoff().size(), 2u);
    EXPECT_DOUBLE_EQ(env.orb_cutoff()[0], 5.0);
    EXPECT_DOUBLE_EQ(env.orb_cutoff()[1], 6.0);
    EXPECT_EQ(env.ks_solver(), "scalapack");
}

TEST_F(DftuFsKTest, EnvReflectsExternalMutation)
{
    // The env holds references, not copies: mutating the caller-side
    // objects must be visible through the env accessors.
    std::vector<double> orb_cutoff = {1.0};
    DftuFsEnv env(*dftu, ucell, *gd, *pv, *fsr, orb_cutoff, "genelpa");

    orb_cutoff[0] = 2.5;
    EXPECT_DOUBLE_EQ(env.orb_cutoff()[0], 2.5);

    ucell.lat0 = 7.0;
    EXPECT_DOUBLE_EQ(env.ucell().lat0, 7.0);
}

} // namespace
} // namespace DFTU_LCAO
