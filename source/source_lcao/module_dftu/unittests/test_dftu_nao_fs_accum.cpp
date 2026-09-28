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

#include "../dftu_nao_fs_accum.h"

#include "source_base/matrix.h"
#include "source_basis/module_ao/parallel_orbitals.h"
#include "gtest/gtest.h"

#include <complex>
#include <vector>

namespace DFTU_LCAO
{
namespace
{

class DftuFsAccumTest : public ::testing::Test
{
  protected:
    void SetUp() override
    {
        // 2 atoms, 3x3 serial local block with identity local<->global map
        pv.set_serial(3, 3);

        // iwt2iat: orbital 0,1 -> atom 0; orbital 2 -> atom 1
        // Must use new[] because Statistics destructor will delete[] it.
        ucell.iwt2iat = new int[3]{0, 0, 1};
    }

    Parallel_Orbitals pv;
    UnitCell ucell;
};

TEST_F(DftuFsAccumTest, DiagForceDouble)
{
    // dm is column-major local block: dm[ic * nrow + ir]
    // Diagonal entries (ir == ic): dm[0], dm[4], dm[8]
    std::vector<double> dm(9, 0.0);
    dm[0] = 1.0; // (0,0) -> atom 0
    dm[4] = 2.0; // (1,1) -> atom 0
    dm[8] = 3.0; // (2,2) -> atom 1
    // Off-diagonal entries should be ignored
    dm[1] = 10.0;
    dm[3] = 20.0;

    ModuleBase::matrix force(2, 3);
    accumulate_diag_force(pv, ucell, dm.data(), 0, force);

    EXPECT_DOUBLE_EQ(force(0, 0), 1.0 + 2.0);
    EXPECT_DOUBLE_EQ(force(1, 0), 3.0);
    EXPECT_DOUBLE_EQ(force(0, 1), 0.0);
    EXPECT_DOUBLE_EQ(force(1, 1), 0.0);
}

TEST_F(DftuFsAccumTest, DiagForceComplex)
{
    std::vector<std::complex<double>> dm(9, std::complex<double>(0.0, 0.0));
    dm[0] = std::complex<double>(1.0, 0.5);
    dm[4] = std::complex<double>(2.0, 1.0);
    dm[8] = std::complex<double>(3.0, -0.5);

    ModuleBase::matrix force(2, 3);
    accumulate_diag_force(pv, ucell, dm.data(), 1, force);

    // Only real part contributes
    EXPECT_DOUBLE_EQ(force(0, 1), 1.0 + 2.0);
    EXPECT_DOUBLE_EQ(force(1, 1), 3.0);
    EXPECT_DOUBLE_EQ(force(0, 0), 0.0);
}

TEST_F(DftuFsAccumTest, DiagStressDouble)
{
    std::vector<double> dm(9, 0.0);
    dm[0] = 1.0;
    dm[4] = 2.0;
    dm[8] = 3.0;

    ModuleBase::matrix stress(3, 3);
    accumulate_diag_stress(pv, dm.data(), 0, 1, 2.0, stress);

    // factor = 2.0 applied to all diagonal entries
    EXPECT_DOUBLE_EQ(stress(0, 1), 2.0 * (1.0 + 2.0 + 3.0));
    EXPECT_DOUBLE_EQ(stress(0, 0), 0.0);
    EXPECT_DOUBLE_EQ(stress(1, 1), 0.0);
}

TEST_F(DftuFsAccumTest, DiagStressComplex)
{
    std::vector<std::complex<double>> dm(9, std::complex<double>(0.0, 0.0));
    dm[0] = std::complex<double>(1.0, 0.5);
    dm[4] = std::complex<double>(2.0, 1.0);
    dm[8] = std::complex<double>(3.0, -0.5);

    ModuleBase::matrix stress(3, 3);
    accumulate_diag_stress(pv, dm.data(), 2, 2, -0.5, stress);

    // factor = -0.5, real part only
    EXPECT_DOUBLE_EQ(stress(2, 2), -0.5 * (1.0 + 2.0 + 3.0));
    EXPECT_DOUBLE_EQ(stress(0, 0), 0.0);
}

} // namespace
} // namespace DFTU_LCAO
