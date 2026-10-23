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

#include "../dftu_nao_folding.h"

#include "source_basis/module_ao/parallel_orbitals.h"
#include "gtest/gtest.h"

#include <string>

// Minimal stub for get_linear_index to keep the link closure small; only the
// index arithmetic is under test. Parallel_Orbitals' ctor/dtor come from the
// real parallel_orbitals.cpp (wired in via CMakeLists) so the class layout
// matches the MPI-built base library.
namespace DFTU_LCAO
{
int get_linear_index(const std::string& ks_solver,
                     const int mu,
                     const int nu,
                     const Parallel_Orbitals& pv)
{
    if (ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER(ks_solver))
    {
        return mu + nu * pv.nrow;
    }
    return mu * pv.ncol + nu;
}
} // namespace DFTU_LCAO

namespace DFTU_LCAO
{

TEST(GetLinearIndexTest, RowMajor)
{
    Parallel_Orbitals pv;
    pv.nrow = 3;
    pv.ncol = 4;

    // row-major: mu * ncol + nu ("cg" is not column-major)
    EXPECT_EQ(get_linear_index("cg", 0, 0, pv), 0);
    EXPECT_EQ(get_linear_index("cg", 0, 1, pv), 1);
    EXPECT_EQ(get_linear_index("cg", 1, 0, pv), 4);
    EXPECT_EQ(get_linear_index("cg", 1, 2, pv), 6);
    EXPECT_EQ(get_linear_index("cg", 2, 3, pv), 11);
}

TEST(GetLinearIndexTest, ColumnMajor)
{
    Parallel_Orbitals pv;
    pv.nrow = 3;
    pv.ncol = 4;

    // column-major: mu + nu * nrow
    EXPECT_EQ(get_linear_index("scalapack_gvx", 0, 0, pv), 0);
    EXPECT_EQ(get_linear_index("scalapack_gvx", 1, 0, pv), 1);
    EXPECT_EQ(get_linear_index("scalapack_gvx", 0, 1, pv), 3);
    EXPECT_EQ(get_linear_index("scalapack_gvx", 2, 1, pv), 5);
    EXPECT_EQ(get_linear_index("scalapack_gvx", 2, 3, pv), 11);
}

} // namespace DFTU_LCAO
