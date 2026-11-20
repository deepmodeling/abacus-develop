// vdw.cpp: the rVV10 entries of vdw::is_rvv10 and vdw::make_vdw. rVV10 is evaluated inside the XC functional
// (nonlocal_xc.cpp), so make_vdw returns no pairwise correction; without libvdwxc it stops with a rebuild hint.
#include "source_cell/unitcell.h"
#include "source_hamilt/module_vdw/vdw.h"
#include "source_io/module_parameter/input_parameter.h"

#include "gtest/gtest.h"
#include <string>

// make_vdw does not read the cell for rVV10; the target links the vdW objects only.
pseudo::pseudo()
{
}
pseudo::~pseudo()
{
}
Atom::Atom()
{
}
Atom::~Atom()
{
}
Atom_pseudo::Atom_pseudo()
{
}
Atom_pseudo::~Atom_pseudo()
{
}
UnitCell::UnitCell()
{
}
UnitCell::~UnitCell()
{
}
Magnetism::Magnetism()
{
}
Magnetism::~Magnetism()
{
}
SepPot::SepPot()
{
}
SepPot::~SepPot()
{
}
Sep_Cell::Sep_Cell() noexcept
{
}
Sep_Cell::~Sep_Cell() noexcept
{
}

TEST(VdwRvv10, IsRvv10SelectsOnlyTheExactName)
{
    EXPECT_TRUE(vdw::is_rvv10("rvv10"));
    EXPECT_FALSE(vdw::is_rvv10("RVV10"));
    EXPECT_FALSE(vdw::is_rvv10("d3_bj"));
    EXPECT_FALSE(vdw::is_rvv10("none"));
}

TEST(VdwRvv10, MakeVdwHasNoPairwiseTerm)
{
    UnitCell ucell;
    Input_para input;
    input.vdw_method = "rvv10";
#ifdef __LIBVDWXC
    EXPECT_TRUE(vdw::make_vdw(ucell, input) == nullptr);
#else
    testing::internal::CaptureStdout();
    EXPECT_EXIT(vdw::make_vdw(ucell, input), ::testing::ExitedWithCode(1), "");
    const std::string output = testing::internal::GetCapturedStdout();
#endif
}
