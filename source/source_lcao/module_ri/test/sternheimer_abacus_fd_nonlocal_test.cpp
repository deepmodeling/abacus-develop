#include "source_lcao/module_ri/sternheimer_abacus_fd_nonlocal.h"

#include "source_base/matrix.h"
#include "source_cell/atom_spec.h"
#include "source_cell/unitcell.h"

#include <stdexcept>
#include <gtest/gtest.h>

Magnetism::Magnetism() {}
Magnetism::~Magnetism() {}

#ifdef __LCAO
InfoNonlocal::InfoNonlocal() {}
InfoNonlocal::~InfoNonlocal() {}
#endif

namespace
{

constexpr double y00()
{
    return 0.28209479177387814347;
}

ModuleRI::SternheimerFDHamiltonian::Grid line_grid()
{
    ModuleRI::SternheimerFDHamiltonian::Grid grid;
    grid.nx = 3;
    grid.ny = 1;
    grid.nz = 1;
    grid.hx = 1.0;
    grid.hy = 1.0;
    grid.hz = 1.0;
    grid.periodic = false;
    return grid;
}

void fill_single_s_projector(Atom& atom)
{
    atom.ncpp.has_so = false;
    atom.ncpp.tvanp = false;
    atom.ncpp.nbeta = 1;
    atom.ncpp.mesh = 3;
    atom.ncpp.r = {0.0, 1.0, 2.0};
    atom.ncpp.lll = {0};
    atom.ncpp.betar.create(1, 3);
    atom.ncpp.betar(0, 0) = 2.0;
    atom.ncpp.betar(0, 1) = 4.0;
    atom.ncpp.betar(0, 2) = 6.0;
    atom.ncpp.dion.create(1, 1);
    atom.ncpp.dion(0, 0) = 1.5;
}

void fill_one_atom_cell(UnitCell& ucell)
{
    ucell.ntype = 1;
    ucell.nat = 1;
    ucell.set_atom_flag = true;
    ucell.atoms = new Atom[1];
    ucell.atoms[0].na = 1;
    ucell.atoms[0].tau = {ModuleBase::Vector3<double>(0.0, 0.0, 0.0)};
    fill_single_s_projector(ucell.atoms[0]);
}

} // namespace

TEST(SternheimerABACUSFDNonlocal, BuildsProjectorFromUnitCellAtomPseudo)
{
    UnitCell ucell;
    fill_one_atom_cell(ucell);

    const auto projector
        = ModuleRI::make_sternheimer_fd_nonlocal_projector_from_unitcell(ucell, line_grid(), 1.0);

    ASSERT_NE(projector, nullptr);
    ASSERT_EQ(projector->blocks().size(), 1);
    ASSERT_EQ(projector->blocks()[0].projectors.size(), 1);
    EXPECT_NEAR(projector->blocks()[0].projectors[0][0].real(), 2.0 * y00(), 1.0e-14);
    EXPECT_NEAR(projector->blocks()[0].projectors[0][1].real(), 4.0 * y00(), 1.0e-14);
    EXPECT_NEAR(projector->blocks()[0].projectors[0][2].real(), 6.0 * y00(), 1.0e-14);
    EXPECT_NEAR(projector->blocks()[0].d_matrix[0][0].real(), 1.5, 1.0e-14);
}

TEST(SternheimerABACUSFDNonlocal, ReturnsNullWhenNoAtomsCarryProjectors)
{
    UnitCell ucell;
    ucell.ntype = 1;
    ucell.nat = 1;
    ucell.set_atom_flag = true;
    ucell.atoms = new Atom[1];
    ucell.atoms[0].na = 1;
    ucell.atoms[0].tau = {ModuleBase::Vector3<double>(0.0, 0.0, 0.0)};

    const auto projector
        = ModuleRI::make_sternheimer_fd_nonlocal_projector_from_unitcell(ucell, line_grid(), 1.0);

    EXPECT_EQ(projector, nullptr);
}

TEST(SternheimerABACUSFDNonlocal, RejectsUnsupportedUltrasoftAndSOCPseudopotentials)
{
    UnitCell ucell;
    fill_one_atom_cell(ucell);
    ucell.atoms[0].ncpp.tvanp = true;
    EXPECT_THROW(ModuleRI::make_sternheimer_fd_nonlocal_projector_from_unitcell(ucell, line_grid(), 1.0),
                 std::invalid_argument);

    ucell.atoms[0].ncpp.tvanp = false;
    ucell.atoms[0].ncpp.has_so = true;
    EXPECT_THROW(ModuleRI::make_sternheimer_fd_nonlocal_projector_from_unitcell(ucell, line_grid(), 1.0),
                 std::invalid_argument);
}
