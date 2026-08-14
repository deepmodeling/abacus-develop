#include "gtest/gtest.h"

#include "source_cell/atom_pseudo.h"
#include "source_cell/atom_spec.h"
#include "source_cell/magnetism.h"
#include "source_cell/pseudo.h"
#include "source_cell/unitcell.h"
#include "source_cell/module_symmetry/little_group.h"

#include "source_base/global_variable.h"
#include <cmath>
#include <vector>

// Linker stubs: the symmetry library needs these symbols (see klist_test.cpp).
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
SepPot::SepPot() {}
SepPot::~SepPot() {}
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
Sep_Cell::Sep_Cell() noexcept {}
Sep_Cell::~Sep_Cell() noexcept {}

// abbreviated from module_symmetry/test/symm_test.cpp
struct atomtype_
{
    std::string atomname;
    std::vector<std::vector<double>> coordinate;
};

struct stru_
{
    int ibrav;
    std::string point_group;    // Schoenflies symbol
    std::string point_group_hm; // Hermann-Mauguin notation.
    std::string space_group;
    std::vector<double> cell;
    std::vector<atomtype_> all_type;
};

// primitive cubic with one atom at the origin -> O_h (48 operations)
std::vector<stru_> stru_lib{stru_{1,
                                  "O_h",
                                  "m-3m",
                                  "Pm-3m",
                                  std::vector<double>{1., 0., 0., 0., 1., 0., 0., 0., 1.},
                                  std::vector<atomtype_>{atomtype_{"C",
                                                                   std::vector<std::vector<double>>{
                                                                       {0., 0., 0.},
                                                                   }}}}};

class LittleGroupTest : public testing::Test
{
  protected:
    UnitCell ucell;
    ModuleSymmetry::Symmetry symm;

    void SetUp() override
    {
        std::vector<atomtype_> coord = stru_lib[0].all_type;
        ucell.a1 = ModuleBase::Vector3<double>(stru_lib[0].cell[0], stru_lib[0].cell[1], stru_lib[0].cell[2]);
        ucell.a2 = ModuleBase::Vector3<double>(stru_lib[0].cell[3], stru_lib[0].cell[4], stru_lib[0].cell[5]);
        ucell.a3 = ModuleBase::Vector3<double>(stru_lib[0].cell[6], stru_lib[0].cell[7], stru_lib[0].cell[8]);
        ucell.ntype = stru_lib[0].all_type.size();
        ucell.atoms = new Atom[ucell.ntype];
        ucell.nat = 0;
        ucell.latvec.e11 = ucell.a1.x;
        ucell.latvec.e12 = ucell.a1.y;
        ucell.latvec.e13 = ucell.a1.z;
        ucell.latvec.e21 = ucell.a2.x;
        ucell.latvec.e22 = ucell.a2.y;
        ucell.latvec.e23 = ucell.a2.z;
        ucell.latvec.e31 = ucell.a3.x;
        ucell.latvec.e32 = ucell.a3.y;
        ucell.latvec.e33 = ucell.a3.z;
        ucell.GT = ucell.latvec.Inverse();
        ucell.G = ucell.GT.Transpose();
        ucell.lat0 = 1.8897261254578281;
        for (int i = 0; i < coord.size(); i++)
        {
            ucell.atoms[i].label = coord[i].atomname;
            ucell.atoms[i].na = coord[i].coordinate.size();
            ucell.atoms[i].tau.resize(ucell.atoms[i].na);
            ucell.atoms[i].taud.resize(ucell.atoms[i].na);
            for (int j = 0; j < ucell.atoms[i].na; j++)
            {
                std::vector<double> this_atom = coord[i].coordinate[j];
                ucell.atoms[i].tau[j] = ModuleBase::Vector3<double>(this_atom[0], this_atom[1], this_atom[2]);
                ucell.atoms[i].taud[j] = ModuleBase::Vector3<double>(0.0, 0.0, 0.0);
            }
            ucell.nat += ucell.atoms[i].na;
        }
        GlobalV::ofs_running.open("tmp_little_group");
        const int cal_symm_repr[2] = {0, 6};
        symm.analy_sys(ucell.lat, ucell.st, ucell.atoms, GlobalV::ofs_running, 1e-6, 1, "scf", cal_symm_repr);
        GlobalV::ofs_running.close();
    }

    void TearDown() override
    {
        delete[] ucell.atoms;
        remove("tmp_little_group");
    }
};

TEST_F(LittleGroupTest, LittleGroupSizeAtKnownPoints)
{
    // the full space group of the primitive cubic cell
    EXPECT_EQ(symm.nrotk, 48);

    // Gamma: all operations keep (0,0,0)
    ModuleSymmetry::LittleGroup lg;
    lg.set_q(ModuleBase::Vector3<double>(0.0, 0.0, 0.0), symm);
    EXPECT_EQ(lg.get_little_group_ops().size(), 48);

    // R point (1/2,1/2,1/2): -q is congruent to q (mod 1), so also all 48
    lg.set_q(ModuleBase::Vector3<double>(0.5, 0.5, 0.5), symm);
    EXPECT_EQ(lg.get_little_group_ops().size(), 48);

    // X point (1/2,0,0): D_4h little group, 16 operations
    lg.set_q(ModuleBase::Vector3<double>(0.5, 0.0, 0.0), symm);
    EXPECT_EQ(lg.get_little_group_ops().size(), 16);

    // M point (1/2,1/2,0): D_4h little group, 16 operations
    lg.set_q(ModuleBase::Vector3<double>(0.5, 0.5, 0.0), symm);
    EXPECT_EQ(lg.get_little_group_ops().size(), 16);

    // generic point: only the identity
    lg.set_q(ModuleBase::Vector3<double>(0.13, 0.27, 0.41), symm);
    EXPECT_EQ(lg.get_little_group_ops().size(), 1);
}

TEST_F(LittleGroupTest, PlaceholderIrrep)
{
    ModuleSymmetry::LittleGroup lg;
    lg.set_q(ModuleBase::Vector3<double>(0.0, 0.0, 0.0), symm);

    EXPECT_EQ(lg.get_nirr(), 1); // fully-symmetric A1 placeholder
    EXPECT_TRUE(lg.get_mode_basis(0).empty());
    EXPECT_EQ(lg.get_q(), ModuleBase::Vector3<double>(0.0, 0.0, 0.0));
}
