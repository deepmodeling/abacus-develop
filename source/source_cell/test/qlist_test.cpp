#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include <iostream>
#include <streambuf>
#define private public
#include "source_cell/atom_pseudo.h"
#include "source_cell/atom_spec.h"
#include "source_cell/pseudo.h"
#include "source_cell/qlist.h"
#include "source_cell/unitcell.h"
#include "source_cell/magnetism.h"
#undef private
#include "source_base/mathzone.h"
#include "source_base/parallel_global.h"
#include "source_base/global_variable.h"

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

/************************************************
 *  unit test of class QList
 ***********************************************/

/**
 * - Tested Functions:
 *   - generate_mesh()
 *     - the Monkhorst-Pack q-point mesh is generated and reduced by star
 *       (time-reversal included)
 *   - get_nq() / get_q()
 *     - access the reduced q-point list
 *   - get_nirr() / get_irrep_modes()
 *     - placeholder irrep data (one fully-symmetric irrep per q-point)
 *   - read_from_file()
 *     - placeholder interface, must not crash
 */

// abbreviated from module_symmetry/test/symm_test.cpp and klist_test.cpp
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

std::vector<stru_> stru_lib{stru_{1,
                                  "O_h",
                                  "m-3m",
                                  "Pm-3m",
                                  std::vector<double>{1., 0., 0., 0., 1., 0., 0., 0., 1.},
                                  std::vector<atomtype_>{atomtype_{"C",
                                                                   std::vector<std::vector<double>>{
                                                                       {0., 0., 0.},
                                                                   }}}}};

class QListTest : public testing::Test
{
  protected:
    ModuleCell::QList qlist;
    std::ifstream ifs;
    std::ofstream ofs;
    std::ofstream ofs_running;
    std::string output;

    UnitCell ucell;
    void construct_ucell(stru_& stru)
    {
        std::vector<atomtype_> coord = stru.all_type;
        ucell.a1 = ModuleBase::Vector3<double>(stru.cell[0], stru.cell[1], stru.cell[2]);
        ucell.a2 = ModuleBase::Vector3<double>(stru.cell[3], stru.cell[4], stru.cell[5]);
        ucell.a3 = ModuleBase::Vector3<double>(stru.cell[6], stru.cell[7], stru.cell[8]);
        ucell.ntype = stru.all_type.size();
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
                ModuleBase::Mathzone::Cartesian_to_Direct(ucell.atoms[i].tau[j].x,
                                                          ucell.atoms[i].tau[j].y,
                                                          ucell.atoms[i].tau[j].z,
                                                          ucell.a1.x,
                                                          ucell.a1.y,
                                                          ucell.a1.z,
                                                          ucell.a2.x,
                                                          ucell.a2.y,
                                                          ucell.a2.z,
                                                          ucell.a3.x,
                                                          ucell.a3.y,
                                                          ucell.a3.z,
                                                          ucell.atoms[i].taud[j].x,
                                                          ucell.atoms[i].taud[j].y,
                                                          ucell.atoms[i].taud[j].z);
            }
            ucell.nat += ucell.atoms[i].na;
        }
    }

    void ClearUcell()
    {
        delete[] ucell.atoms;
    }
};

TEST_F(QListTest, GenerateMeshFullSymmetry)
{
    construct_ucell(stru_lib[0]);
    GlobalV::ofs_running.open("tmp_qlist_1");
    ModuleSymmetry::Symmetry symm;
    const int cal_symm_repr[2] = {0, 6};
    symm.analy_sys(ucell.lat, ucell.st, ucell.atoms, GlobalV::ofs_running, 1e-6, 1, "scf", cal_symm_repr);

    qlist.generate_mesh(ucell, symm, {8, 8, 8}, true);

    // full mesh 512 -> irreducible q-points of the primitive cubic lattice
    EXPECT_EQ(qlist.nkstot_full, 512);
    EXPECT_EQ(qlist.get_nq(), 35);
    EXPECT_EQ(qlist.get_nq(), qlist.nkstot);
    EXPECT_TRUE(qlist.is_mp);

    // weights must sum to 1 after normalization inside generate_mesh
    double sum = 0.0;
    for (int i = 0; i < qlist.get_nq(); ++i)
    {
        sum += qlist.wk[i];
    }
    EXPECT_NEAR(sum, 1.0, 1e-10);

    // q-points must be unique
    for (int i = 0; i < qlist.get_nq(); ++i)
    {
        for (int j = i + 1; j < qlist.get_nq(); ++j)
        {
            EXPECT_FALSE(qlist.get_q(i) == qlist.get_q(j));
        }
    }

    GlobalV::ofs_running.close();
    ClearUcell();
    remove("tmp_qlist_1");
}

TEST_F(QListTest, GenerateMeshSmallGrid)
{
    construct_ucell(stru_lib[0]);
    GlobalV::ofs_running.open("tmp_qlist_2");
    ModuleSymmetry::Symmetry symm;
    const int cal_symm_repr[2] = {0, 6};
    symm.analy_sys(ucell.lat, ucell.st, ucell.atoms, GlobalV::ofs_running, 1e-6, 1, "scf", cal_symm_repr);

    qlist.generate_mesh(ucell, symm, {2, 2, 2}, true);

    // {0,0.5}^3 under O_h folds to Gamma + X + M + R
    EXPECT_EQ(qlist.nkstot_full, 8);
    EXPECT_EQ(qlist.get_nq(), 4);

    // the first irreducible q-point must be Gamma (0,0,0)
    EXPECT_DOUBLE_EQ(qlist.get_q(0).x, 0.0);
    EXPECT_DOUBLE_EQ(qlist.get_q(0).y, 0.0);
    EXPECT_DOUBLE_EQ(qlist.get_q(0).z, 0.0);

    GlobalV::ofs_running.close();
    ClearUcell();
    remove("tmp_qlist_2");
}

TEST_F(QListTest, GammaOnlyGrid)
{
    construct_ucell(stru_lib[0]);
    GlobalV::ofs_running.open("tmp_qlist_3");
    ModuleSymmetry::Symmetry symm;
    const int cal_symm_repr[2] = {0, 6};
    symm.analy_sys(ucell.lat, ucell.st, ucell.atoms, GlobalV::ofs_running, 1e-6, 1, "scf", cal_symm_repr);

    qlist.generate_mesh(ucell, symm, {1, 1, 1}, true);

    EXPECT_EQ(qlist.nkstot_full, 1);
    EXPECT_EQ(qlist.get_nq(), 1);
    EXPECT_DOUBLE_EQ(qlist.wk[0], 1.0);
    EXPECT_DOUBLE_EQ(qlist.get_q(0).x, 0.0);

    GlobalV::ofs_running.close();
    ClearUcell();
    remove("tmp_qlist_3");
}

TEST_F(QListTest, IrrepPlaceholder)
{
    construct_ucell(stru_lib[0]);
    GlobalV::ofs_running.open("tmp_qlist_4");
    ModuleSymmetry::Symmetry symm;
    const int cal_symm_repr[2] = {0, 6};
    symm.analy_sys(ucell.lat, ucell.st, ucell.atoms, GlobalV::ofs_running, 1e-6, 1, "scf", cal_symm_repr);

    qlist.generate_mesh(ucell, symm, {2, 2, 2}, true);

    // placeholder: one fully-symmetric irrep per q-point, empty mode list
    for (int i = 0; i < qlist.get_nq(); ++i)
    {
        EXPECT_EQ(qlist.get_nirr(i), 1);
        EXPECT_TRUE(qlist.get_irrep_modes(i, 0).empty());
    }

    // out-of-range access must return an empty list instead of crashing
    EXPECT_TRUE(qlist.get_irrep_modes(-1, 0).empty());
    EXPECT_TRUE(qlist.get_irrep_modes(qlist.get_nq(), 0).empty());
    EXPECT_TRUE(qlist.get_irrep_modes(0, 5).empty());

    GlobalV::ofs_running.close();
    ClearUcell();
    remove("tmp_qlist_4");
}

TEST_F(QListTest, ReadFromFilePlaceholder)
{
    qlist.read_from_file("nonexistent_qpoints", ucell);
    EXPECT_EQ(qlist.get_nq(), 0);
}
