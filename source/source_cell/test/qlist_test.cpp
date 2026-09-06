#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include <fstream>
#include <iostream>
#include <iterator>
#include <streambuf>
#include <string>
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
 *     - ReadFromFileDirect: explicit Direct list with weights
 *     - ReadFromFileCartesian: explicit Cartesian list with weights
 *     - ReadFromFileMonkhorstPack: auto mesh (nkstot == 0)
 *     - ReadFromFileLinePath / ReadFromFileLineCartesian: line interpolation
 *     - ReadFromFileNegativeNkstot: negative count is rejected cleanly
 *     - ReadFromFileLineRejectsZeroCount: non-positive line count quits
 *     - ReadFromFileMissing / ReadFromFileBadHeader: must not crash
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
    ofs_running.open("tmp_qlist_1");
    ModuleSymmetry::Symmetry symm;
    const int cal_symm_repr[2] = {0, 6};
    symm.analy_sys(ucell.lat, ucell.st, ucell.atoms, ofs_running, 1e-6, 1, "scf", cal_symm_repr);

    qlist.generate_mesh(ucell, symm, {8, 8, 8}, true);

    // full mesh 512 -> irreducible q-points of the primitive cubic lattice
    EXPECT_EQ(qlist.nkstot_nospin, 512);
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

    ofs_running.close();
    ClearUcell();
    remove("tmp_qlist_1");
}

TEST_F(QListTest, GenerateMeshSmallGrid)
{
    construct_ucell(stru_lib[0]);
    ofs_running.open("tmp_qlist_2");
    ModuleSymmetry::Symmetry symm;
    const int cal_symm_repr[2] = {0, 6};
    symm.analy_sys(ucell.lat, ucell.st, ucell.atoms, ofs_running, 1e-6, 1, "scf", cal_symm_repr);

    qlist.generate_mesh(ucell, symm, {2, 2, 2}, true);

    // {0,0.5}^3 under O_h folds to Gamma + X + M + R
    EXPECT_EQ(qlist.nkstot_nospin, 8);
    EXPECT_EQ(qlist.get_nq(), 4);

    // the first irreducible q-point must be Gamma (0,0,0)
    EXPECT_DOUBLE_EQ(qlist.get_q(0).x, 0.0);
    EXPECT_DOUBLE_EQ(qlist.get_q(0).y, 0.0);
    EXPECT_DOUBLE_EQ(qlist.get_q(0).z, 0.0);

    ofs_running.close();
    ClearUcell();
    remove("tmp_qlist_2");
}

TEST_F(QListTest, GammaOnlyGrid)
{
    construct_ucell(stru_lib[0]);
    ofs_running.open("tmp_qlist_3");
    ModuleSymmetry::Symmetry symm;
    const int cal_symm_repr[2] = {0, 6};
    symm.analy_sys(ucell.lat, ucell.st, ucell.atoms, ofs_running, 1e-6, 1, "scf", cal_symm_repr);

    qlist.generate_mesh(ucell, symm, {1, 1, 1}, true);

    EXPECT_EQ(qlist.nkstot_nospin, 1);
    EXPECT_EQ(qlist.get_nq(), 1);
    EXPECT_DOUBLE_EQ(qlist.wk[0], 1.0);
    EXPECT_DOUBLE_EQ(qlist.get_q(0).x, 0.0);

    ofs_running.close();
    ClearUcell();
    remove("tmp_qlist_3");
}

TEST_F(QListTest, IrrepPlaceholder)
{
    construct_ucell(stru_lib[0]);
    ofs_running.open("tmp_qlist_4");
    ModuleSymmetry::Symmetry symm;
    const int cal_symm_repr[2] = {0, 6};
    symm.analy_sys(ucell.lat, ucell.st, ucell.atoms, ofs_running, 1e-6, 1, "scf", cal_symm_repr);

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

    ofs_running.close();
    ClearUcell();
    remove("tmp_qlist_4");
}

TEST_F(QListTest, CartesianCoordinatesComputed)
{
    construct_ucell(stru_lib[0]);
    ofs_running.open("tmp_qlist_cart");
    ModuleSymmetry::Symmetry symm;
    const int cal_symm_repr[2] = {0, 6};
    symm.analy_sys(ucell.lat, ucell.st, ucell.atoms, ofs_running, 1e-6, 1, "scf", cal_symm_repr);

    qlist.generate_mesh(ucell, symm, {2, 2, 2}, true);

    // after generate_mesh the Cartesian coordinates must be available and
    // consistent with the direct ones: kvec_c = kvec_d * G
    EXPECT_TRUE(qlist.kc_done);
    for (int i = 0; i < qlist.get_nq(); ++i)
    {
        ModuleBase::Vector3<double> qc = qlist.kvec_d[i] * ucell.G;
        EXPECT_DOUBLE_EQ(qlist.kvec_c[i].x, qc.x);
        EXPECT_DOUBLE_EQ(qlist.kvec_c[i].y, qc.y);
        EXPECT_DOUBLE_EQ(qlist.kvec_c[i].z, qc.z);
    }
    // Gamma is the first irreducible q-point
    EXPECT_DOUBLE_EQ(qlist.kvec_c[0].x, 0.0);
    EXPECT_DOUBLE_EQ(qlist.kvec_c[0].y, 0.0);
    EXPECT_DOUBLE_EQ(qlist.kvec_c[0].z, 0.0);

    ofs_running.close();
    ClearUcell();
    remove("tmp_qlist_cart");
}

TEST_F(QListTest, UseIrrepsSwitch)
{
    construct_ucell(stru_lib[0]);
    ofs_running.open("tmp_qlist_irreps");
    ModuleSymmetry::Symmetry symm;
    const int cal_symm_repr[2] = {0, 6};
    symm.analy_sys(ucell.lat, ucell.st, ucell.atoms, ofs_running, 1e-6, 1, "scf", cal_symm_repr);

    // use_irreps = false: the q mesh is still reduced, but no irrep data
    qlist.generate_mesh(ucell, symm, {2, 2, 2}, false);
    EXPECT_EQ(qlist.get_nq(), 4);
    EXPECT_EQ(qlist.get_nirr(0), 0); // no irrep data was computed
    EXPECT_TRUE(qlist.get_irrep_modes(0, 0).empty());

    ofs_running.close();
    ClearUcell();
    remove("tmp_qlist_irreps");
}

TEST_F(QListTest, PrintQlists)
{
    construct_ucell(stru_lib[0]);
    ofs_running.open("tmp_qlist_print");
    ModuleSymmetry::Symmetry symm;
    const int cal_symm_repr[2] = {0, 6};
    symm.analy_sys(ucell.lat, ucell.st, ucell.atoms, ofs_running, 1e-6, 1, "scf", cal_symm_repr);

    qlist.generate_mesh(ucell, symm, {1, 1, 1}, false);

    std::ofstream ofs("tmp_qlist_print_out");
    qlist.print_qlists(ofs);
    ofs.close();

    // the printed table must contain both coordinate frames
    std::ifstream ifs("tmp_qlist_print_out");
    std::string content((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_NE(content.find("Q-POINTS CARTESIAN COORDINATES"), std::string::npos);
    EXPECT_NE(content.find("Q-POINTS DIRECT COORDINATES"), std::string::npos);

    ofs_running.close();
    ClearUcell();
    remove("tmp_qlist_print");
    remove("tmp_qlist_print_out");
}

TEST_F(QListTest, ReadFromFileDirect)
{
    construct_ucell(stru_lib[0]);

    // explicit direct-coordinate q-point list with weights
    const char* fname = "tmp_qpoints_direct";
    std::ofstream ofs(fname);
    ofs << "Q_POINTS\n2\nDirect\n0.0 0.0 0.0 1.0\n0.5 0.0 0.0 1.0\n";
    ofs.close();

    qlist.read_from_file(fname, ucell);
    EXPECT_EQ(qlist.get_nq(), 2);
    EXPECT_TRUE(qlist.kd_done);
    EXPECT_TRUE(qlist.kc_done);
    EXPECT_DOUBLE_EQ(qlist.get_q(0).x, 0.0);
    EXPECT_DOUBLE_EQ(qlist.get_q(1).x, 0.5);
    // weights normalized to sum 1
    EXPECT_NEAR(qlist.wk[0] + qlist.wk[1], 1.0, 1e-10);
    // Cartesian = direct * G
    ModuleBase::Vector3<double> qc = qlist.kvec_d[1] * ucell.G;
    EXPECT_DOUBLE_EQ(qlist.kvec_c[1].x, qc.x);

    remove(fname);
    ClearUcell();
}

TEST_F(QListTest, ReadFromFileMonkhorstPack)
{
    construct_ucell(stru_lib[0]);

    const char* fname = "tmp_qpoints_mp";
    std::ofstream ofs(fname);
    ofs << "Q_POINTS\n0\nGamma\n2 2 2 0 0 0\n";
    ofs.close();

    qlist.read_from_file(fname, ucell);
    EXPECT_EQ(qlist.get_nq(), 8); // full MP mesh, no symmetry reduction here
    EXPECT_TRUE(qlist.is_mp);
    double sum = 0.0;
    for (int i = 0; i < qlist.get_nq(); ++i)
    {
        sum += qlist.wk[i];
    }
    EXPECT_NEAR(sum, 1.0, 1e-10);

    remove(fname);
    ClearUcell();
}

TEST_F(QListTest, ReadFromFileLinePath)
{
    construct_ucell(stru_lib[0]);

    const char* fname = "tmp_qpoints_line";
    std::ofstream ofs(fname);
    // G -> X segment with 4 points plus the final special point (5 total)
    ofs << "Q_POINTS\n2\nLine_Direct\n0.0 0.0 0.0 4\n0.5 0.0 0.0 1\n";
    ofs.close();

    qlist.read_from_file(fname, ucell);
    EXPECT_EQ(qlist.get_nq(), 5);
    EXPECT_TRUE(qlist.kd_done);
    EXPECT_TRUE(qlist.kc_done);
    // segment points
    EXPECT_DOUBLE_EQ(qlist.get_q(0).x, 0.0);
    EXPECT_DOUBLE_EQ(qlist.get_q(1).x, 0.5 / 4.0);
    EXPECT_DOUBLE_EQ(qlist.get_q(2).x, 2.0 * 0.5 / 4.0);
    EXPECT_DOUBLE_EQ(qlist.get_q(3).x, 3.0 * 0.5 / 4.0);
    EXPECT_DOUBLE_EQ(qlist.get_q(4).x, 0.5);
    // line weights are not normalized
    EXPECT_DOUBLE_EQ(qlist.wk[0], 1.0);

    remove(fname);
    ClearUcell();
}

TEST_F(QListTest, ReadFromFileCartesian)
{
    construct_ucell(stru_lib[0]);

    const char* fname = "tmp_qpoints_cart";
    std::ofstream ofs(fname);
    ofs << "Q_POINTS\n2\nCartesian\n0.0 0.0 0.0 1.0\n0.5 0.0 0.0 1.0\n";
    ofs.close();

    qlist.read_from_file(fname, ucell);
    EXPECT_EQ(qlist.get_nq(), 2);
    EXPECT_TRUE(qlist.kc_done);
    EXPECT_TRUE(qlist.kd_done);
    EXPECT_DOUBLE_EQ(qlist.kvec_c[1].x, 0.5);
    // direct coordinates complemented from the Cartesian ones (G = I here)
    EXPECT_DOUBLE_EQ(qlist.get_q(1).x, 0.5);
    // weights normalized to sum 1
    EXPECT_NEAR(qlist.wk[0] + qlist.wk[1], 1.0, 1e-10);

    remove(fname);
    ClearUcell();
}

TEST_F(QListTest, ReadFromFileLineCartesian)
{
    construct_ucell(stru_lib[0]);

    const char* fname = "tmp_qpoints_line_cart";
    std::ofstream ofs(fname);
    // G -> X segment with 4 points plus the final special point (5 total)
    ofs << "Q_POINTS\n2\nLine_Cartesian\n0.0 0.0 0.0 4\n0.5 0.0 0.0 1\n";
    ofs.close();

    qlist.read_from_file(fname, ucell);
    EXPECT_EQ(qlist.get_nq(), 5);
    EXPECT_TRUE(qlist.kc_done);
    EXPECT_TRUE(qlist.kd_done);
    EXPECT_DOUBLE_EQ(qlist.get_q(1).x, 0.5 / 4.0);
    EXPECT_DOUBLE_EQ(qlist.get_q(4).x, 0.5);
    // line weights are not normalized
    EXPECT_DOUBLE_EQ(qlist.wk[0], 1.0);

    remove(fname);
    ClearUcell();
}

TEST_F(QListTest, ReadFromFileNegativeNkstot)
{
    // a negative count must be rejected cleanly instead of crashing in renew()
    const char* fname = "tmp_qpoints_negative";
    std::ofstream ofs(fname);
    ofs << "Q_POINTS\n-3\nDirect\n0.0 0.0 0.0 1.0\n";
    ofs.close();

    qlist.read_from_file(fname, ucell);
    EXPECT_EQ(qlist.get_nq(), 0);

    remove(fname);
}

TEST_F(QListTest, ReadFromFileLineRejectsZeroCount)
{
    const char* fname = "tmp_qpoints_zero_count";
    std::ofstream ofs(fname);
    ofs << "Q_POINTS\n2\nLine_Direct\n0.0 0.0 0.0 0\n0.5 0.0 0.0 1\n";
    ofs.close();

    EXPECT_EXIT(qlist.read_from_file(fname, ucell), ::testing::ExitedWithCode(1), "");

    remove(fname);
}

TEST_F(QListTest, ReadFromFileMissing)
{
    // a nonexistent file yields an empty q-point list, not a crash
    qlist.read_from_file("nonexistent_qpoints", ucell);
    EXPECT_EQ(qlist.get_nq(), 0);
    // out-of-range access returns the zero vector instead of crashing
    const ModuleBase::Vector3<double> q0 = qlist.get_q(0);
    EXPECT_DOUBLE_EQ(q0.x, 0.0);
}

TEST_F(QListTest, ReadFromFileBadHeader)
{
    construct_ucell(stru_lib[0]);

    const char* fname = "tmp_qpoints_bad";
    std::ofstream ofs(fname);
    ofs << "not a q-points file\n";
    ofs.close();

    qlist.read_from_file(fname, ucell);
    EXPECT_EQ(qlist.get_nq(), 0);

    remove(fname);
    ClearUcell();
}
