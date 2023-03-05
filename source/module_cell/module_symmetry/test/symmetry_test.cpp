#include "module_base/mathzone.h"
#include "../symmetry.h"
#include "module_base/matrix3.h"
#include "module_base/vector3.h"
#include "mpi.h"

#include "gtest/gtest.h"

/************************************************
 *  unit test of class Symmetry
 ***********************************************/

// mock the useless functions
void output::printM3(std::ofstream &ofs, const std::string &description, const ModuleBase::Matrix3 &m)
{
}
pseudo_nc::pseudo_nc()
{
}
pseudo_nc::~pseudo_nc()
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

/**
    switch(ibrav)
    {
        case 1: return "01. Cubic P (simple)";
        case 2: return "02. Cubic I (body-centered)";
        case 3: return "03. Cubic F (face-centered)";
        case 4: return "04. Hexagonal cell";
        case 5: return "05. Tetrogonal P (simple)";
        case 6: return "06. Tetrogonal I (body-centered)";
        case 7: return "07. Rhombohedral (Trigonal) cell";
        case 8: return "08. Orthorhombic P(simple)";
        case 9: return "09. Orthorhombic I (body-centered)";
        case 10: return "10. Orthorhombic F (face-centered)";
        case 11: return "11. Orthorhombic C (base-centered)";
        case 12: return "12. Monoclinic P (simple)";
        case 13: return "13. Monoclinic A (base-center)";
        case 14: return "14. Triclinic cell";
        case 15: return "wrong !! ";
    }
*/

struct atomtype_
{
    std::string atomname;
    std::vector<std::vector<double>> coordinate;
};

struct stru_
{
    int ibrav;
    std::string point_group; // Schoenflies symbol
    std::string point_group_hm; // Hermann–Mauguin notation.
    std::string space_group;
    std::vector<double> cell;
    std::vector<atomtype_> all_type;
};

std::vector<stru_> stru_lib{
    stru_{1,
          "O_h",
          "m-3m",
          "Pm-3m",
          std::vector<double>{1., 0., 0., 0., 1., 0., 0., 0., 1.},
          std::vector<atomtype_>{atomtype_{"C",
                                           std::vector<std::vector<double>>{
                                               {0., 0., 0.},
                                           }}}},
    stru_{2,
          "O_h",
          "m-3m",
          "Im-3m",
          std::vector<double>{-0.5, 0.5, 0.5, 0.5, -0.5, 0.5, 0.5, 0.5, -0.5},
          std::vector<atomtype_>{atomtype_{"C",
                                           std::vector<std::vector<double>>{
                                               {0., 0., 0.},
                                           }}}},
    stru_{3,
          "O_h",
          "m-3m",
          "Fm-3m",
          std::vector<double>{0., 0.5, 0.5, 0.5, 0., 0.5, 0.5, 0.5, 0.},
          std::vector<atomtype_>{atomtype_{"C",
                                           std::vector<std::vector<double>>{
                                               {0., 0., 0.},
                                           }}}},
    stru_{4,
          "D_6h",
          "6/mmm",
          "P6/mmm",
          std::vector<double>{1., 0., 0., -0.5, 0.8660254, 0., 0., 0., 2.},
          std::vector<atomtype_>{atomtype_{"C",
                                           std::vector<std::vector<double>>{
                                               {0., 0., 0.},
                                           }}}},
    stru_{5,
          "D_4h",
          "4/mmm",
          "P4/mmm",
          std::vector<double>{1., 0., 0., 0., 1., 0., 0., 0., 2.},
          std::vector<atomtype_>{atomtype_{"C",
                                           std::vector<std::vector<double>>{
                                               {0., 0., 0.},
                                           }}}},
    stru_{6,
          "D_4h",
          "4/mmm",
          "I4/mmm",
          std::vector<double>{-0.35355339, 0.35355339, 1., 0.35355339, -0.35355339, 1., 0.35355339, 0.35355339, -1.},
          std::vector<atomtype_>{atomtype_{"C",
                                           std::vector<std::vector<double>>{
                                               {0., 0., 0.},
                                           }}}},
    stru_{7,
          "D_3d",
          "-3m",
          "R-3m",
          std::vector<double>{0.57357644,
                              0.33115451,
                              0.74923078,
                              -0.57357644,
                              0.33115451,
                              0.74923078,
                              0.,
                              -0.66230902,
                              0.74923078},
          std::vector<atomtype_>{atomtype_{"C",
                                           std::vector<std::vector<double>>{
                                               {-0., 0., 0.},
                                           }}}},
    stru_{8,
          "D_2h",
          "mmm",
          "Pmmm",
          std::vector<double>{1., 0., 0., 0., 2., 0., 0., 0., 3.},
          std::vector<atomtype_>{atomtype_{"C",
                                           std::vector<std::vector<double>>{
                                               {0., 0., 0.},
                                           }}}},
    stru_{9,
          "D_2h",
          "mmm",
          "Immm",
          std::vector<double>{-0.25, 0.75, 1., 0.25, -0.75, 1., 0.25, 0.75, -1.},
          std::vector<atomtype_>{atomtype_{"C",
                                           std::vector<std::vector<double>>{
                                               {0., 0., 0.},
                                           }}}},
    stru_{10,
          "D_2h",
          "mmm",
          "Fmmm",
          std::vector<double>{0., 1., 1.5, 0.5, 0., 1.5, 0.5, 1., 0.},
          std::vector<atomtype_>{atomtype_{"C",
                                           std::vector<std::vector<double>>{
                                               {0., 0., 0.},
                                           }}}},
    stru_{11,
          "D_2h",
          "mmm",
          "Cmmm",
          std::vector<double>{0.5, -1.5, 0., 0.5, 1.5, 0., 0., 0., 2.},
          std::vector<atomtype_>{atomtype_{"C",
                                           std::vector<std::vector<double>>{
                                               {0., 0., 0.},
                                           }}}},
    stru_{12,
          "C_2h",
          "2/m",
          "P2/m",
          std::vector<double>{1., 0., 0., 0., 2., 0., -0.02606043, 0., 2.81907786},
          std::vector<atomtype_>{atomtype_{"C",
                                           std::vector<std::vector<double>>{
                                               {0., 0., 0.},
                                           }}}},
//    stru_{13,
//          "C_2h",
//          "2/m",
//          "C2/m",
//          std::vector<double>{0.5, -1., 0., 0.5, 1., 0., -0.40192379, 0., 1.5},
//          std::vector<atomtype_>{atomtype_{"C",
//                                           std::vector<std::vector<double>>{
//                                               {0., 0., 0.},
//                                           }}}},
//    stru_{14,
//          "C_i",
//          "-1",
//          "P-1",
//          std::vector<double>{1., 0., 0., -0.28989928, 1.53691386, 0., -0.31595971, -0.66789914, 1.75670135},
//          std::vector<atomtype_>{atomtype_{"C",
//                                           std::vector<std::vector<double>>{
//                                               {0., 0., 0.},
//                                           }}}},

};

class SymmetryTest : public testing::Test
{
  protected:
    UnitCell ucell;
    std::ofstream ofs_running;

    void construct_ucell(stru_ &stru)
    {
        std::vector<atomtype_> coord = stru.all_type;
        ucell.a1 = ModuleBase::Vector3<double>(stru.cell[0], stru.cell[1], stru.cell[2]);
        ucell.a2 = ModuleBase::Vector3<double>(stru.cell[3], stru.cell[4], stru.cell[5]);
        ucell.a3 = ModuleBase::Vector3<double>(stru.cell[6], stru.cell[7], stru.cell[8]);
        ucell.ntype = stru.all_type.size();
        ucell.atoms = new Atom[ucell.ntype];
        ucell.nat = 0;

        for (int i = 0; i < coord.size(); i++)
        {
            ucell.atoms[i].label = coord[i].atomname;
            ucell.atoms[i].na = coord[i].coordinate.size();
            ucell.atoms[i].tau = new ModuleBase::Vector3<double>[ucell.atoms[i].na];
            ucell.atoms[i].taud = new ModuleBase::Vector3<double>[ucell.atoms[i].na];
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
        for (int i = 0; i < ucell.ntype; i++)
        {
            delete[] ucell.atoms[i].tau;
            delete[] ucell.atoms[i].taud;
        }
        delete[] ucell.atoms;
    }
};

TEST_F(SymmetryTest, AnalySys)
{
    for (int i = 0; i < stru_lib.size(); i++)
    {
        ModuleSymmetry::Symmetry symm;
        construct_ucell(stru_lib[i]);
        symm.analy_sys(ucell, ofs_running);
        std::string ref_point_group = stru_lib[i].point_group;
        std::string cal_point_group = symm.pgname;
        int ref_ibrav = stru_lib[i].ibrav;
        int cal_ibrav = symm.real_brav;
        EXPECT_EQ(cal_ibrav, ref_ibrav);
        EXPECT_EQ(cal_point_group, ref_point_group) << "ibrav=" << stru_lib[i].ibrav;
        ClearUcell();
    }
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();
    MPI_Finalize();
    return result;
}
