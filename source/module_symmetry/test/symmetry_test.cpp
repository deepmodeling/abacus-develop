#include "module_base/mathzone.h"
#include "../symmetry.h"
#include "module_base/matrix3.h"
#include "module_base/vector3.h"
#include "mpi.h"

#include "gtest/gtest.h"

#define DOUBLETHRESHOLD 1e-8

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

        //1. ibrav
        std::string ref_point_group = stru_lib[i].point_group;
        std::string cal_point_group = symm.pgname;
        int ref_ibrav = stru_lib[i].ibrav;
        int cal_ibrav = symm.real_brav;
        EXPECT_EQ(cal_ibrav, ref_ibrav);
        EXPECT_EQ(cal_point_group, ref_point_group) << "ibrav=" << stru_lib[i].ibrav;
        
        //2. input and optimized lattice, gtrans_convert and veccon 
        //input lattice
        EXPECT_EQ(symm.s1, ucell.a1);
        EXPECT_EQ(symm.s2, ucell.a2);
        EXPECT_EQ(symm.s3, ucell.a3);
        //optimized lattice
        EXPECT_EQ(symm.a1, ModuleBase::Vector3<double>(symm.optlat.e11, symm.optlat.e12, symm.optlat.e13));
        EXPECT_EQ(symm.a2, ModuleBase::Vector3<double>(symm.optlat.e21, symm.optlat.e22, symm.optlat.e23));
        EXPECT_EQ(symm.a3, ModuleBase::Vector3<double>(symm.optlat.e31, symm.optlat.e32, symm.optlat.e33));
        //gtrans_convert
        std::vector<ModuleBase::Vector3<double>> gtrans_optconf(symm.nrotk);
        double* gtrans_veccon=new double [symm.nrotk*3];
        for (int i=0;i<symm.nrotk;++i)
        {
            gtrans_veccon[3*i]=symm.gtrans[i].x;
            gtrans_veccon[3*i+1]=symm.gtrans[i].y;
            gtrans_veccon[3*i+2]=symm.gtrans[i].z;
        }
        double* gtrans_optconf_veccon=new double [symm.nrotk*3];
        symm.gtrans_convert(symm.gtrans, gtrans_optconf.data(), symm.nrotk, ucell.latvec, symm.optlat);
        symm.veccon(gtrans_veccon, gtrans_optconf_veccon, symm.nrotk, symm.s1, symm.s2, symm.s3, symm.a1, symm.a2, symm.a3);
        for(int i=0;i<symm.nrotk;++i)
            EXPECT_EQ(gtrans_optconf[i], ModuleBase::Vector3<double>(gtrans_optconf_veccon[i*3], 
            gtrans_optconf_veccon[i*3+1], gtrans_optconf_veccon[i*3+2]));
        delete[] gtrans_veccon;
        delete[] gtrans_optconf_veccon;

        //3. invmap
        int* ivmp=new int[symm.nrotk];
        symm.gmatrix_invmap(symm.gmatrix, symm.nrotk, ivmp);
        ModuleBase::Matrix3 test;

        for (int i=0;i<symm.nrotk;++i)
        {
            test=symm.gmatrix[i]*symm.gmatrix[ivmp[i]];
            EXPECT_NEAR(test.e11, 1, DOUBLETHRESHOLD);
            EXPECT_NEAR(test.e22, 1, DOUBLETHRESHOLD);
            EXPECT_NEAR(test.e33, 1, DOUBLETHRESHOLD);
            EXPECT_NEAR(test.e12, 0, DOUBLETHRESHOLD);
            EXPECT_NEAR(test.e21, 0, DOUBLETHRESHOLD);
            EXPECT_NEAR(test.e13, 0, DOUBLETHRESHOLD);
            EXPECT_NEAR(test.e31, 0, DOUBLETHRESHOLD);
            EXPECT_NEAR(test.e23, 0, DOUBLETHRESHOLD);
            EXPECT_NEAR(test.e32, 0, DOUBLETHRESHOLD);

        }
        delete[] ivmp;

        //4. gmatrix_convert : input <-> opt <-> reciprocal
        ModuleBase::Matrix3* gmatrix_inputconf=new ModuleBase::Matrix3[symm.nrotk];
        ModuleBase::Matrix3* gmatrix_optconf=new ModuleBase::Matrix3[symm.nrotk];
        ModuleBase::Matrix3* kgmatrix=new ModuleBase::Matrix3[symm.nrotk];
        symm.gmatrix_convert(symm.gmatrix, gmatrix_optconf, symm.nrotk, ucell.latvec, symm.optlat);
        symm.gmatrix_convert(gmatrix_optconf, gmatrix_inputconf, symm.nrotk, symm.optlat, ucell.latvec);
        symm.gmatrix_convert(gmatrix_optconf, kgmatrix, symm.nrotk, symm.optlat, ucell.G);
        for (int i=0;i<symm.nrotk;++i)
        {

            EXPECT_NEAR(symm.gmatrix[i].e11, gmatrix_inputconf[i].e11, DOUBLETHRESHOLD);
            EXPECT_NEAR(symm.gmatrix[i].e22, gmatrix_inputconf[i].e22, DOUBLETHRESHOLD);
            EXPECT_NEAR(symm.gmatrix[i].e33, gmatrix_inputconf[i].e33, DOUBLETHRESHOLD);
            EXPECT_NEAR(symm.gmatrix[i].e12, gmatrix_inputconf[i].e12, DOUBLETHRESHOLD);
            EXPECT_NEAR(symm.gmatrix[i].e21, gmatrix_inputconf[i].e21, DOUBLETHRESHOLD);
            EXPECT_NEAR(symm.gmatrix[i].e13, gmatrix_inputconf[i].e13, DOUBLETHRESHOLD);
            EXPECT_NEAR(symm.gmatrix[i].e31, gmatrix_inputconf[i].e31, DOUBLETHRESHOLD);
            EXPECT_NEAR(symm.gmatrix[i].e23, gmatrix_inputconf[i].e23, DOUBLETHRESHOLD);
            EXPECT_NEAR(symm.gmatrix[i].e32, gmatrix_inputconf[i].e32, DOUBLETHRESHOLD);
            
            ModuleBase::Matrix3 tmpA=ucell.latvec*ucell.latvec*kgmatrix[i];
            ModuleBase::Matrix3 tmpB=symm.gmatrix[i]*ucell.latvec*ucell.latvec;
            EXPECT_NEAR(tmpA.e11, tmpB.e11, DOUBLETHRESHOLD);
            EXPECT_NEAR(tmpA.e22, tmpB.e22, DOUBLETHRESHOLD);
            EXPECT_NEAR(tmpA.e33, tmpB.e33, DOUBLETHRESHOLD);
            EXPECT_NEAR(tmpA.e12, tmpB.e12, DOUBLETHRESHOLD);
            EXPECT_NEAR(tmpA.e21, tmpB.e21, DOUBLETHRESHOLD);
            EXPECT_NEAR(tmpA.e13, tmpB.e13, DOUBLETHRESHOLD);
            EXPECT_NEAR(tmpA.e31, tmpB.e31, DOUBLETHRESHOLD);
            EXPECT_NEAR(tmpA.e23, tmpB.e23, DOUBLETHRESHOLD);
            EXPECT_NEAR(tmpA.e32, tmpB.e32, DOUBLETHRESHOLD);
        }
            
        delete[] gmatrix_inputconf;
        delete[] gmatrix_optconf;
        delete[] kgmatrix;

        ClearUcell();
    }
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();
    MPI_Finalize();
    return 0;
}
