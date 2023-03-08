#include "module_base/mathzone.h"
#include "../symmetry.h"
#include "module_base/matrix3.h"
#include "module_base/vector3.h"
#include "mpi.h"

#include "gtest/gtest.h"

#define DOUBLETHRESHOLD 1e-8

/************************************************
 *  unit test of class Symmetry
 * 1. function: `analy_sys`:
 * 1-1. test if the bravis lattice analysis is right;
 * 1-2. check if matrix-type and vector3-type;
 *     input and optimized lattice vectors are right;
 * 1-3. double-check for  if `gtrans_convert` 
 *     gives the same results as `veccon`;
 * 1-4. check `invmap` function gives the right result;
 * 1-5 test if `gmatrix_convert` and `gmatrix_convert_int`
 *     gives the right result;
 * 2. function: `atom_ordering_new:
 * test the new atom-sort algorithm gives the right result;
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
   stru_{13,
         "C_2h",
         "2/m",
         "C2/m",
         std::vector<double>{0.5, -1., 0., 0.5, 1., 0., -0.40192379, 0., 1.5},
         std::vector<atomtype_>{atomtype_{"C",
                                          std::vector<std::vector<double>>{
                                              {0., 0., 0.},
                                          }}}},
   stru_{14,
         "S_2",
         "-1",
         "P-1",
         std::vector<double>{1., 0., 0., -0.28989928, 1.53691386, 0., -0.31595971, -0.66789914, 1.75670135},
         std::vector<atomtype_>{atomtype_{"C",
                                          std::vector<std::vector<double>>{
                                              {0., 0., 0.},
                                          }}}},

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
        ucell.latvec.e11=ucell.a1.x;
        ucell.latvec.e12=ucell.a1.y;
        ucell.latvec.e13=ucell.a1.z;
        ucell.latvec.e21=ucell.a2.x;
        ucell.latvec.e22=ucell.a2.y;
        ucell.latvec.e23=ucell.a2.z;
        ucell.latvec.e31=ucell.a3.x;
        ucell.latvec.e32=ucell.a3.y;
        ucell.latvec.e33=ucell.a3.z;
        ucell.GT = ucell.latvec.Inverse();
        ucell.G = ucell.GT.Transpose();
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
    for (int stru = 0; stru < stru_lib.size(); stru++)
    {
        ModuleSymmetry::Symmetry symm;
        construct_ucell(stru_lib[stru]);
        symm.analy_sys(ucell, ofs_running);

        //1. ibrav
        std::string ref_point_group = stru_lib[stru].point_group;
        std::string cal_point_group = symm.pgname;
        int ref_ibrav = stru_lib[stru].ibrav;
        int cal_ibrav = symm.real_brav;
        EXPECT_EQ(cal_ibrav, ref_ibrav);
        EXPECT_EQ(cal_point_group, ref_point_group) << "ibrav=" << stru_lib[stru].ibrav;
        
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

        //4. gmatrix_convert : input(gmatrix) -> opt(gmatrix_opt) ->input(gmatrix_input_back)
        //-> opt(gmatrix_optback)   <-> reciprocal(int or non-int)
        ModuleBase::Matrix3* gmatrix_input_back=new ModuleBase::Matrix3[symm.nrotk];//3
        ModuleBase::Matrix3* gmatrix_opt=new ModuleBase::Matrix3[symm.nrotk];//2
        ModuleBase::Matrix3* gmatrix_opt_back=new ModuleBase::Matrix3[symm.nrotk];//4
        ModuleBase::Matrix3* kgmatrix_nonint=new ModuleBase::Matrix3[symm.nrotk];
        symm.gmatrix_convert_int(symm.gmatrix, gmatrix_opt, symm.nrotk, ucell.latvec, symm.optlat); //1->2
        symm.gmatrix_convert_int(gmatrix_opt, gmatrix_input_back, symm.nrotk, symm.optlat, ucell.latvec);   //2->3
        symm.gmatrix_convert_int(gmatrix_input_back, gmatrix_opt_back, symm.nrotk, ucell.latvec, symm.optlat); //3->4
        
        symm.gmatrix_convert(symm.gmatrix, kgmatrix_nonint, symm.nrotk, ucell.latvec, ucell.G);
        for (int i=0;i<symm.nrotk;++i)
        {
            EXPECT_NEAR(symm.gmatrix[i].e11, gmatrix_input_back[i].e11, DOUBLETHRESHOLD);
            EXPECT_NEAR(symm.gmatrix[i].e22, gmatrix_input_back[i].e22, DOUBLETHRESHOLD);
            EXPECT_NEAR(symm.gmatrix[i].e33, gmatrix_input_back[i].e33, DOUBLETHRESHOLD);
            EXPECT_NEAR(symm.gmatrix[i].e12, gmatrix_input_back[i].e12, DOUBLETHRESHOLD);
            EXPECT_NEAR(symm.gmatrix[i].e21, gmatrix_input_back[i].e21, DOUBLETHRESHOLD);
            EXPECT_NEAR(symm.gmatrix[i].e13, gmatrix_input_back[i].e13, DOUBLETHRESHOLD);
            EXPECT_NEAR(symm.gmatrix[i].e31, gmatrix_input_back[i].e31, DOUBLETHRESHOLD);
            EXPECT_NEAR(symm.gmatrix[i].e23, gmatrix_input_back[i].e23, DOUBLETHRESHOLD);
            EXPECT_NEAR(symm.gmatrix[i].e32, gmatrix_input_back[i].e32, DOUBLETHRESHOLD);
            EXPECT_NEAR(gmatrix_opt[i].e11, gmatrix_opt_back[i].e11, DOUBLETHRESHOLD);
            EXPECT_NEAR(gmatrix_opt[i].e22, gmatrix_opt_back[i].e22, DOUBLETHRESHOLD);
            EXPECT_NEAR(gmatrix_opt[i].e33, gmatrix_opt_back[i].e33, DOUBLETHRESHOLD);
            EXPECT_NEAR(gmatrix_opt[i].e12, gmatrix_opt_back[i].e12, DOUBLETHRESHOLD);
            EXPECT_NEAR(gmatrix_opt[i].e21, gmatrix_opt_back[i].e21, DOUBLETHRESHOLD);
            EXPECT_NEAR(gmatrix_opt[i].e13, gmatrix_opt_back[i].e13, DOUBLETHRESHOLD);
            EXPECT_NEAR(gmatrix_opt[i].e31, gmatrix_opt_back[i].e31, DOUBLETHRESHOLD);
            EXPECT_NEAR(gmatrix_opt[i].e23, gmatrix_opt_back[i].e23, DOUBLETHRESHOLD);
            EXPECT_NEAR(gmatrix_opt[i].e32, gmatrix_opt_back[i].e32, DOUBLETHRESHOLD);
            
            ModuleBase::Matrix3 tmpA=symm.optlat.Inverse()*gmatrix_opt[i]*symm.optlat; //A^-1*SA*A
            ModuleBase::Matrix3 tmpB=ucell.latvec.Inverse()*symm.gmatrix[i]*ucell.latvec;//B^-1*SB*B
            ModuleBase::Matrix3 tmpG_int=ucell.G.Inverse()*symm.kgmatrix[i]*ucell.G;//G^-1*SG*G
            ModuleBase::Matrix3 tmpG=ucell.G.Inverse()*kgmatrix_nonint[i]*ucell.G;//G^-1*SG*G
            EXPECT_NEAR(tmpA.e11, tmpB.e11, DOUBLETHRESHOLD);
            EXPECT_NEAR(tmpA.e22, tmpB.e22, DOUBLETHRESHOLD);
            EXPECT_NEAR(tmpA.e33, tmpB.e33, DOUBLETHRESHOLD);
            EXPECT_NEAR(tmpA.e12, tmpB.e12, DOUBLETHRESHOLD);
            EXPECT_NEAR(tmpA.e21, tmpB.e21, DOUBLETHRESHOLD);
            EXPECT_NEAR(tmpA.e13, tmpB.e13, DOUBLETHRESHOLD);
            EXPECT_NEAR(tmpA.e31, tmpB.e31, DOUBLETHRESHOLD);
            EXPECT_NEAR(tmpA.e23, tmpB.e23, DOUBLETHRESHOLD);
            EXPECT_NEAR(tmpA.e32, tmpB.e32, DOUBLETHRESHOLD);

            if(!symm.equal(tmpG.e13, tmpG_int.e13) || !symm.equal(tmpG.e23, tmpG_int.e23) || !symm.equal(tmpG.e12, tmpG_int.e12))
            {
                std::cout<<"stru_ibrav:"<<stru_lib[stru].ibrav<<std::endl;
                std::cout<<"isymm: "<<i<<std::endl;
                std::cout<<"kgmatrix[i], int:"<<std::endl;
                std::cout<<symm.kgmatrix[i].e11<<" "<<symm.kgmatrix[i].e12<<" "<<symm.kgmatrix[i].e13<<std::endl;
                std::cout<<symm.kgmatrix[i].e21<<" "<<symm.kgmatrix[i].e22<<" "<<symm.kgmatrix[i].e23<<std::endl;
                std::cout<<symm.kgmatrix[i].e31<<" "<<symm.kgmatrix[i].e32<<" "<<symm.kgmatrix[i].e33<<std::endl;
                std::cout<<"kgmatrix[i], nonint:"<<std::endl;
                std::cout<<kgmatrix_nonint[i].e11<<" "<<kgmatrix_nonint[i].e12<<" "<<kgmatrix_nonint[i].e13<<std::endl;
                std::cout<<kgmatrix_nonint[i].e21<<" "<<kgmatrix_nonint[i].e22<<" "<<kgmatrix_nonint[i].e23<<std::endl;
                std::cout<<kgmatrix_nonint[i].e31<<" "<<kgmatrix_nonint[i].e32<<" "<<kgmatrix_nonint[i].e33<<std::endl;
            }
            EXPECT_NEAR(tmpA.e11, tmpG.e11, DOUBLETHRESHOLD);
            EXPECT_NEAR(tmpA.e22, tmpG.e22, DOUBLETHRESHOLD);
            EXPECT_NEAR(tmpA.e33, tmpG.e33, DOUBLETHRESHOLD);
            EXPECT_NEAR(tmpA.e12, tmpG.e12, DOUBLETHRESHOLD);
            EXPECT_NEAR(tmpA.e21, tmpG.e21, DOUBLETHRESHOLD);
            EXPECT_NEAR(tmpA.e13, tmpG.e13, DOUBLETHRESHOLD);
            EXPECT_NEAR(tmpA.e31, tmpG.e31, DOUBLETHRESHOLD);
            EXPECT_NEAR(tmpA.e23, tmpG.e23, DOUBLETHRESHOLD);
            EXPECT_NEAR(tmpA.e32, tmpG.e32, DOUBLETHRESHOLD);
        }
            
        delete[] gmatrix_input_back;
        delete[] gmatrix_opt;
        delete[] gmatrix_opt_back;
        delete[] kgmatrix_nonint;

        ClearUcell();
    }
}

TEST_F(SymmetryTest, AtomOrderingNew)
{
    // the old function `atom_ordering` has bugs 
    // so here I do not compare with its results
    ModuleSymmetry::Symmetry symm;
    symm.epsilon=1e-5;
    int nat=20; //number of atoms
    double* old_pos=new double[nat*3];
    double* new_pos=new double[nat*3];
    int* subindex=new int[nat];
    //generate random number and restrict to [-0.5,0.5)
    srand(time(NULL));
    for(int i=0;i<3*nat;++i)
    {
        old_pos[i]=double(rand())/double(RAND_MAX)-0.5;
        new_pos[i]=old_pos[i];
    }
    //ordering
    symm.test_atom_ordering(new_pos, nat, subindex);
    //check 
    for (int i=0;i<nat-1;++i)
    {
        //x[i]<=x[i+1]
        EXPECT_LE(new_pos[3*i], new_pos[3*(i+1)] +symm.epsilon);
        if (symm.equal(new_pos[3*i], new_pos[3*(i+1)]))
        {   //y[i]<=y[i+1]
            EXPECT_LE(new_pos[3*i+1], new_pos[3*(i+1)+1] + symm.epsilon);
            if(symm.equal(new_pos[3*i+1], new_pos[3*(i+1)+1]))
            {   //z[i]<=z[i+1]
                EXPECT_LE(new_pos[3*i+2], new_pos[3*(i+1)+2] +symm.epsilon);
            }
        }
    }
    //output old_pos and new_pos
    std::cout<<"random direct coords: "<<std::endl;
    std::cout<<"iatom     x      y     z"<<std::endl;
    for (int i=0;i<nat;++i)
        std::cout<<i<<" "<<old_pos[3*i]<<" "<<old_pos[3*i+1]<<" "<<old_pos[3*i+2]<<std::endl;
    std::cout<<"sorted direct coords: "<<std::endl;
    std::cout<<"iatom     x      y     z"<<std::endl;
    for (int i=0;i<nat;++i)
        std::cout<<i<<" "<<new_pos[3*i]<<" "<<new_pos[3*i+1]<<" "<<new_pos[3*i+2]<<std::endl;
    delete[] old_pos;
    delete[] new_pos;
    delete[] subindex;

}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();
    MPI_Finalize();
    return result;
}
