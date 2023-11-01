#include "symmetry_test_cases.h"
#include "mpi.h"
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
 *3. function: `pricell`:
 * test if the number of primitive cells are right, 
 * using cases whose space group
 * is different from its point group.
 ***********************************************/
// mock the useless functions
void output::printM3(std::ofstream &ofs, const std::string &description, const ModuleBase::Matrix3 &m){}
pseudo::pseudo()
{
}
pseudo::~pseudo()
{
}
Atom::Atom(){}
Atom::~Atom(){}
Atom_pseudo::Atom_pseudo(){}
Atom_pseudo::~Atom_pseudo(){}
UnitCell::UnitCell(){}
UnitCell::~UnitCell(){}
Magnetism::Magnetism(){}
Magnetism::~Magnetism() {}

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

    // 1. Random Test
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

    //2. Special Case Test
    //(1). z-to-x
    new_pos[0]=0.0; new_pos[1]=1./3.; new_pos[2]=0.1;
    new_pos[3]=symm.epsilon*0.9; new_pos[4]=1./3.; new_pos[5]=0.0;
    symm.test_atom_ordering(new_pos, 2, subindex);
    EXPECT_NEAR(new_pos[5], 0.1, DOUBLETHRESHOLD);
    EXPECT_NEAR(new_pos[2], 0.0, DOUBLETHRESHOLD);

    delete[] old_pos;
    delete[] new_pos;
    delete[] subindex;

}

TEST_F(SymmetryTest, SG_Pricell)
{
    for (int stru = 0; stru < supercell_lib.size(); stru++)
    {
        ModuleSymmetry::Symmetry symm;
        symm.epsilon = 1e-5;
        construct_ucell(supercell_lib[stru]);
        symm.analy_sys(ucell, ofs_running);

        std::string ref_point_group = supercell_lib[stru].point_group;
        std::string cal_point_group = symm.pgname;
        std::string ref_space_group = supercell_lib[stru].space_group;
        std::string cal_space_group = symm.spgname;
        
        int ref_ncells = supercell_lib[stru].ibrav;
        EXPECT_EQ(symm.ncell, ref_ncells);
        EXPECT_EQ(cal_point_group, ref_point_group);
        EXPECT_EQ(cal_space_group, ref_space_group);

        ClearUcell();
    }
}

TEST_F(SymmetryTest, SubGroup)
{
    ModuleSymmetry::Symmetry symm;
    EXPECT_EQ(symm.subgroup(30, 1, 6, 8, 0, 0, 6, 0, 0, 8), 29);//24, T_h
    EXPECT_EQ(symm.subgroup(32, 1, 9, 0, 6, 0, 9, 0, 6, 0), 20);//16, D_4h
    EXPECT_EQ(symm.subgroup(17, 0, 5, 2, 0, 2, 7, 0, 0, 0), 25);//12, C_6v
    EXPECT_EQ(symm.subgroup(10, 1, 3, 0, 1, 0, 3, 0, 1, 0), 8);//8, D_2h
    EXPECT_EQ(symm.subgroup(10, 1, 1, 2, 0, 0, 3, 1, 0, 1), 12);//6, C_3v
    EXPECT_EQ(symm.subgroup(7, 1, 1, 0, 1, 0, 1, 0, 2, 0), 15);//4, S_4
    EXPECT_EQ(symm.subgroup(5, 0, 2, 2, 0, 0, 0, 0, 0, 0), 9);///3, C_3
    EXPECT_EQ(symm.subgroup(4, 0, 0, 1, 0, 0, 1, 1, 0, 0), 4);//C_1h
    EXPECT_EQ(symm.subgroup(4, 0, 0, 1, 0, 0, 0, 2, 0, 0), 1);//C_1
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();
    MPI_Finalize();
    return result;
}
