#include "../sltk_atom_arrange.h"

#include <iostream>
#include <string>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "prepare_unitcell.h"

#ifdef __LCAO
InfoNonlocal::InfoNonlocal()
{
}
InfoNonlocal::~InfoNonlocal()
{
}
LCAO_Orbitals::LCAO_Orbitals()
{
}
LCAO_Orbitals::~LCAO_Orbitals()
{
}
#endif
Magnetism::Magnetism()
{
    this->tot_magnetization = 0.0;
    this->abs_magnetization = 0.0;
    this->start_magnetization = nullptr;
}
Magnetism::~Magnetism()
{
    delete[] this->start_magnetization;
}

/************************************************
 *  unit test of atom_arrange
 ***********************************************/

/**
 * - Tested Functions:
 *   - atom_arrange::delete_vector(void)
 *     - delete vector
 *   - atom_arrange::set_sr_NL
 * 	   - set the sr: search radius including nonlocal beta
 */

void SetGlobalV()
{
    GlobalV::test_grid = 0;
    GlobalV::test_grid_driver = 0;
    GlobalV::test_deconstructor = 0;
}

class SltkAtomArrangeTest : public testing::Test
{
  protected:
    UnitCell* ucell;
    UcellTestPrepare utp = UcellTestLib["Si"];
    std::ofstream ofs;
    std::ifstream ifs;
    bool pbc = true;
    double radius = ((8 + 5.01) * 2.0 + 0.01) / 10.2;
    int test_atom_in = 0;
    std::string output;
    void SetUp()
    {
        SetGlobalV();
        ucell = utp.SetUcellInfo();
    }
    void TearDown()
    {
        delete ucell;
    }
};

TEST_F(SltkAtomArrangeTest, setsrNL)
{
    atom_arrange test;
    std::string teststring = "m";
    double rcutmax_Phi = 1;
    double rcutmax_Beta = 2;
    bool gamma_only_local = true;
    double test_sr = 0;
    std::ofstream ofs;
    ofs.open("./to_test_arrange.txt");
    test_sr = test.set_sr_NL(ofs, teststring, rcutmax_Phi, rcutmax_Beta, gamma_only_local);
    EXPECT_DOUBLE_EQ(test_sr, 2.01);

    gamma_only_local = false;
    test_sr = test.set_sr_NL(ofs, teststring, rcutmax_Phi, rcutmax_Beta, gamma_only_local);
    EXPECT_DOUBLE_EQ(test_sr, 6.01);

    teststring = "no";
    test_sr = test.set_sr_NL(ofs, teststring, rcutmax_Phi, rcutmax_Beta, gamma_only_local);
    ofs.close();
    std::ifstream ifs;
    std::string test2, s;
    ifs.open("./to_test_arrange.txt");
    std::string str((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_THAT(str, testing::HasSubstr("longest orb rcut (Bohr) = 1"));
    EXPECT_THAT(str, testing::HasSubstr("longest nonlocal projector rcut (Bohr) = 2"));
    ifs.close();
    remove("./to_test_arrange");
}

TEST_F(SltkAtomArrangeTest, Search)
{
    ucell->check_dtau();
    Grid_Driver grid_d(GlobalV::test_deconstructor, GlobalV::test_grid_driver, GlobalV::test_grid);
    ofs.open("test.out");
    bool test_only = true;
    atom_arrange::search(pbc, ofs, grid_d, *ucell, radius, test_atom_in, test_only);
    EXPECT_EQ(grid_d.getType(0),0);
    EXPECT_EQ(grid_d.getNatom(0), 1); // adjacent atom is 1
    ofs.close();
    ifs.open("test.out");
    std::string str((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_THAT(str, testing::HasSubstr("search neighboring atoms done."));
    remove("test.out");
}

int main(int argc, char **argv)
{
        MPI_Init(&argc, &argv);
        testing::InitGoogleTest(&argc, argv);

        int nproc;
        int my_rank;

        MPI_Comm_size(MPI_COMM_WORLD,&nproc);
        MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

        int result = RUN_ALL_TESTS();

        MPI_Finalize();

        return result;
}