#include <gtest/gtest.h>
#include "module_io/orb_io.h"

#ifdef __MPI
#include <mpi.h>
#endif

#include "module_base/constants.h"
#include "module_base/global_variable.h"

class OrbIOTest : public testing::Test
{
  protected:
    void SetUp();
    void TearDown(){};

    const std::string file = "../../../../../tests/PP_ORB/Ti_gga_10au_100Ry_4s2p2d1f.orb";
    const double tol = 1e-12;
};

void OrbIOTest::SetUp()
{
#ifdef __MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &GlobalV::MY_RANK);
#endif
}

TEST_F(OrbIOTest, ReadAbacusOrb)
{
    std::ifstream ifs(file);
    std::string elem;
    double ecut, dr;
    int nr;
    std::vector<int> nzeta;
    std::vector<std::vector<double>> radials;
    ModuleIO::read_abacus_orb(ifs, elem, ecut, nr, dr, nzeta, radials);
    EXPECT_EQ(elem, "Ti");
    EXPECT_DOUBLE_EQ(ecut, 100.0);
    EXPECT_EQ(nr, 1001);
    EXPECT_DOUBLE_EQ(dr, 0.01);
    EXPECT_EQ(nzeta.size(), 4); // l from 0 to 3
    EXPECT_EQ(nzeta[0], 4);
    EXPECT_EQ(nzeta[1], 2);
    EXPECT_EQ(nzeta[2], 2);
    EXPECT_EQ(nzeta[3], 1);
    EXPECT_EQ(radials.size(), 9); // 4 + 2 + 2 + 1
    for(auto& radial: radials)
    {
        EXPECT_EQ(radial.size(), 1001);
    }
    EXPECT_EQ(radials[0][0], -1.581711853170e-01);
    EXPECT_EQ(radials[0][4], -1.583907030513e-01);
    EXPECT_EQ(radials[0][996], -4.183526380009e-05);
    EXPECT_EQ(radials[0][1000], 0);
    EXPECT_EQ(radials[3][0], -1.166292682541e+00);
    EXPECT_EQ(radials[3][4], -1.164223359672e+00);
    EXPECT_EQ(radials[3][996], -3.183325576529e-04);
    EXPECT_EQ(radials[3][1000], 0);
    EXPECT_EQ(radials[8][0], 0);
    EXPECT_EQ(radials[8][4], 3.744878535962e-05);
    EXPECT_EQ(radials[8][996], 7.495357740660e-05);
    EXPECT_EQ(radials[8][1000], 0);
    ifs.close();
}

TEST_F(OrbIOTest, WriteAbacusOrb)
{
    std::ifstream ifs(file);
    std::string elem;
    double ecut, dr;
    int nr;
    std::vector<int> nzeta;
    std::vector<std::vector<double>> radials;
    ModuleIO::read_abacus_orb(ifs, elem, ecut, nr, dr, nzeta, radials);

    const std::string ftmp = "tmp.orb";
    std::ofstream ofs(ftmp);
    ModuleIO::write_abacus_orb(ofs, elem, ecut, nr, dr, nzeta, radials);
    ofs.close();

    std::ifstream ifs1(ftmp);
    std::string elem1;
    double ecut1, dr1;
    int nr1;
    std::vector<int> nzeta1;
    std::vector<std::vector<double>> radials1;
    ModuleIO::read_abacus_orb(ifs1, elem1, ecut1, nr1, dr1, nzeta1, radials1);
    EXPECT_EQ(elem, elem1);
    EXPECT_DOUBLE_EQ(ecut, ecut1);
    EXPECT_EQ(nr, nr1);
    EXPECT_DOUBLE_EQ(dr, dr1);
    EXPECT_EQ(nzeta.size(), nzeta1.size());
    for (int i = 0; i < nzeta.size(); ++i)
    {
        EXPECT_EQ(nzeta[i], nzeta1[i]);
    }
    EXPECT_EQ(radials.size(), radials1.size());
    for (int i = 0; i < radials.size(); ++i)
    {
        EXPECT_EQ(radials[i].size(), radials1[i].size());
        for (int j = 0; j < radials[i].size(); ++j)
        {
            EXPECT_NEAR(radials[i][j], radials1[i][j], tol);
        }
    }
    ifs1.close();
    remove(ftmp.c_str());
}

int main(int argc, char** argv)
{

#ifdef __MPI
    MPI_Init(&argc, &argv);
#endif

    testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();

#ifdef __MPI
    MPI_Finalize();
#endif

    return result;
}
