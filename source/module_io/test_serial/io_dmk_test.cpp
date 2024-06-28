#include "module_io/io_dmk.h"

#include "module_base/global_variable.h"
#include "prepare_unitcell.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

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
 *  unit test of read_dmk and write_dmk
 ***********************************************/

/**
 * - Tested Functions:
 *   - read_dmk()
 *     - the function to read density matrix K from file
 *     - the serial version without MPI
 *   - write_dmk()
 *     - the function to write density matrix K to file
 *     - the serial version without MPI
 */

class DMIOTest : public ::testing::Test
{
  protected:
    int nspin = 1;
    int lgd = 26;
    int nnrg = 26 * 26;
    double*** DM;
    double** DM_R;
    UnitCell* ucell;
    void SetUp()
    {
        DM = new double**[nspin];
        DM_R = new double*[nspin];
        ucell = new UnitCell;
        for (int is = 0; is < nspin; ++is)
        {
            DM[is] = new double*[lgd];
            DM_R[is] = new double[nnrg];
            for (int ig = 0; ig < lgd; ++ig)
            {
                DM[is][ig] = new double[lgd];
            }
        }
    }
    void TearDown()
    {
        for (int is = 0; is < nspin; ++is)
        {
            for (int ig = 0; ig < lgd; ++ig)
            {
                delete[] DM[is][ig];
            }
            delete[] DM[is];
            delete[] DM_R[is];
        }
        delete[] DM;
        delete[] DM_R;
        delete ucell;
    }
};

TEST_F(DMIOTest, Read)
{
    GlobalV::MY_RANK = 0;
    GlobalV::NLOCAL = lgd;
    GlobalV::GAMMA_ONLY_LOCAL = true;
    int is = 0;
    std::string fn = "./support/SPIN1_DM";
    double ef;
    UcellTestPrepare utp = UcellTestLib["Si"];
    ucell = utp.SetUcellInfo();
    ModuleIO::read_dmk(GlobalV::GAMMA_ONLY_LOCAL, GlobalV::NLOCAL, GlobalV::NSPIN, is, fn, DM, DM_R, ef, ucell);
    EXPECT_DOUBLE_EQ(ef, 0.570336288802337);
    EXPECT_NEAR(DM[0][0][0], 3.904e-01, 1e-6);
    EXPECT_NEAR(DM[0][25][25], 3.445e-02, 1e-6);
}

TEST_F(DMIOTest, Write)
{
    // first read
    GlobalV::MY_RANK = 0;
    GlobalV::NLOCAL = lgd;
    GlobalV::GAMMA_ONLY_LOCAL = true;
    int is = 0;
    std::string fn = "./support/SPIN1_DM";
    double ef;
    UcellTestPrepare utp = UcellTestLib["Si"];
    ucell = utp.SetUcellInfo();
    ModuleIO::read_dmk(GlobalV::GAMMA_ONLY_LOCAL, GlobalV::NLOCAL, GlobalV::NSPIN, is, fn, DM, DM_R, ef, ucell);
    EXPECT_DOUBLE_EQ(ef, 0.570336288802337);
    EXPECT_NEAR(DM[0][0][0], 3.904e-01, 1e-6);
    EXPECT_NEAR(DM[0][25][25], 3.445e-02, 1e-6);
    // then write
    int precision = 3;
    std::vector<std::vector<double>> dmk(nspin, std::vector<double>(lgd * lgd, 0.0));
    for (int is = 0; is < nspin; ++is)
    {
        for (int ig = 0; ig < lgd; ++ig)
        {
            for (int jg = 0; jg < lgd; ++jg)
            {
                dmk[is][ig * lgd + jg] = DM[is][ig][jg];
            }
        }
    }
    Parallel_2D pv;
    pv.nrow = lgd;
    pv.ncol = lgd;
    ModuleIO::write_dmk(dmk, precision, std::vector<double>(nspin, ef), ucell, pv);
    std::ifstream ifs;
    ifs.open("SPIN1_DM");
    std::string str((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_THAT(str, testing::HasSubstr("0.570336288802337 (fermi energy)"));
    ifs.close();
    // remove("SPIN1_DM");
}
