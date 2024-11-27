#include "module_io/io_dmk.h"

#define private public
#include "module_parameter/parameter.h"
#undef private
#include "module_base/global_variable.h"
#include "prepare_unitcell.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#ifdef __LCAO
InfoNonlocal::InfoNonlocal() {}
InfoNonlocal::~InfoNonlocal() {}
LCAO_Orbitals::LCAO_Orbitals() {}
LCAO_Orbitals::~LCAO_Orbitals() {}
#endif
Magnetism::Magnetism() {
    this->tot_magnetization = 0.0;
    this->abs_magnetization = 0.0;
    this->start_magnetization = nullptr;
}
Magnetism::~Magnetism() { delete[] this->start_magnetization; }

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

TEST(DMKTest, GenFileName) {
    std::string fname = ModuleIO::dmk_gen_fname(true, 0, 0);
    EXPECT_EQ(fname, "SPIN1_DM");
    fname = ModuleIO::dmk_gen_fname(true, 1, 1);
    EXPECT_EQ(fname, "SPIN2_DM");

    // do not support non-gamma-only case now
    std::string output;
    testing::internal::CaptureStdout();
    EXPECT_EXIT(ModuleIO::dmk_gen_fname(false, 2, 0),
                ::testing::ExitedWithCode(1),
                "");
    output = testing::internal::GetCapturedStdout();
};

class DMKIOTest : public ::testing::Test {
  protected:
    int nspin = 2;
    int nk = 1;
    int nlocal = 20;
    std::vector<std::vector<double>> dmk;
    Parallel_2D pv;
    std::vector<double> efs;

    void SetUp() {
        dmk.resize(nspin * nk, std::vector<double>(nlocal * nlocal, 0.0));
        for (int i = 0; i < nspin * nk; i++) {
            for (int j = 0; j < nlocal * nlocal; j++) {
                dmk[i][j] = 1.0 * i + 0.1 * j;
            }
        }

        efs.resize(nspin, 0.0);
        for (int i = 0; i < nspin; i++) {
            efs[i] = 0.1 * i;
        }

        pv.nrow = nlocal;
        pv.ncol = nlocal;

        PARAM.sys.global_out_dir = "./";
    }
};

TEST_F(DMKIOTest, WriteDMK) {
    UnitCell* ucell;
    UcellTestPrepare utp = UcellTestLib["Si"];
    ucell = utp.SetUcellInfo();
    ModuleIO::write_dmk(dmk, 3, efs, ucell, pv);
    std::ifstream ifs;

    std::string fn = "SPIN1_DM";
    ifs.open(fn);
    std::string str((std::istreambuf_iterator<char>(ifs)),
                    std::istreambuf_iterator<char>());
    EXPECT_THAT(str, testing::HasSubstr("0.00000 (fermi energy)"));
    EXPECT_THAT(str, testing::HasSubstr("20 20"));
    EXPECT_THAT(
        str,
        testing::HasSubstr("0.000e+00 1.000e-01 2.000e-01 3.000e-01 4.000e-01 "
                           "5.000e-01 6.000e-01 7.000e-01\n"));
    EXPECT_THAT(
        str,
        testing::HasSubstr("8.000e-01 9.000e-01 1.000e+00 1.100e+00 1.200e+00 "
                           "1.300e+00 1.400e+00 1.500e+00\n"));
    EXPECT_THAT(
        str,
        testing::HasSubstr("1.600e+00 1.700e+00 1.800e+00 1.900e+00\n"));
    ifs.close();

    fn = "SPIN2_DM";
    ifs.open(fn);
    str = std::string((std::istreambuf_iterator<char>(ifs)),
                      std::istreambuf_iterator<char>());
    EXPECT_THAT(str, testing::HasSubstr("0.10000 (fermi energy)"));
    EXPECT_THAT(str, testing::HasSubstr("20 20"));
    EXPECT_THAT(
        str,
        testing::HasSubstr("1.000e+00 1.100e+00 1.200e+00 1.300e+00 1.400e+00 "
                           "1.500e+00 1.600e+00 1.700e+00\n"));
    EXPECT_THAT(
        str,
        testing::HasSubstr("1.800e+00 1.900e+00 2.000e+00 2.100e+00 2.200e+00 "
                           "2.300e+00 2.400e+00 2.500e+00\n"));
    EXPECT_THAT(
        str,
        testing::HasSubstr("2.600e+00 2.700e+00 2.800e+00 2.900e+00\n"));
    ifs.close();

    delete ucell;
    // remove the generated files
    remove("SPIN1_DM");
    remove("SPIN2_DM");
};

TEST_F(DMKIOTest, ReadDMK) {
    pv.nrow = 26;
    pv.ncol = 26;
    EXPECT_TRUE(ModuleIO::read_dmk(1, 1, pv, "./support/", dmk));
    EXPECT_EQ(dmk.size(), 1);
    EXPECT_EQ(dmk[0].size(), 26 * 26);
    EXPECT_NEAR(dmk[0][0], 3.904e-01, 1e-6);
    EXPECT_NEAR(dmk[0][25 * 26 + 25], 3.445e-02, 1e-6);
}
