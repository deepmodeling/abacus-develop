#include "gmock/gmock.h"
#include "gtest/gtest.h"

#define private public
#define protected public
#include "module_elecstate/module_charge/charge_extra.h"
#include "prepare_unitcell.h"

// mock functions for UnitCell
#ifdef __LCAO
InfoNonlocal::InfoNonlocal()
{
}
InfoNonlocal::~InfoNonlocal()
{
}
#endif
Magnetism::Magnetism()
{
    this->start_magnetization = nullptr;
}
Magnetism::~Magnetism()
{
    delete[] this->start_magnetization;
}

// mock functions for Charge
Charge::Charge()
{
    rhopw = new ModulePW::PW_Basis;
    rhopw->nrxx = 8;
    rhopw->nx = 2;
    rhopw->ny = 2;
    rhopw->nz = 2;
    rho = new double*[1];
    rho[0] = new double[rhopw->nrxx];
    ModuleBase::GlobalFunc::ZEROS(rho[0], rhopw->nrxx);
    for (int i = 0; i < rhopw->nrxx; ++i)
    {
        rho[0][i] = i + 1;
    }
}
Charge::~Charge()
{
    delete[] rho[0];
    delete[] rho;
    delete rhopw;
}
void Charge::atomic_rho(const int spin_number_need,
                        const double& omega,
                        double** rho_in,
                        const ModuleBase::ComplexMatrix& strucFac,
                        const UnitCell& ucell) const
{
}

// mock functions for PW_Basis
namespace ModulePW
{
PW_Basis::PW_Basis()
{
}
PW_Basis::~PW_Basis()
{
}
FFT::FFT()
{
}
FFT::~FFT()
{
}
void PW_Basis::initgrids(const double lat0_in, const ModuleBase::Matrix3 latvec_in, const double gridecut)
{
}
void PW_Basis::initgrids(const double lat0_in,
                         const ModuleBase::Matrix3 latvec_in,
                         const int nx_in,
                         int ny_in,
                         int nz_in)
{
}
void PW_Basis::distribute_r()
{
}
} // namespace ModulePW

// mock functions for Structure_Factor
Structure_Factor::Structure_Factor()
{
}
Structure_Factor::~Structure_Factor()
{
}
void Structure_Factor::setup_structure_factor(UnitCell* Ucell, const ModulePW::PW_Basis* rho_basis)
{
}

/************************************************
 *  unit test of module_charge/charge_extra.cpp
 ***********************************************/

/**
 * - Tested Functions:
 *   - Charge_Extra::Init_CE()
 *     - Initialization of viriables used in charge extrapolation methods
 *   - Charge_Extra::extrapolate_charge()
 *     - charge extrapolation
 *   - Charge_Extra::update_all_dis()
 *     - update displacements
 *   - Charge_Extra::save_files()
 *     - save the difference of the convergent charge density and the initial atomic charge density
 *   - Charge_Extra::find_alpha_and_beta()
 *     - determine alpha and beta
 *   - Charge_Extra::read_files()
 *     - read cube files containing the charge difference of previous steps
 */

class ChargeExtraTest : public ::testing::Test
{
  protected:
    Charge_Extra CE;
    UcellTestPrepare utp = UcellTestLib["Si"];
    std::unique_ptr<UnitCell> ucell;
    Charge charge;
    Structure_Factor sf;
    void SetUp() override
    {
        GlobalV::NSPIN = 1;
        GlobalV::global_out_dir = "./support/";
        ucell = utp.SetUcellInfo();
        ucell->omega = 1.0;
    }
    void TearDown() override
    {
    }
};

TEST_F(ChargeExtraTest, InitCECase1)
{
    GlobalV::chg_extrap = "none";
    CE.Init_CE(ucell->nat);
    EXPECT_EQ(CE.pot_order, 0);
}

TEST_F(ChargeExtraTest, InitCECase2)
{
    GlobalV::chg_extrap = "atomic";
    CE.Init_CE(ucell->nat);
    EXPECT_EQ(CE.pot_order, 1);
}

TEST_F(ChargeExtraTest, InitCECase3)
{
    GlobalV::chg_extrap = "first-order";
    CE.Init_CE(ucell->nat);
    EXPECT_EQ(CE.pot_order, 2);
}

TEST_F(ChargeExtraTest, InitCECase4)
{
    GlobalV::chg_extrap = "second-order";
    CE.Init_CE(ucell->nat);
    EXPECT_EQ(CE.pot_order, 3);
    EXPECT_DOUBLE_EQ(CE.alpha, 1.0);
    EXPECT_DOUBLE_EQ(CE.beta, 0.0);
    EXPECT_NE(CE.dis_old1, nullptr);
    EXPECT_NE(CE.dis_old2, nullptr);
    EXPECT_NE(CE.dis_now, nullptr);
}

TEST_F(ChargeExtraTest, InitCEWarningQuit)
{
    GlobalV::chg_extrap = "error";
    testing::internal::CaptureStdout();
    EXPECT_EXIT(CE.Init_CE(ucell->nat), ::testing::ExitedWithCode(0), "");
    std::string output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("charge extrapolation method is not available"));
}

TEST_F(ChargeExtraTest, SaveFiles)
{
    int istep = 2;
    CE.pot_order = 3;

    CE.save_files(istep, *ucell.get(), &charge, &sf);

    std::ifstream input1("./support/NOW_SPIN1_CHG.cube");
    EXPECT_TRUE(input1.good());
    input1.close();

    std::ifstream input2("./support/OLD1_SPIN1_CHG.cube");
    EXPECT_TRUE(input2.good());
    input2.close();

    std::ifstream input3("./support/OLD2_SPIN1_CHG.cube");
    EXPECT_TRUE(input3.good());
    input3.close();
}

TEST_F(ChargeExtraTest, ReadFiles)
{
    CE.read_files(2, 2, 2, ucell.get(), "OLD1", charge.rho);
    EXPECT_DOUBLE_EQ(charge.rho[0][0], 1.0);
    EXPECT_DOUBLE_EQ(charge.rho[0][4], 5.0);
}

TEST_F(ChargeExtraTest, ExtrapolateChargeCase1)
{
    CE.istep = 0;
    CE.pot_order = 3;

    GlobalV::ofs_running.open("log");
    CE.extrapolate_charge(*ucell.get(), &charge, &sf);
    GlobalV::ofs_running.close();

    // Check the results
    std::ifstream ifs("log");
    std::string expected_output = " charge density from previous step !\n";
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    ifs.close();
    std::remove("log");

    EXPECT_EQ(output, expected_output);
    EXPECT_EQ(CE.rho_extr, 0);
}

TEST_F(ChargeExtraTest, ExtrapolateChargeCase2)
{
    CE.istep = 1;
    CE.pot_order = 3;

    GlobalV::ofs_running.open("log");
    CE.extrapolate_charge(*ucell.get(), &charge, &sf);
    GlobalV::ofs_running.close();

    // Check the results
    std::ifstream ifs("log");
    std::string expected_output = " Find the file, try to read charge from file.\n read in fermi energy = 0\n Read "
                                  "SPIN = 1 charge now.\n NEW-OLD atomic charge density approx. for the potential !\n";
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    ifs.close();
    std::remove("log");

    EXPECT_EQ(output, expected_output);
    EXPECT_EQ(CE.rho_extr, 1);
}

TEST_F(ChargeExtraTest, ExtrapolateChargeCase3)
{
    CE.istep = 2;
    CE.pot_order = 3;

    GlobalV::ofs_running.open("log");
    CE.extrapolate_charge(*ucell.get(), &charge, &sf);
    GlobalV::ofs_running.close();

    // Check the results
    std::ifstream ifs("log");
    std::string expected_output
        = " Find the file, try to read charge from file.\n read in fermi energy = 0\n Read SPIN = 1 charge now.\n "
          "first order charge density extrapolation !\n Find the file, try to read charge from file.\n read in fermi "
          "energy = 0\n Read SPIN = 1 charge now.\n";
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    ifs.close();
    std::remove("log");

    EXPECT_EQ(output, expected_output);
    EXPECT_EQ(CE.rho_extr, 2);
}

TEST_F(ChargeExtraTest, ExtrapolateChargeCase4)
{
    GlobalV::chg_extrap = "second-order";
    CE.Init_CE(ucell->nat);
    CE.istep = 3;

    GlobalV::ofs_running.open("log");
    CE.extrapolate_charge(*ucell.get(), &charge, &sf);
    GlobalV::ofs_running.close();

    // Check the results
    std::ifstream ifs("log");
    std::string expected_output
        = " Find the file, try to read charge from file.\n read in fermi energy = 0\n Read SPIN = 1 charge now.\n "
          "second order charge density extrapolation !\n alpha = 0\n beta = 0\n Find the file, try to read charge from "
          "file.\n read in fermi energy = 0\n Read SPIN = 1 charge now.\n Find the file, try to read charge from "
          "file.\n read in fermi energy = 0\n Read SPIN = 1 charge now.\n";
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    ifs.close();
    std::remove("log");

    EXPECT_EQ(output, expected_output);
    EXPECT_EQ(CE.rho_extr, 3);
    std::remove("./support/OLD2_SPIN1_CHG.cube");
}

TEST_F(ChargeExtraTest, UpdateAllDis)
{
    GlobalV::chg_extrap = "second-order";
    CE.Init_CE(ucell->nat);
    CE.istep = 3;
    for (int i = 0; i < ucell->nat; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            CE.dis_old1[i][j] = i;
            CE.dis_now[i][j] = j;
        }
    }

    CE.update_all_dis(*ucell.get());

    EXPECT_EQ(CE.istep, 4);
    EXPECT_DOUBLE_EQ(CE.dis_old2[0][2], 0.0);
    EXPECT_DOUBLE_EQ(CE.dis_old1[0][2], 2.0);
    EXPECT_DOUBLE_EQ(CE.dis_now[0][2], 0.0);
}

TEST_F(ChargeExtraTest, FindAlphaAndBeta)
{
    GlobalV::chg_extrap = "second-order";
    CE.Init_CE(ucell->nat);
    CE.istep = 3;
    for (int i = 0; i < ucell->nat; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            CE.dis_old1[i][j] = i;
            CE.dis_now[i][j] = j;
        }
    }

    CE.find_alpha_and_beta(ucell->nat);

    EXPECT_DOUBLE_EQ(CE.alpha, 1.0);
    EXPECT_DOUBLE_EQ(CE.beta, 0.0);
}
