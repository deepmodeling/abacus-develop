#include "module_io/io_dmk.h"

#define private public
#include "module_parameter/parameter.h"
#undef private
#include "module_base/global_variable.h"
#include "prepare_unitcell.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "module_base/scalapack_connector.h"

#ifdef __MPI
#include "mpi.h"
#endif

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

void init_pv(int nlocal, Parallel_2D& pv)
{
#ifdef __MPI
        pv.init(nlocal, nlocal, 1, MPI_COMM_WORLD);
#else
        pv.nrow = nlocal;
        pv.ncol = nlocal;  
#endif             
}

void gen_dmk(std::vector<std::vector<double>>& dmk, std::vector<double>& efs,  int nspin, int nk, int nlocal, Parallel_2D& pv)
{
    int myrank = 0;
#ifdef __MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
#endif
    std::vector<std::vector<double>> dmk_global(nspin * nk, std::vector<double>(nlocal * nlocal, 0.0));
    if (myrank == 0)
    {
        for (int i = 0; i < nspin * nk; i++)
        {
            for (int j = 0; j < nlocal * nlocal; j++)
            {
                dmk_global[i][j] = 1.0 * i + 0.1 * j;
            }
        }
    }
#ifdef __MPI
    Parallel_2D pv_global;
    pv_global.init(nlocal, nlocal, nlocal, MPI_COMM_WORLD);
    dmk.resize(nspin * nk, std::vector<double>(pv.get_local_size(), 0.0));
    for (int i = 0; i < nspin * nk; i++)
    {
        Cpxgemr2d(nlocal,
                  nlocal,
                  dmk_global[i].data(),
                  1,
                  1,
                  pv_global.desc,
                  dmk[i].data(),
                  1,
                  1,
                  pv.desc,
                  pv.blacs_ctxt);
    }
#else
    dmk = dmk_global;
#endif

    efs.resize(nspin, 0.0);
    for (int i = 0; i < nspin; i++)
    {
        efs[i] = 0.1 * i;
    }
}


TEST(DMKTest, GenFileName) {
    std::string fname = ModuleIO::dmk_gen_fname(true, 0, 0);
    EXPECT_EQ(fname, "SPIN1_DM");
    fname = ModuleIO::dmk_gen_fname(true, 1, 1);
    EXPECT_EQ(fname, "SPIN2_DM");

    std::string output;
    testing::internal::CaptureStdout();
    EXPECT_EXIT(ModuleIO::dmk_gen_fname(false, 2, 0),
                ::testing::ExitedWithCode(1),
                "");
    output = testing::internal::GetCapturedStdout();
};


TEST(DMKTest,WriteDMK) {
    UnitCell* ucell;
    UcellTestPrepare utp = UcellTestLib["Si"];
    ucell = utp.SetUcellInfo();

    int nspin = 2;
    int nk = 1;
    int nlocal = 20;
    std::vector<std::vector<double>> dmk;
    Parallel_2D pv;
    std::vector<double> efs;
    init_pv(nlocal, pv);

    gen_dmk(dmk, efs, nspin, nk, nlocal, pv);
    PARAM.sys.global_out_dir = "./";

    ModuleIO::write_dmk(dmk, 3, efs, ucell, pv);
    std::ifstream ifs;

    int pass = 0;
    if (GlobalV::MY_RANK == 0)
    {
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
        remove("SPIN1_DM");
        remove("SPIN2_DM");
    }

    delete ucell;
    // remove the generated files
    
};



TEST(DMKTest, ReadDMK) {
    int nlocal = 26;
    std::vector<std::vector<double>> dmk;
    Parallel_2D pv;
    std::vector<double> efs;
    PARAM.sys.global_out_dir = "./";

    init_pv(nlocal, pv);
    EXPECT_TRUE(ModuleIO::read_dmk(1, 1, pv, "./support/", dmk));
    EXPECT_EQ(dmk.size(), 1);
    EXPECT_EQ(dmk[0].size(), pv.get_local_size());
    if (GlobalV::MY_RANK == 0)
    {
        EXPECT_NEAR(dmk[0][0], 3.904e-01, 1e-6);
    }
}


#ifdef __MPI
int main(int argc, char** argv)
{
    GlobalV::MY_RANK = 0;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &GlobalV::NPROC);
    MPI_Comm_rank(MPI_COMM_WORLD, &GlobalV::MY_RANK);

    testing::InitGoogleTest(&argc, argv);
    ::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();

    if (GlobalV::MY_RANK != 0) {
        delete listeners.Release(listeners.default_result_printer());
    }

    int result = RUN_ALL_TESTS();
    MPI_Bcast(&result, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (GlobalV::MY_RANK == 0 && result != 0)
    {
        std::cout << "ERROR:some tests are not passed" << std::endl;
	}

    MPI_Finalize();
    return result;
}
#endif