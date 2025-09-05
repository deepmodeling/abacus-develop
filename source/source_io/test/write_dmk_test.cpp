#include "source_io/write_dmk.h"

#define private public
#include "source_io/module_parameter/parameter.h"
#undef private
#include "source_base/global_variable.h"
#include "prepare_unitcell.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "source_base/module_external/scalapack_connector.h"

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
    this->tot_mag = 0.0;
    this->abs_mag = 0.0;
    this->start_mag = nullptr;
}
Magnetism::~Magnetism() { delete[] this->start_mag; }

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

template <typename T>
void gen_dmk(std::vector<std::vector<T>>& dmk, std::vector<double>& efs,  int nspin, int nk, int nlocal, Parallel_2D& pv)
{
    int myrank = 0;
#ifdef __MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
#endif
    std::vector<std::vector<T>> dmk_global(nspin * nk, std::vector<T>(nlocal * nlocal, T(0.0)));
    if (myrank == 0)
    {
        for (int i = 0; i < nspin * nk; i++)
        {
            for (int j = 0; j < nlocal * nlocal; j++)
            {
                std::complex<double> value = std::complex<double>(1.0 * i + 0.1 * j, 0.1 * i + 1.0 * j);
                dmk_global[i][j] = reinterpret_cast<T*>(&value)[0];
            }
        }
    }
#ifdef __MPI
    Parallel_2D pv_global;
    pv_global.init(nlocal, nlocal, nlocal, MPI_COMM_WORLD);
    dmk.resize(nspin * nk, std::vector<T>(pv.get_local_size(), 0.0));
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
    bool gamma_only = true;
    int ispin = 0;
    int nspin = 2;
    int ik = 0;
    int istep = 0;
    std::string fname = ModuleIO::dmk_gen_fname(gamma_only, ispin, nspin, ik, istep);
    EXPECT_EQ(fname, "dms1g1_nao.txt");

    ispin = 1;

    fname = ModuleIO::dmk_gen_fname(gamma_only, ispin, nspin, ik, istep);
    EXPECT_EQ(fname, "dms2g1_nao.txt");

    ispin = 0;
    gamma_only = false;    

    fname = ModuleIO::dmk_gen_fname(gamma_only, ispin, nspin, ik, istep);
    EXPECT_EQ(fname, "dmk1s1g1_nao.txt");

    ispin = 1;
    ik = 1;

    fname = ModuleIO::dmk_gen_fname(gamma_only, ispin, nspin, ik, istep);
    EXPECT_EQ(fname, "dmk2s2g1_nao.txt");
};


TEST(DMKTest,WriteDMK) {
    UnitCell* ucell;
    UcellTestPrepare utp = UcellTestLib["Si"];
    ucell = utp.SetUcellInfo();

    int nspin = 2;
    int nk = 1;
    int nk_multik = 2;
    int nlocal = 20;
    std::vector<std::vector<double>> dmk;
    std::vector<std::vector<std::complex<double>>> dmk_multik;
    Parallel_2D pv;
    std::vector<double> efs;
    init_pv(nlocal, pv);

    gen_dmk(dmk, efs, nspin, nk, nlocal, pv);
    gen_dmk(dmk_multik, efs, nspin, nk_multik, nlocal, pv);
    PARAM.sys.global_out_dir = "./";

    const int istep = -1;
    ModuleIO::write_dmk(dmk, 3, efs, ucell, pv, istep);
    ModuleIO::write_dmk(dmk_multik, 3, efs, ucell, pv, istep);
    std::ifstream ifs;

    int pass = 0;
    if (GlobalV::MY_RANK == 0)
    {
        std::string fn = "dms1_nao.txt";
        ifs.open(fn);
        std::string str((std::istreambuf_iterator<char>(ifs)),
                        std::istreambuf_iterator<char>());
        EXPECT_THAT(str, testing::HasSubstr("0.00000 (fermi energy)"));
        EXPECT_THAT(str, testing::HasSubstr("20 (number of basis)"));
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

        fn = "dms2_nao.txt";
        ifs.open(fn);
        str = std::string((std::istreambuf_iterator<char>(ifs)),
                          std::istreambuf_iterator<char>());
        EXPECT_THAT(str, testing::HasSubstr("0.10000 (fermi energy)"));
        EXPECT_THAT(str, testing::HasSubstr("20 (number of basis)"));
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

        fn = "dmk1s1_nao.txt";
        ifs.open(fn);
        str = std::string((std::istreambuf_iterator<char>(ifs)),
                          std::istreambuf_iterator<char>());
        EXPECT_THAT(str, testing::HasSubstr("0.00000 (fermi energy)"));
        EXPECT_THAT(str, testing::HasSubstr("20 (number of basis)"));
        EXPECT_THAT(
            str,
            testing::HasSubstr("(0.000e+00,0.000e+00) (1.000e-01,1.000e+00) (2.000e-01,2.000e+00) "
                               "(3.000e-01,3.000e+00)\n (4.000e-01,4.000e+00) (5.000e-01,5.000e+00) "
                               "(6.000e-01,6.000e+00) (7.000e-01,7.000e+00)\n"));
        EXPECT_THAT(
            str,
            testing::HasSubstr("(8.000e-01,8.000e+00) (9.000e-01,9.000e+00) (1.000e+00,1.000e+01) "
                               "(1.100e+00,1.100e+01)\n (1.200e+00,1.200e+01) (1.300e+00,1.300e+01) "
                               "(1.400e+00,1.400e+01) (1.500e+00,1.500e+01)\n"));
        EXPECT_THAT(
            str,
            testing::HasSubstr("(1.600e+00,1.600e+01) (1.700e+00,1.700e+01) (1.800e+00,1.800e+01) (1.900e+00,1.900e+01)\n"));
        ifs.close();

        fn = "dmk2s1_nao.txt";
        ifs.open(fn);
        str = std::string((std::istreambuf_iterator<char>(ifs)),
                          std::istreambuf_iterator<char>());
        EXPECT_THAT(str, testing::HasSubstr("0.00000 (fermi energy)"));
        EXPECT_THAT(str, testing::HasSubstr("20 (number of basis)"));
        EXPECT_THAT(
            str,
            testing::HasSubstr("(1.000e+00,1.000e-01) (1.100e+00,1.100e+00) (1.200e+00,2.100e+00) "
                               "(1.300e+00,3.100e+00)\n (1.400e+00,4.100e+00) (1.500e+00,5.100e+00) "
                               "(1.600e+00,6.100e+00) (1.700e+00,7.100e+00)\n"));
        EXPECT_THAT(
            str,
            testing::HasSubstr("(1.800e+00,8.100e+00) (1.900e+00,9.100e+00) (2.000e+00,1.010e+01) "
                               "(2.100e+00,1.110e+01)\n (2.200e+00,1.210e+01) (2.300e+00,1.310e+01) "
                               "(2.400e+00,1.410e+01) (2.500e+00,1.510e+01)\n"));
        EXPECT_THAT(
            str,
            testing::HasSubstr("(2.600e+00,1.610e+01) (2.700e+00,1.710e+01) (2.800e+00,1.810e+01) (2.900e+00,1.910e+01)\n"));
        ifs.close();

        fn = "dmk1s2_nao.txt";
        ifs.open(fn);
        str = std::string((std::istreambuf_iterator<char>(ifs)),
                          std::istreambuf_iterator<char>());
        EXPECT_THAT(str, testing::HasSubstr("0.10000 (fermi energy)"));
        EXPECT_THAT(str, testing::HasSubstr("20 (number of basis)"));
        EXPECT_THAT(
            str,
            testing::HasSubstr("(2.000e+00,2.000e-01) (2.100e+00,1.200e+00) (2.200e+00,2.200e+00) "
                               "(2.300e+00,3.200e+00)\n (2.400e+00,4.200e+00) (2.500e+00,5.200e+00) "
                               "(2.600e+00,6.200e+00) (2.700e+00,7.200e+00)\n"));
        EXPECT_THAT(
            str,
            testing::HasSubstr("(2.800e+00,8.200e+00) (2.900e+00,9.200e+00) (3.000e+00,1.020e+01) "
                               "(3.100e+00,1.120e+01)\n (3.200e+00,1.220e+01) (3.300e+00,1.320e+01) "
                               "(3.400e+00,1.420e+01) (3.500e+00,1.520e+01)\n"));
        EXPECT_THAT(
            str,
            testing::HasSubstr("(3.600e+00,1.620e+01) (3.700e+00,1.720e+01) (3.800e+00,1.820e+01) (3.900e+00,1.920e+01)\n"));
        ifs.close();

        fn = "dmk2s2_nao.txt";
        ifs.open(fn);
        str = std::string((std::istreambuf_iterator<char>(ifs)),
                          std::istreambuf_iterator<char>());
        EXPECT_THAT(str, testing::HasSubstr("0.10000 (fermi energy)"));
        EXPECT_THAT(str, testing::HasSubstr("20 (number of basis)"));
        EXPECT_THAT(
            str,
            testing::HasSubstr("(3.000e+00,3.000e-01) (3.100e+00,1.300e+00) (3.200e+00,2.300e+00) "
                               "(3.300e+00,3.300e+00)\n (3.400e+00,4.300e+00) (3.500e+00,5.300e+00) "
                               "(3.600e+00,6.300e+00) (3.700e+00,7.300e+00)\n"));
        EXPECT_THAT(
            str,
            testing::HasSubstr("(3.800e+00,8.300e+00) (3.900e+00,9.300e+00) (4.000e+00,1.030e+01) "
                               "(4.100e+00,1.130e+01)\n (4.200e+00,1.230e+01) (4.300e+00,1.330e+01) "
                               "(4.400e+00,1.430e+01) (4.500e+00,1.530e+01)\n"));
        EXPECT_THAT(
            str,
            testing::HasSubstr("(4.600e+00,1.630e+01) (4.700e+00,1.730e+01) (4.800e+00,1.830e+01) (4.900e+00,1.930e+01)\n"));
        ifs.close();
        remove("dms1_nao.txt");
        remove("dms2_nao.txt");
        remove("dmk1s1_nao.txt");
        remove("dmk2s1_nao.txt");
        remove("dmk1s2_nao.txt");
        remove("dmk2s2_nao.txt");
    }

    delete ucell;
    // remove the generated files
    
};


// no function in the main code calls read_dmk??? mohan note 2025-05-25
TEST(DMKTest, ReadDMK) {
    int nlocal = 26;
    std::vector<std::vector<double>> dmk;
    std::vector<std::vector<std::complex<double>>> dmk_multik;
    Parallel_2D pv;
    std::vector<double> efs;
    PARAM.sys.global_out_dir = "./";

    init_pv(nlocal, pv);

    std::ofstream ofs_running("running_log.txt");

    GlobalV::ofs_warning.open("warning.log");

    EXPECT_TRUE(ModuleIO::read_dmk(1, 1, pv, "./support/", dmk, ofs_running));
    ModuleIO::read_dmk(1, 1, pv, "./support/", dmk_multik, ofs_running);
    EXPECT_TRUE(ModuleIO::read_dmk(1, 1, pv, "./support/", dmk_multik, ofs_running));
    EXPECT_EQ(dmk.size(), 1);
    EXPECT_EQ(dmk_multik.size(), 1);
    EXPECT_EQ(dmk[0].size(), pv.get_local_size());
    EXPECT_EQ(dmk_multik[0].size(), pv.get_local_size());
    if (GlobalV::MY_RANK == 0)
    {
        EXPECT_NEAR(dmk[0][0], 3.904e-01, 1e-6);
        EXPECT_NEAR(dmk_multik[0][1].real(), -4.479e-03, 1e-6);
        EXPECT_NEAR(dmk_multik[0][1].imag(), 3.208e-04, 1e-6);
    }

    ofs_running.close();
    GlobalV::ofs_warning.close();
    remove("running_log.txt");
    remove("warning.log");
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
