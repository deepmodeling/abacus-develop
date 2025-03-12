#include "../write_wfc_nao.h"

#define private public
#include "module_parameter/parameter.h"
#undef private
#include "../binstream.h"
#include "module_base/global_variable.h"
#include "module_base/scalapack_connector.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#ifdef __MPI
#include "mpi.h"
#endif

TEST(GenWfcLcaoFnameTest, OutType1GammaOnlyOutAppFlagTrue)
{
    int out_type = 1;
    bool gamma_only = true;
    bool out_app_flag = true;
    int ik = 0;
    int istep = 0;

    std::string expected_output = "WFC_NAO_GAMMA1.txt";
    std::string result = ModuleIO::wfc_nao_gen_fname(out_type, gamma_only, out_app_flag, ik, istep);

    EXPECT_EQ(result, expected_output);
}

TEST(GenWfcLcaoFnameTest, OutType2GammaOnlyOutAppFlagFalse)
{
    int out_type = 2;
    bool gamma_only = true;
    bool out_app_flag = false;
    int ik = 1;
    int istep = 2;

    std::string expected_output = "WFC_NAO_GAMMA2_ION3.dat";
    std::string result = ModuleIO::wfc_nao_gen_fname(out_type, gamma_only, out_app_flag, ik, istep);

    EXPECT_EQ(result, expected_output);
}

TEST(GenWfcLcaoFnameTest, OutTypeInvalid)
{
    int out_type = 3;
    bool gamma_only = false;
    bool out_app_flag = true;
    int ik = 2;
    int istep = 3;

    std::string expected_output = "WFC_NAO_K3.txt";
    // catch the screen output
    testing::internal::CaptureStdout();
    std::string result = ModuleIO::wfc_nao_gen_fname(out_type, gamma_only, out_app_flag, ik, istep);
    std::string output = testing::internal::GetCapturedStdout();

    EXPECT_EQ(result, expected_output);
}

template <typename T>
void read_bin(const std::string& name, std::vector<T>& data)
{
    std::ifstream ifs(name, std::ios::binary);
    ifs.seekg(0, std::ios::beg);
    int nbands, nbasis;
    if (std::is_same<T, std::complex<double>>::value)
    {
        ifs.ignore(sizeof(int));
        ifs.ignore(sizeof(double) * 3);
    }
    ifs.read(reinterpret_cast<char*>(&nbands), sizeof(int));
    ifs.read(reinterpret_cast<char*>(&nbasis), sizeof(int));
    data.resize(nbands * nbasis);

    for (int i = 0; i < nbands; i++)
    {
        ifs.ignore(sizeof(int));
        ifs.ignore(sizeof(double) * 2);
        for (int j = 0; j < nbasis; j++)
        {
            ifs.read(reinterpret_cast<char*>(&data[i * nbasis + j]), sizeof(T));
        }
    }
    ifs.close();
}

class WriteWfcLcaoTest : public testing::Test
{
  protected:
    ModuleBase::matrix ekb;
    ModuleBase::matrix wg;
    std::vector<ModuleBase::Vector3<double>> kvec_c;
    std::vector<double> psi_init_double;
    std::vector<std::complex<double>> psi_init_complex;
    std::vector<double> psi_local_double;
    std::vector<std::complex<double>> psi_local_complex;
    Parallel_Orbitals pv;
    int nk = 2;
    int nbands = 3;
    int nbasis = 4;
    int nbands_local;
    int nbasis_local;

    void SetUp() override
    {
        PARAM.input.out_app_flag = true;
        ekb.create(nk, nbands); // in this test the value of ekb and wg is not important and not used.
        wg.create(nk, nbands);
        kvec_c.resize(nk, ModuleBase::Vector3<double>(0.0, 0.0, 0.0));
        psi_init_double.resize(nk * nbands * nbasis);
        psi_init_complex.resize(nk * nbands * nbasis);
        for (int i = 0; i < nk * nbands * nbasis; i++)
        {
            psi_init_double[i] = i * 0.1;
            psi_init_complex[i] = std::complex<double>(i * 0.2, i * 0.3);
        }
#ifdef __MPI
        pv.init(nbasis, nbands, 1, MPI_COMM_WORLD);
        for (int i = 0; i < 9; i++)
        {
            pv.desc_wfc[i] = pv.desc[i];
        }

        // create a global PV, and distribute the psi_init to local
        Parallel_2D pv_glb;
        pv_glb.init(nbasis, nbands, std::max(nbasis, nbands), MPI_COMM_WORLD);
        psi_local_double.resize(nk * pv.get_row_size() * pv.get_col_size());
        psi_local_complex.resize(nk * pv.get_row_size() * pv.get_col_size());
        for (int ik = 0; ik < nk; ik++)
        {
            Cpxgemr2d(nbasis,
                      nbands,
                      psi_init_double.data() + ik * nbands * nbasis,
                      1,
                      1,
                      pv_glb.desc,
                      psi_local_double.data() + ik * pv.get_row_size() * pv.get_col_size(),
                      1,
                      1,
                      pv.desc,
                      pv.blacs_ctxt);
            Cpxgemr2d(nbasis,
                      nbands,
                      psi_init_complex.data() + ik * nbands * nbasis,
                      1,
                      1,
                      pv_glb.desc,
                      psi_local_complex.data() + ik * pv.get_row_size() * pv.get_col_size(),
                      1,
                      1,
                      pv.desc,
                      pv.blacs_ctxt);
        }
        nbands_local = pv.get_col_size();
        nbasis_local = pv.get_row_size();
#else
        psi_local_double = psi_init_double;
        psi_local_complex = psi_init_complex;
        nbands_local = nbands;
        nbasis_local = nbasis;
#endif
    }
};

TEST_F(WriteWfcLcaoTest, WriteWfcLcao)
{
    // create a psi object
    psi::Psi<double> my_psi(psi_local_double.data(), nk, nbands_local, nbasis_local, nbasis_local, true);
    PARAM.sys.global_out_dir = "./";
    ModuleIO::write_wfc_nao(2, my_psi, ekb, wg, kvec_c, pv, -1);

    // check the output
    if (GlobalV::MY_RANK == 0)
    {
        for (int ik = 0; ik < nk; ik++)
        {
            std::string fname = ModuleIO::wfc_nao_gen_fname(2, true, PARAM.input.out_app_flag, ik, -1);
            std::ifstream file1(fname);
            EXPECT_TRUE(file1.good());
            std::vector<double> data;
            read_bin(fname, data);

            EXPECT_EQ(data.size(), nbands * nbasis);

            for (int i = 0; i < nbands * nbasis; i++)
            {
                EXPECT_DOUBLE_EQ(data[i], psi_init_double[nbands * nbasis * ik + i]);
            }
            // remove the output files
            std::remove(fname.c_str());
        }
    }
}

TEST_F(WriteWfcLcaoTest, WriteWfcLcaoComplex)
{
    psi::Psi<std::complex<double>> my_psi(psi_local_complex.data(), nk, nbands_local, nbasis_local, true);
    PARAM.sys.global_out_dir = "./";
    ModuleIO::write_wfc_nao(2, my_psi, ekb, wg, kvec_c, pv, -1);

    // check the output file
    if (GlobalV::MY_RANK == 0)
    {
        for (int ik = 0; ik < nk; ik++)
        {
            std::string fname = ModuleIO::wfc_nao_gen_fname(2, false, PARAM.input.out_app_flag, ik, -1);
            std::ifstream file1(fname);
            EXPECT_TRUE(file1.good());
            std::vector<std::complex<double>> data;
            read_bin(fname, data);

            EXPECT_EQ(data.size(), nbands * nbasis);
            for (int i = 0; i < nbands * nbasis; i++)
            {
                EXPECT_DOUBLE_EQ(data[i].real(), psi_init_complex[nbands * nbasis * ik + i].real());
                EXPECT_DOUBLE_EQ(data[i].imag(), psi_init_complex[nbands * nbasis * ik + i].imag());
            }
            // remove the output files
            std::remove(fname.c_str());
        }
    }
}

TEST(ModuleIOTest, WriteWfcNao)
{
    if (GlobalV::MY_RANK == 0)
    {
        // Set up GlobalV
        GlobalV::DRANK = 0;
        PARAM.input.nbands = 2;
        PARAM.sys.nlocal = 2;
        PARAM.input.out_app_flag = true;

        // Set up test data
        std::string filename = "test_wfc_nao.txt";
        std::vector<double> ctot = {0.1, 0.2, 0.3, 0.4};
        ModuleBase::matrix ekb(2, 2);
        ekb(0, 0) = 0.5;
        ekb(1, 0) = 0.6;
        ekb(0, 1) = 0.7;
        ekb(1, 1) = 0.8;
        ModuleBase::matrix wg(2, 2);
        wg(0, 0) = 0.9;
        wg(1, 0) = 1.0;
        wg(0, 1) = 1.1;
        wg(1, 1) = 1.2;

        // Call the function to be tested
        ModuleIO::wfc_nao_write2file(filename, ctot.data(), PARAM.sys.nlocal, 0, ekb, wg, false);

        // Check the output file
        std::ifstream ifs(filename);
        std::string str((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
        EXPECT_THAT(str, testing::HasSubstr("2 (number of bands)"));
        EXPECT_THAT(str, testing::HasSubstr("2 (number of orbitals)"));
        EXPECT_THAT(str, testing::HasSubstr("1 (band)"));
        EXPECT_THAT(str, testing::HasSubstr("5.00000000e-01 (Ry)"));
        EXPECT_THAT(str, testing::HasSubstr("9.00000000e-01 (Occupations)"));
        EXPECT_THAT(str, testing::HasSubstr("1.00000000e-01 2.00000000e-01"));
        EXPECT_THAT(str, testing::HasSubstr("2 (band)"));
        EXPECT_THAT(str, testing::HasSubstr("7.00000000e-01 (Ry)"));
        EXPECT_THAT(str, testing::HasSubstr("1.10000000e+00 (Occupations)"));
        EXPECT_THAT(str, testing::HasSubstr("3.00000000e-01 4.00000000e-01"));
        ifs.close();

        // clean up
        std::remove(filename.c_str()); // remove the test file
    }
}

TEST(ModuleIOTest, WriteWfcNaoBinary)
{
    if (GlobalV::MY_RANK == 0)
    {
        // Set up GlobalV
        GlobalV::DRANK = 0;
        PARAM.input.nbands = 2;
        PARAM.sys.nlocal = 2;
        PARAM.input.out_app_flag = true;

        // Set up test data
        std::string filename = "test_wfc_nao.dat";
        std::vector<double> ctot = {0.1, 0.2, 0.3, 0.4};
        ModuleBase::matrix ekb(2, 2);
        ekb(0, 0) = 0.5;
        ekb(1, 0) = 0.6;
        ekb(0, 1) = 0.7;
        ekb(1, 1) = 0.8;
        ModuleBase::matrix wg(2, 2);
        wg(0, 0) = 0.9;
        wg(1, 0) = 1.0;
        wg(0, 1) = 1.1;
        wg(1, 1) = 1.2;

        // Call the function to be tested
        ModuleIO::wfc_nao_write2file(filename, ctot.data(), PARAM.sys.nlocal, 0, ekb, wg, true);

        // Check the output file
        Binstream wfc(filename, "r");
        int nbands, nlocal;
        wfc >> nbands;
        wfc >> nlocal;
        EXPECT_EQ(nbands, 2);
        EXPECT_EQ(nlocal, 2);
        for (int i = 0; i < nbands; i++)
        {
            int band_index;
            double ekb, wg;
            wfc >> band_index;
            wfc >> ekb;
            wfc >> wg;
            EXPECT_EQ(band_index, i + 1);
            EXPECT_DOUBLE_EQ(ekb, 0.5 + i * 0.2);
            EXPECT_DOUBLE_EQ(wg, 0.9 + i * 0.2);
            for (int j = 0; j < nlocal; j++)
            {
                double ctot;
                wfc >> ctot;
                EXPECT_DOUBLE_EQ(ctot, 0.1 + i * 0.2 + j * 0.1);
            }
        }
        wfc.close();

        // clean up
        std::remove(filename.c_str()); // remove the test file
    }
}

TEST(ModuleIOTest, WriteWfcNaoComplex)
{
    if (GlobalV::MY_RANK == 0)
    {
        // Set up GlobalV
        PARAM.input.nbands = 2;
        PARAM.sys.nlocal = 3;
        PARAM.input.out_app_flag = true;
        // set up test data
        std::string name = "test_wfc_nao_complex.txt";
        int ik = 0;
        ModuleBase::Vector3<double> kvec_c{0.0, 0.0, 0.0};
        ModuleBase::matrix ekb(1, 2);
        ekb(0, 0) = 0.9;
        ekb(0, 1) = 1.1;
        ModuleBase::matrix wg(1, 2);
        wg(0, 0) = 0.11;
        wg(0, 1) = 0.22;
        std::vector<std::complex<double>> ctot = {std::complex<double>(1.0, 0.0),
                                                  std::complex<double>(2.0, 0.0),
                                                  std::complex<double>(3.0, 0.0),
                                                  std::complex<double>(0.0, 1.0),
                                                  std::complex<double>(0.0, 2.0),
                                                  std::complex<double>(0.0, 3.0)};

        // Call the function
        ModuleIO::wfc_nao_write2file_complex(name, ctot.data(), PARAM.sys.nlocal, ik, kvec_c, ekb, wg);
        // Check the output file
        std::ifstream ifs(name);
        std::string str((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
        EXPECT_THAT(str, testing::HasSubstr("1 (index of k points)"));
        EXPECT_THAT(str, testing::HasSubstr("0 0 0"));
        EXPECT_THAT(str, testing::HasSubstr("2 (number of bands)"));
        EXPECT_THAT(str, testing::HasSubstr("3 (number of orbitals)"));
        EXPECT_THAT(str, testing::HasSubstr("1 (band)"));
        EXPECT_THAT(str, testing::HasSubstr("(Ry)"));
        EXPECT_THAT(str, testing::HasSubstr("(Occupations)"));
        EXPECT_THAT(str, testing::HasSubstr("2 (band)"));
        // ifs.close();
        //  Clean up
        std::remove(name.c_str());
    }
}

TEST(ModuleIOTest, WriteWfcNaoComplexBinary)
{
    if (GlobalV::MY_RANK == 0)
    {
        // Set up GlobalV
        PARAM.input.nbands = 2;
        PARAM.sys.nlocal = 3;
        PARAM.input.out_app_flag = true;
        // set up test data
        std::string name = "test_wfc_nao_complex.dat";
        int ik = 0;
        ModuleBase::Vector3<double> kvec_c{0.0, 0.0, 0.0};
        ModuleBase::matrix ekb(1, 2);
        ekb(0, 0) = 0.9;
        ekb(0, 1) = 1.1;
        ModuleBase::matrix wg(1, 2);
        wg(0, 0) = 0.11;
        wg(0, 1) = 0.22;
        std::vector<std::complex<double>> ctot = {std::complex<double>(1.0, 0.0),
                                                  std::complex<double>(2.0, 2.0),
                                                  std::complex<double>(3.0, 4.0),
                                                  std::complex<double>(4.0, 4.0),
                                                  std::complex<double>(5.0, 6.0),
                                                  std::complex<double>(6.0, 8.0)};

        // Call the function
        ModuleIO::wfc_nao_write2file_complex(name, ctot.data(), PARAM.sys.nlocal, ik, kvec_c, ekb, wg, true);
        // Check the output file

        Binstream wfc(name, "r");
        int ik1, nbands, nlocal;
        double kx, ky, kz;
        wfc >> ik1;
        wfc >> kx;
        wfc >> ky;
        wfc >> kz;
        wfc >> nbands;
        wfc >> nlocal;
        EXPECT_EQ(ik1, 1);
        EXPECT_DOUBLE_EQ(kx, 0.0);
        EXPECT_DOUBLE_EQ(ky, 0.0);
        EXPECT_DOUBLE_EQ(kz, 0.0);
        EXPECT_EQ(nbands, 2);
        EXPECT_EQ(nlocal, 3);
        for (int i = 0; i < nbands; i++)
        {
            int band_index;
            double ekb, wg;
            wfc >> band_index;
            wfc >> ekb;
            wfc >> wg;
            EXPECT_EQ(band_index, i + 1);
            EXPECT_DOUBLE_EQ(ekb, 0.9 + i * 0.2);
            EXPECT_DOUBLE_EQ(wg, 0.11 + i * 0.11);
            for (int j = 0; j < nlocal; j++)
            {
                std::complex<double> ctot;
                wfc >> ctot;
                EXPECT_DOUBLE_EQ(ctot.real(), 1.0 + i * 3.0 + j * 1.0);
                EXPECT_DOUBLE_EQ(ctot.imag(), 0.0 + i * 4.0 + j * 2.0);
            }
        }
        wfc.close();
        // Clean up
        std::remove(name.c_str());
    }
}

int main(int argc, char** argv)
{
    GlobalV::MY_RANK = 0;
#ifdef __MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &GlobalV::NPROC);
    MPI_Comm_rank(MPI_COMM_WORLD, &GlobalV::MY_RANK);
#endif

    testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();
#ifdef __MPI
    MPI_Finalize();
#endif
    return result;
}
