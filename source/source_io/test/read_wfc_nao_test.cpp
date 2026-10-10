#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "source_io/module_wf/read_wfc_nao.h"
#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_io/module_wf/write_wfc_nao.h"
#include "source_base/module_out/filename.h"
#include "source_base/global_function.h"

#include <cstdio>

/************************************************
 *  unit test of functions in read_wfc_nao.cpp
 ***********************************************/

/**
 * - Tested Functions:
 *   - distri_wfc_nao()
 *     - calculate memory required.
 *   - read_wfc_nao()
 *     - read wave functions from file.
 */

class ReadWfcNaoTest : public ::testing::Test
{
protected:
    int my_rank = 0;
    int nproc = 1;
    std::string binary_test_dir;

    void SetUp() override
    {
#ifdef __MPI
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &nproc);
#endif
        binary_test_dir = "./read_wfc_nao_binary_np" + std::to_string(nproc) + "/";
        ModuleBase::GlobalFunc::MAKE_DIR(binary_test_dir);
#ifdef __MPI
        MPI_Barrier(MPI_COMM_WORLD);
#endif
    }

    void TearDown() override
    {
#ifdef __MPI
        MPI_Barrier(MPI_COMM_WORLD);
#endif
        if (my_rank == 0)
        {
            std::remove((binary_test_dir + "wf_nao.dat").c_str());
            std::remove((binary_test_dir + "wf_nao.txt").c_str());
            std::remove((binary_test_dir + "wfk1_nao.dat").c_str());
            std::remove(binary_test_dir.substr(0, binary_test_dir.size() - 1).c_str());
        }
#ifdef __MPI
        MPI_Barrier(MPI_COMM_WORLD);
#endif
    }

    void initialize_parallel_orbitals(Parallel_Orbitals& para, const int nlocal, const int nbands)
    {
#ifdef __MPI
        std::ofstream ofs_running, ofs_warning;
        para.init(nlocal, nlocal, 1, MPI_COMM_WORLD);
        para.set_nloc_wfc_Eij(nbands, ofs_running, ofs_warning);
        para.set_desc_wfc_Eij(nlocal, nbands, para.nrow);
#else
        para.set_serial(nlocal, nlocal);
        para.nrow_bands = nlocal;
        para.ncol_bands = nbands;
#endif
    }
};


TEST_F(ReadWfcNaoTest,ReadWfcNao)
{
      //Global variables
      const int nbands = 3;
      const int nlocal = 3;
      const std::string readin_dir = "./support/";
      const int nks = 1;
      const int nspin = 1;
      int my_rank = 0;

      Parallel_Orbitals ParaV;
#ifdef __MPI
      MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
      std::ofstream ofs_running, ofs_warning;
      ParaV.init(nlocal, nlocal, 1, MPI_COMM_WORLD);   
      ParaV.set_nloc_wfc_Eij(nbands, ofs_running, ofs_warning);   
      ParaV.set_desc_wfc_Eij(nlocal, nbands, ParaV.nrow);
#else
      ParaV.set_serial(nlocal, nlocal);
      ParaV.nrow_bands = nlocal;
      ParaV.ncol_bands = nbands;  
#endif 

      psi::Psi<double> psid;
      ModuleBase::matrix ekb;
      ModuleBase::matrix wg;
      ekb.create(nks,nbands);
      wg.create(nks,nbands);

      std::vector<int> ik2iktot = {0};
      const int nkstot = 1;

      // Act
	  ModuleIO::read_wfc_nao(readin_dir, ParaV, psid,
			  ekb, wg, ik2iktot, nkstot, nspin, false);
      // Assert
      EXPECT_NEAR(ekb(0,1),0.31482195194888534794941393,1e-5);
      EXPECT_NEAR(wg(0,1),0.0,1e-5);
      if (my_rank == 0)
      {
            EXPECT_NEAR(psid(0,0,0),5.3759239842e-01,1e-5);
      }
}

TEST_F(ReadWfcNaoTest, ReadWfcNaoPart)
{
    //Global variables
    const int nbands = 2;
    const int skip_band = 1;
    const int nlocal = 3;
    const std::string readin_dir = "./support/";
    const int nks = 1;
    const int nspin = 1;
    const int nstep = -1;
    int my_rank = 0;

    Parallel_Orbitals ParaV;
#ifdef __MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    std::ofstream ofs_running, ofs_warning;
    ParaV.init(nlocal, nlocal, 1, MPI_COMM_WORLD);
    ParaV.set_nloc_wfc_Eij(nbands, ofs_running, ofs_warning); 
    ParaV.set_desc_wfc_Eij(nlocal, nbands, ParaV.nrow);
#else
    ParaV.set_serial(nlocal, nlocal);
    ParaV.nrow_bands = nlocal;
    ParaV.ncol_bands = nbands;
#endif 

    psi::Psi<double> psid;
    ModuleBase::matrix ekb;
    ModuleBase::matrix wg;
    ekb.create(nks, nbands);
    wg.create(nks, nbands);

	std::vector<int> ik2iktot = {0};
	const int nkstot = 1;

	// Act
	ModuleIO::read_wfc_nao(readin_dir, ParaV, psid,
			ekb, wg, ik2iktot, nkstot, nspin, false, skip_band, nstep);

    // Assert
    EXPECT_NEAR(ekb(0, 1), 7.4141254894954844445464914e-01, 1e-5);
    if (my_rank == 0)
    {
        EXPECT_NEAR(psid(0, 0, 0), 1.8587183851, 1e-5);
    }
}

TEST_F(ReadWfcNaoTest, ReadBinaryGamma)
{
    const int nbands = 2;
    const int nlocal = 3;
    Parallel_Orbitals para;
    initialize_parallel_orbitals(para, nlocal, nbands);

    const std::vector<double> coefficients = {0.1, 0.2, 0.3, 1.1, 1.2, 1.3};
    ModuleBase::matrix ekb_source(1, nbands);
    ModuleBase::matrix wg_source(1, nbands);
    ekb_source(0, 0) = -0.5;
    ekb_source(0, 1) = 0.7;
    wg_source(0, 0) = 2.0;
    wg_source(0, 1) = 0.0;
    if (my_rank == 0)
    {
        ModuleIO::wfc_nao_write2file(binary_test_dir + "wf_nao.dat",
                                     coefficients.data(),
                                     nlocal,
                                     0,
                                     ekb_source,
                                     wg_source,
                                     true,
                                     false);
    }
#ifdef __MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    psi::Psi<double> psi_read;
    ModuleBase::matrix ekb(1, nbands);
    ModuleBase::matrix wg(1, nbands);
    const std::vector<int> ik2iktot = {0};
    EXPECT_TRUE(ModuleIO::read_wfc_nao(binary_test_dir,
                                       para,
                                       psi_read,
                                       ekb,
                                       wg,
                                       ik2iktot,
                                       1,
                                       1,
                                       true));
    EXPECT_DOUBLE_EQ(ekb(0, 0), ekb_source(0, 0));
    EXPECT_DOUBLE_EQ(ekb(0, 1), ekb_source(0, 1));
    EXPECT_DOUBLE_EQ(wg(0, 0), wg_source(0, 0));
    EXPECT_DOUBLE_EQ(wg(0, 1), wg_source(0, 1));
    if (my_rank == 0)
    {
        EXPECT_DOUBLE_EQ(psi_read(0, 0, 0), coefficients[0]);
    }
}

TEST_F(ReadWfcNaoTest, ReadBinaryGammaFloat)
{
    const int nbands = 2;
    const int nlocal = 3;
    Parallel_Orbitals para;
    initialize_parallel_orbitals(para, nlocal, nbands);

    const std::vector<double> coefficients = {0.1, 0.2, 0.3, 1.1, 1.2, 1.3};
    ModuleBase::matrix ekb_source(1, nbands);
    ModuleBase::matrix wg_source(1, nbands);
    ekb_source(0, 0) = -0.5;
    ekb_source(0, 1) = 0.7;
    wg_source(0, 0) = 2.0;
    wg_source(0, 1) = 0.0;
    if (my_rank == 0)
    {
        ModuleIO::wfc_nao_write2file(binary_test_dir + "wf_nao.dat",
                                     coefficients.data(),
                                     nlocal,
                                     0,
                                     ekb_source,
                                     wg_source,
                                     true,
                                     false);
    }
#ifdef __MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    psi::Psi<float> psi_read;
    ModuleBase::matrix ekb(1, nbands);
    ModuleBase::matrix wg(1, nbands);
    const std::vector<int> ik2iktot = {0};
    EXPECT_TRUE(ModuleIO::read_wfc_nao(binary_test_dir,
                                       para,
                                       psi_read,
                                       ekb,
                                       wg,
                                       ik2iktot,
                                       1,
                                       1,
                                       true));
    EXPECT_DOUBLE_EQ(ekb(0, 0), ekb_source(0, 0));
    EXPECT_DOUBLE_EQ(ekb(0, 1), ekb_source(0, 1));
    EXPECT_DOUBLE_EQ(wg(0, 0), wg_source(0, 0));
    EXPECT_DOUBLE_EQ(wg(0, 1), wg_source(0, 1));
    if (my_rank == 0)
    {
        EXPECT_FLOAT_EQ(psi_read(0, 0, 0), static_cast<float>(coefficients[0]));
    }
}

TEST_F(ReadWfcNaoTest, ReadBinaryComplex)
{
    const int nbands = 2;
    const int nlocal = 3;
    Parallel_Orbitals para;
    initialize_parallel_orbitals(para, nlocal, nbands);

    const std::vector<std::complex<double>> coefficients
        = {{0.1, -0.1}, {0.2, -0.2}, {0.3, -0.3}, {1.1, 0.4}, {1.2, 0.5}, {1.3, 0.6}};
    ModuleBase::matrix ekb_source(1, nbands);
    ModuleBase::matrix wg_source(1, nbands);
    ekb_source(0, 0) = -0.4;
    ekb_source(0, 1) = 0.8;
    wg_source(0, 0) = 1.0;
    wg_source(0, 1) = 0.0;
    if (my_rank == 0)
    {
        ModuleIO::wfc_nao_write2file_complex(binary_test_dir + "wfk1_nao.dat",
                                             coefficients.data(),
                                             nlocal,
                                             0,
                                             ModuleBase::Vector3<double>(0.25, 0.0, 0.0),
                                             ekb_source,
                                             wg_source,
                                             true,
                                             false);
    }
#ifdef __MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    psi::Psi<std::complex<double>> psi_read;
    ModuleBase::matrix ekb(1, nbands);
    ModuleBase::matrix wg(1, nbands);
    const std::vector<int> ik2iktot = {0};
    EXPECT_TRUE(ModuleIO::read_wfc_nao(binary_test_dir,
                                       para,
                                       psi_read,
                                       ekb,
                                       wg,
                                       ik2iktot,
                                       1,
                                       1,
                                       true));
    EXPECT_DOUBLE_EQ(ekb(0, 0), ekb_source(0, 0));
    EXPECT_DOUBLE_EQ(ekb(0, 1), ekb_source(0, 1));
    EXPECT_DOUBLE_EQ(wg(0, 0), wg_source(0, 0));
    EXPECT_DOUBLE_EQ(wg(0, 1), wg_source(0, 1));
    if (my_rank == 0)
    {
        EXPECT_DOUBLE_EQ(psi_read(0, 0, 0).real(), coefficients[0].real());
        EXPECT_DOUBLE_EQ(psi_read(0, 0, 0).imag(), coefficients[0].imag());
    }
}

TEST_F(ReadWfcNaoTest, ReadBinaryPart)
{
    const int nbands_file = 3;
    const int nbands = 2;
    const int skip_band = 1;
    const int nlocal = 3;
    Parallel_Orbitals para;
    initialize_parallel_orbitals(para, nlocal, nbands);

    const std::vector<double> coefficients
        = {0.1, 0.2, 0.3, 1.1, 1.2, 1.3, 2.1, 2.2, 2.3};
    ModuleBase::matrix ekb_source(1, nbands_file);
    ModuleBase::matrix wg_source(1, nbands_file);
    for (int ib = 0; ib < nbands_file; ++ib)
    {
        ekb_source(0, ib) = -0.5 + ib;
        wg_source(0, ib) = 2.0 - ib;
    }
    if (my_rank == 0)
    {
        ModuleIO::wfc_nao_write2file(binary_test_dir + "wf_nao.dat",
                                     coefficients.data(),
                                     nlocal,
                                     0,
                                     ekb_source,
                                     wg_source,
                                     true,
                                     false);
    }
#ifdef __MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    psi::Psi<double> psi_read;
    ModuleBase::matrix ekb(1, nbands);
    ModuleBase::matrix wg(1, nbands);
    const std::vector<int> ik2iktot = {0};
    EXPECT_TRUE(ModuleIO::read_wfc_nao(binary_test_dir,
                                       para,
                                       psi_read,
                                       ekb,
                                       wg,
                                       ik2iktot,
                                       1,
                                       1,
                                       true,
                                       skip_band));
    EXPECT_DOUBLE_EQ(ekb(0, 0), ekb_source(0, 1));
    EXPECT_DOUBLE_EQ(ekb(0, 1), ekb_source(0, 2));
    EXPECT_DOUBLE_EQ(wg(0, 0), wg_source(0, 1));
    EXPECT_DOUBLE_EQ(wg(0, 1), wg_source(0, 2));
    if (my_rank == 0)
    {
        EXPECT_DOUBLE_EQ(psi_read(0, 0, 0), coefficients[nlocal]);
    }
}

TEST_F(ReadWfcNaoTest, BinaryDoesNotFallBackToText)
{
    const int nbands = 2;
    const int nlocal = 2;
    Parallel_Orbitals para;
    initialize_parallel_orbitals(para, nlocal, nbands);

    const std::vector<double> coefficients = {0.25, 0.5, 0.75, 1.0};
    ModuleBase::matrix ekb_source(1, nbands);
    ModuleBase::matrix wg_source(1, nbands);
    if (my_rank == 0)
    {
        ModuleIO::wfc_nao_write2file(binary_test_dir + "wf_nao.txt",
                                     coefficients.data(),
                                     nlocal,
                                     0,
                                     ekb_source,
                                     wg_source,
                                     false,
                                     false);
    }
#ifdef __MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    psi::Psi<double> psi_read;
    ModuleBase::matrix ekb(1, nbands);
    ModuleBase::matrix wg(1, nbands);
    const std::vector<int> ik2iktot = {0};
    EXPECT_FALSE(ModuleIO::read_wfc_nao(binary_test_dir,
                                        para,
                                        psi_read,
                                        ekb,
                                        wg,
                                        ik2iktot,
                                        1,
                                        1,
                                        true));
}

TEST_F(ReadWfcNaoTest, RejectTruncatedBinary)
{
    const int nbands = 2;
    const int nlocal = 2;
    Parallel_Orbitals para;
    initialize_parallel_orbitals(para, nlocal, nbands);

    if (my_rank == 0)
    {
        std::ofstream ofs(binary_test_dir + "wf_nao.dat", std::ios::binary);
        ofs.write(reinterpret_cast<const char*>(&nbands), sizeof(nbands));
        ofs.write(reinterpret_cast<const char*>(&nlocal), sizeof(nlocal));
    }
#ifdef __MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    psi::Psi<double> psi_read;
    ModuleBase::matrix ekb(1, nbands);
    ModuleBase::matrix wg(1, nbands);
    const std::vector<int> ik2iktot = {0};
    EXPECT_FALSE(ModuleIO::read_wfc_nao(binary_test_dir,
                                        para,
                                        psi_read,
                                        ekb,
                                        wg,
                                        ik2iktot,
                                        1,
                                        1,
                                        true));
}



TEST_F(ReadWfcNaoTest, CompleteSpinPopulations)
{
    const int nbands = 4;
    const int nlocal = 2;
    const int nspin = 2;
    const std::vector<int> ik2iktot = {0, 1};
    ModuleBase::matrix energies(2, nbands);
    ModuleBase::matrix occupations(2, nbands);
    const std::vector<double> real_coefficients(nbands * nlocal, 0.25);
    const std::vector<std::complex<double>> complex_coefficients(nbands * nlocal, {0.25, 0.125});
    for (int ib = 0; ib < nbands; ++ib)
    {
        occupations(0, ib) = ib < 3 ? 1.0 : 0.0;
        occupations(1, ib) = ib < 1 ? 1.0 : 0.0;
    }
    for (bool gamma_only : {true, false})
    {
        for (bool binary : {false, true})
        {
            const int file_type = binary ? 2 : 1;
            std::vector<std::string> files;
            for (int ik = 0; ik < 2; ++ik)
            {
                const std::string file = ModuleIO::filename_output(binary_test_dir, "wf", "nao",
                    ik, ik2iktot, nspin, 2, file_type, false, gamma_only, -1);
                files.push_back(file);
                if (my_rank == 0)
                {
                    if (gamma_only)
                    {
                        ModuleIO::wfc_nao_write2file(file, real_coefficients.data(), nlocal,
                            ik, energies, occupations, binary, false);
                    }
                    else
                    {
                        const ModuleBase::Vector3<double> kpoint(0.25, 0.0, 0.0);
                        ModuleIO::wfc_nao_write2file_complex(file, complex_coefficients.data(), nlocal,
                            ik, kpoint, energies, occupations, binary, false);
                    }
                }
            }
            std::vector<double> populations;
            ASSERT_TRUE(ModuleIO::read_wfc_nao_spin_populations(binary_test_dir, 2, nspin, gamma_only, binary, my_rank, populations));
            EXPECT_DOUBLE_EQ(populations[0], 3.0);
            EXPECT_DOUBLE_EQ(populations[1], 1.0);
            if (my_rank == 0)
            {
                for (const auto& file : files) { std::remove(file.c_str()); }
            }
        }
    }
}

TEST_F(ReadWfcNaoTest, CompleteGlobalKPointPopulations)
{
    const int nks = 4; // two k points per spin, independent of the reader's pool
    const int nbands = 4;
    const int nlocal = 2;
    const std::vector<int> global_indices = {0, 1, 2, 3};
    ModuleBase::matrix energies(nks, nbands);
    ModuleBase::matrix occupations(nks, nbands);
    const std::vector<std::complex<double>> coefficients(nbands * nlocal, {0.25, 0.125});
    std::vector<std::string> files;
    for (int ik = 0; ik < nks; ++ik)
    {
        for (int ib = 0; ib < nbands; ++ib)
        {
            const int occupied = ik < 2 ? 3 : 1;
            occupations(ik, ib) = ib < occupied ? 0.5 : 0.0;
        }
        const std::string file = ModuleIO::filename_output(binary_test_dir, "wf", "nao",
            ik, global_indices, 2, nks, 1, false, false, -1);
        files.push_back(file);
        if (my_rank == 0)
        {
            const ModuleBase::Vector3<double> kpoint(0.25 * ik, 0.0, 0.0);
            ModuleIO::wfc_nao_write2file_complex(file, coefficients.data(), nlocal,
                ik, kpoint, energies, occupations, false, false);
        }
    }
    std::vector<double> populations;
    ASSERT_TRUE(ModuleIO::read_wfc_nao_spin_populations(binary_test_dir, nks,
        2, false, false, my_rank, populations));
    EXPECT_DOUBLE_EQ(populations[0], 3.0);
    EXPECT_DOUBLE_EQ(populations[1], 1.0);
    EXPECT_FALSE(ModuleIO::read_wfc_nao_spin_populations(binary_test_dir, 0,
        2, false, false, my_rank, populations));
    if (my_rank == 0)
    {
        for (const auto& file : files) { std::remove(file.c_str()); }
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
