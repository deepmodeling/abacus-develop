#include "gtest/gtest.h"
#include "gmock/gmock.h"
#define private public
#include "source_io/module_parameter/parameter.h"
#undef private
#include "source_io/read_wfc_nao.h"
#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_io/write_wfc_nao.h"

//define a mock derived class of class ElecState

namespace elecstate
{
      const double* ElecState::getRho(int spin) const{return &(this->eferm.ef);}//just for mock
}


namespace ModuleIO
{
// mock filename_output
std::string filename_output(
            const std::string &directory,
            const std::string &property,
            const std::string &basis,
            const int ik,
            const std::vector<int> &ik2iktot,
            const int nspin,
            const int nkstot,
            const int out_type,
            const bool out_app_flag,
            const bool gamma_only,
            const int istep,
            const int iter)
{
      return "./support/wfs1_nao.txt";
}
}

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
};


TEST_F(ReadWfcNaoTest,ReadWfcNao)
{
      //Global variables
      const int nbands = 3;
      const int nlocal = 3;
      PARAM.sys.global_readin_dir = "./support/";
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
      elecstate::ElecState pelec;
      pelec.ekb.create(nks,nbands);
      pelec.wg.create(nks,nbands);

      std::vector<int> ik2iktot = {0};
      const int nkstot = 1;

      // Act
	  ModuleIO::read_wfc_nao(PARAM.sys.global_readin_dir, ParaV, psid, 
			  &(pelec), ik2iktot, nkstot, nspin);
      // Assert
      EXPECT_NEAR(pelec.ekb(0,1),0.31482195194888534794941393,1e-5);
      EXPECT_NEAR(pelec.wg(0,1),0.0,1e-5);
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
    PARAM.sys.global_readin_dir = "./support/";
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
    elecstate::ElecState pelec;
    pelec.ekb.create(nks, nbands);
    pelec.wg.create(nks, nbands);

	std::vector<int> ik2iktot = {0};
	const int nkstot = 1;

	// Act
	ModuleIO::read_wfc_nao(PARAM.sys.global_readin_dir, ParaV, psid, 
			&(pelec), ik2iktot, nkstot, nspin, skip_band, nstep);

    // Assert
    EXPECT_NEAR(pelec.ekb(0, 1), 7.4141254894954844445464914e-01, 1e-5);
    if (my_rank == 0)
    {
        EXPECT_NEAR(psid(0, 0, 0), 1.8587183851, 1e-5);
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
