#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "module_io/read_wfc_nao.h"
#include "module_basis/module_ao/parallel_orbitals.h"
#include "module_io/write_wfc_nao.h"

//define a mock derived class of class ElecState

namespace elecstate
{
      const double* ElecState::getRho(int spin) const{return &(this->eferm.ef);}//just for mock
      void ElecState::calculate_weights(){}
}

// mock wfc_lcao_gen_fname
std::string ModuleIO::wfc_nao_gen_fname(const int out_type,
                                         const bool gamma_only,
                                         const bool out_app_flag,
                                         const int ik,
                                         const int istep)
{
      return "WFC_NAO_GAMMA1.txt";
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
      int nbands = 3;
      int nlocal = 3;
      GlobalV::global_readin_dir = "./support/";
      int nks = 1;

      Parallel_Orbitals ParaV;
#ifdef __MPI
      ParaV.init(nlocal, nlocal, 1);      
      ParaV.set_desc_wfc_Eij(nlocal, nbands, pv->nrow);
#else
      ParaV.set_serial(nlocal, nlocal);
      ParaV.nrow_bands = nlocal;
      ParaV.ncol_bands = nbands;  
#endif 

      psi::Psi<double> psid;
      elecstate::ElecState pelec;
      pelec.ekb.create(nks,nbands);
      pelec.wg.create(nks,nbands);
      // Act
      ModuleIO::read_wfc_nao(GlobalV::global_readin_dir, ParaV, psid, &(pelec));
      // Assert
      EXPECT_NEAR(pelec.ekb(0,1),0.31482195194888534794941393,1e-5);
      EXPECT_NEAR(pelec.wg(0,1),0.0,1e-5);
      EXPECT_NEAR(psid(0,1,1),0.59595655928,1e-5);
}
TEST_F(ReadWfcNaoTest, ReadWfcNaoPart)
{
    //Global variables
    const int nbands = 2;
    const int skip_band = 1;
    const int nlocal = 3;
    GlobalV::global_readin_dir = "./support/";
    const int nks = 1;

    Parallel_Orbitals ParaV;
#ifdef __MPI
    ParaV.init(nlocal, nlocal, 1);
    ParaV.set_desc_wfc_Eij(nlocal, nbands, pv->nrow);
#else
    ParaV.set_serial(nlocal, nlocal);
    ParaV.nrow_bands = nlocal;
    ParaV.ncol_bands = nbands;
#endif 

    psi::Psi<double> psid;
    elecstate::ElecState pelec;
    pelec.ekb.create(nks, nbands);
    pelec.wg.create(nks, nbands);
    // Act
    ModuleIO::read_wfc_nao(GlobalV::global_readin_dir, ParaV, psid, &(pelec), /*skip_band=*/1);
    // Assert
    EXPECT_NEAR(pelec.ekb(0, 1), 7.4141254894954844445464914e-01, 1e-5);
    EXPECT_NEAR(psid(0, 1, 1), 5.6401881891e-01, 1e-5);
}
