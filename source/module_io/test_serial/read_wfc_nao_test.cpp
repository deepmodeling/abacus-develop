#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "module_io/read_wfc_nao.h"
#include "module_basis/module_ao/parallel_orbitals.h"

//write mock function for Parallel_Orbitals
Parallel_2D::Parallel_2D(){}
Parallel_2D::~Parallel_2D(){}
Parallel_Orbitals::Parallel_Orbitals() {}
Parallel_Orbitals::~Parallel_Orbitals(){}
//define a mock derived class of class ElecState

namespace elecstate
{
      const double* ElecState::getRho(int spin) const{return &(this->eferm.ef);}//just for mock
      void ElecState::calculate_weights(){}
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

TEST_F(ReadWfcNaoTest,DistriWfcNao)
{
      // Arrange
      int is = 0;
      int nks = 1;
      GlobalV::NBANDS = 2;
      GlobalV::NLOCAL = 3;
      int nband = GlobalV::NBANDS;
      int nlocal = GlobalV::NLOCAL;
      int ngk[1] = {10};
      double** ctot = new double*[nband];
      for (int i=0; i<nband; i++)
      {
                        ctot[i] = new double[nlocal];
                        for (int j=0; j<nlocal; j++)
                        {
                                          ctot[i][j] = i*nlocal + j;
                        }
      }
      Parallel_Orbitals* ParaV = new Parallel_Orbitals;
      psi::Psi<double>* psid = new psi::Psi<double>(nks, nband, nlocal, &ngk[0]);
      // Act
      ModuleIO::distri_wfc_nao(ctot, is, ParaV, psid);
      // Assert
      for (int i=0; i<nband; i++)
      {
                        for (int j=0; j<nlocal; j++)
                        {
                                          EXPECT_DOUBLE_EQ(ctot[i][j], psid[0](is, i, j));
                        }
      }
      // Clean up
      for (int i=0; i<nband; i++)
      {
                        delete[] ctot[i];
      }
      delete[] ctot;
      delete ParaV;
      delete psid;
}

TEST_F(ReadWfcNaoTest,ReadWfcNao)
{
      //Global variables
      GlobalV::NBANDS = 3;
      GlobalV::NLOCAL = 3;
      GlobalV::GAMMA_ONLY_LOCAL = true;
      GlobalV::global_readin_dir = "./support/";
      GlobalV::DRANK = 0;
      GlobalV::MY_RANK = 0;
      // Arrange
      int is = 0;
      int nks = 1;
      int nband = GlobalV::NBANDS;
      int nlocal = GlobalV::NLOCAL;
      int ngk[1] = {1};
      double** ctot;
      Parallel_Orbitals* ParaV = new Parallel_Orbitals;
      psi::Psi<double>* psid = new psi::Psi<double>(nks, nband, nlocal, &ngk[0]);
      elecstate::ElecState* pelec = new elecstate::ElecState;
      //elecstate::ElecState* pelec = new elecstate::MockElecState;
      pelec->ekb.create(nks,nband);
      pelec->wg.create(nks,nband);
      // Act
      ModuleIO::read_wfc_nao(ctot, is, ParaV, psid, pelec);
      // Assert
      EXPECT_NEAR(pelec->ekb(0,1),0.314822,1e-5);
      EXPECT_NEAR(pelec->wg(0,1),0.0,1e-5);
      EXPECT_NEAR(psid[0](0,1,1),0.595957,1e-5);
      // clean up
      delete ParaV;
      delete psid;
      delete pelec;
}
