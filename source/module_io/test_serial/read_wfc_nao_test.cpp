#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "module_io/read_wfc_nao.h"
#include "module_basis/module_ao/parallel_orbitals.h"

//write mock function for Parallel_Orbitals
Parallel_Orbitals::Parallel_Orbitals(){}
Parallel_Orbitals::~Parallel_Orbitals(){}
/************************************************
 *  unit test of functions in read_wfc_nao.cpp
 ***********************************************/

/**
 * - Tested Functions:
 *   - distri_wfc_nao()
 *     - calculate memory required.
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
