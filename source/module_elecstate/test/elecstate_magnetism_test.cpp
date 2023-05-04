#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include <cmath>

/************************************************
 *  unit test of magnetism.cpp
 ***********************************************/

/**
 * - Tested Functions:
 *   - Magnetism::Magnetism()
 *   - Magnetism::~Magnetism()
 *   - Magnetism::judge_parallel()
 *   - Magnetism::compute_magnetization()
 *      - compute magnetization for spin-polarized system when GlobalV::NSPIN = 2
 *      - and non-collinear case with GlobalV::NSPIN = 4 
*/

#define private public
#include "module_elecstate/magnetism.h"

Charge::Charge(){}
Charge::~Charge(){}

int Magnetism::get_rhopw_nrxx() const { return 100; }
int Magnetism::get_rhopw_nxyz() const { return 1000; }
double Magnetism::get_ucell_omega() const { return 500.0; }

class MagnetismTest : public ::testing::Test
{
protected:
    Magnetism* magnetism;
    virtual void SetUp()
    {
        magnetism = new Magnetism;
    }
    virtual void TearDown()
    {
        delete magnetism;
    }
};

TEST_F(MagnetismTest, Magnetism)
{
                  EXPECT_EQ(0.0, magnetism->tot_magnetization);
                  EXPECT_EQ(0.0, magnetism->abs_magnetization);
                  EXPECT_EQ(nullptr, magnetism->start_magnetization);
}

TEST_F(MagnetismTest, GlobalInfo)
{
                  EXPECT_EQ(100, magnetism->get_rhopw_nrxx());
                  EXPECT_EQ(1000, magnetism->get_rhopw_nxyz());
                  EXPECT_EQ(500.0, magnetism->get_ucell_omega());
}

TEST_F(MagnetismTest, JudgeParallel)
{
                  double a[3] = {1.0, 0.0, 0.0};
                  ModuleBase::Vector3<double> b(1.0, 0.0, 0.0);
                  EXPECT_TRUE(magnetism->judge_parallel(a, b));
                  b = ModuleBase::Vector3<double>(0.0, 1.0, 0.0);
                  EXPECT_FALSE(magnetism->judge_parallel(a, b));
}

TEST_F(MagnetismTest, ComputeMagnetizationS2)
{
                  GlobalV::NSPIN = 2;
                  GlobalV::TWO_EFERMI = false;
                  GlobalV::nelec = 10.0;
                  Charge* chr = new Charge;
                  chr->rho = new double*[GlobalV::NSPIN];
                  for (int i=0; i< GlobalV::NSPIN; i++)
                  {
                                    chr->rho[i] = new double[magnetism->get_rhopw_nrxx()];
                  }
                  for (int ir=0; ir< magnetism->get_rhopw_nrxx(); ir++)
                  {
                                    chr->rho[0][ir] = 1.00;
                                    chr->rho[1][ir] = 1.01;
                  }
                  double* nelec_spin = new double[2];
                  magnetism->compute_magnetization(chr, nelec_spin);
                  EXPECT_DOUBLE_EQ(-0.5, magnetism->tot_magnetization);
                  EXPECT_DOUBLE_EQ(0.5, magnetism->abs_magnetization);
                  EXPECT_DOUBLE_EQ(4.75, nelec_spin[0]);
                  EXPECT_DOUBLE_EQ(5.25, nelec_spin[1]);
                  delete[] nelec_spin;
                  for (int i=0; i< GlobalV::NSPIN; i++)
                  {
                                    delete[] chr->rho[i];
                  }
                  delete[] chr->rho;
                  delete chr;
}


TEST_F(MagnetismTest, ComputeMagnetizationS4)
{
                    GlobalV::NSPIN = 4;
                    Charge* chr = new Charge;
                    chr->rho = new double*[GlobalV::NSPIN];
                    for (int i=0; i< GlobalV::NSPIN; i++)
                    {
                                        chr->rho[i] = new double[magnetism->get_rhopw_nrxx()];
                    }
                    for (int ir=0; ir< magnetism->get_rhopw_nrxx(); ir++)
                    {
                                        chr->rho[0][ir] = 1.00;
                                        chr->rho[1][ir] = std::sqrt(2.0);
                                        chr->rho[2][ir] = 1.00;
                                        chr->rho[3][ir] = 1.00;
                    }
                    double* nelec_spin = new double[4];
                    magnetism->compute_magnetization(chr, nelec_spin);
                    EXPECT_DOUBLE_EQ(100.0, magnetism->abs_magnetization);
                    EXPECT_DOUBLE_EQ(50.0*std::sqrt(2.0), magnetism->tot_magnetization_nc[0]);
                    EXPECT_DOUBLE_EQ(50.0, magnetism->tot_magnetization_nc[1]);
                    EXPECT_DOUBLE_EQ(50.0, magnetism->tot_magnetization_nc[2]);
                    delete[] nelec_spin;
                    for (int i=0; i< GlobalV::NSPIN; i++)
                    {
                                        delete[] chr->rho[i];
                    }
                    delete[] chr->rho;
                    delete chr;
}


#ifdef __MPI
#include <mpi.h>
int main(int argc, char **argv)
{

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD,&GlobalV::NPROC);
    MPI_Comm_rank(MPI_COMM_WORLD,&GlobalV::MY_RANK);

    testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();

    MPI_Finalize();

    return result;
}
#endif

#undef private
