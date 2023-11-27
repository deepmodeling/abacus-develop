#include <algorithm>
#include <random>

#include "../sc_lambda_lcao.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "module_cell/klist.h"
#include "module_hamilt_lcao/module_deltaspin/spin_constrain.h"

// mockes
K_Vectors::K_Vectors()
{
}
K_Vectors::~K_Vectors()
{
}

template <>
void hamilt::OperatorLCAO<std::complex<double>, std::complex<double>>::init(const int ik_in)
{
}

template <>
void hamilt::OperatorLCAO<std::complex<double>, double>::init(const int ik_in)
{
}

template <>
void hamilt::OperatorLCAO<double, double>::init(const int ik_in)
{
}

template <>
void hamilt::OperatorLCAO<std::complex<double>, std::complex<double>>::contributeHk(const int ik_in)
{
}

template <>
void hamilt::OperatorLCAO<std::complex<double>, double>::contributeHk(const int ik_in)
{
}

template <>
void hamilt::OperatorLCAO<double, double>::contributeHk(const int ik_in)
{
}
// mocks

/************************************************
 *  unit test of class OperatorScLambda
 ***********************************************/

/**
 * - Tested Functions:
 *   -
 */

LCAO_Matrix::LCAO_Matrix()
{
}

LCAO_Matrix::~LCAO_Matrix()
{
}

class ScLambdaLCAOTest : public ::testing::Test
{
  protected:
    std::vector<ModuleBase::Vector3<double>> kvec_d;
    std::vector<int> isk;
    void SetUp()
    {
        kvec_d.push_back(ModuleBase::Vector3<double>(0.0, 0.0, 0.0));
        kvec_d.push_back(ModuleBase::Vector3<double>(0.5, 0.5, 0.5));
        isk.push_back(0);
        isk.push_back(1);
    }
};

TEST_F(ScLambdaLCAOTest, ContributeHk)
{
    // set paraV
    Parallel_Orbitals paraV;
    std::ofstream ofs("test.log");
    int nrow = 2;
    int ncol = 2;
    paraV.set_global2local(nrow, ncol, false, ofs);
    EXPECT_EQ(paraV.nloc, 4);
    // set LM
    LCAO_Matrix LM;
    LM.ParaV = &paraV;
    EXPECT_EQ(LM.ParaV->nloc, 4);
    LM.Sloc2 = {
        std::complex<double>{1.0, 0.0},
        std::complex<double>{0.0, 0.0},
        std::complex<double>{0.0, 0.0},
        std::complex<double>{1.0, 0.0}
    };
    LM.Hloc2.resize(LM.ParaV->nloc, 0.0);
    // set sc
    SpinConstrain<std::complex<double>, psi::DEVICE_CPU>& sc
        = SpinConstrain<std::complex<double>, psi::DEVICE_CPU>::getScInstance();
    sc.set_ParaV(&paraV);
    EXPECT_EQ(sc.ParaV->nloc, 4);
    std::map<int, int> atomCounts = {
        {0, 1}
    };
    std::map<int, int> orbitalCounts = {
        {0, 1}
    };
    sc.set_atomCounts(atomCounts);
    sc.set_orbitalCounts(orbitalCounts);
    sc.set_nspin(4);
    sc.set_npol(2);
    ModuleBase::Vector3<double>* sc_lambda = new ModuleBase::Vector3<double>[1];
    sc_lambda[0][0] = 1.0;
    sc_lambda[0][1] = 1.0;
    sc_lambda[0][2] = 1.0;
    sc.set_sc_lambda(sc_lambda, 1);
    // set KS_SOLVER, which determines IS_COLUMN_MAJOR_KS_SOLVER()
    GlobalV::KS_SOLVER = "genelpa";
    EXPECT_TRUE(ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER());
    // set sc_lambda_op
    auto sc_lambda_op
        = hamilt::OperatorScLambda<hamilt::OperatorLCAO<std::complex<double>, std::complex<double>>>(&LM,
                                                                                                     this->kvec_d,
                                                                                                     nullptr,
                                                                                                     nullptr,
                                                                                                     isk);
    sc_lambda_op.contributeHk(0);
    std::vector<std::complex<double>> columnMajor_h_lambda = {
        std::complex<double>{-1.0, 0.0},
        std::complex<double>{-1.0, 1.0},
        std::complex<double>{-1.0, -1.0},
        std::complex<double>{1.0,  0.0}
    };
    EXPECT_DOUBLE_EQ(LM.Hloc2[0].real(), columnMajor_h_lambda[0].real());
    EXPECT_DOUBLE_EQ(LM.Hloc2[0].imag(), columnMajor_h_lambda[0].imag());
    EXPECT_DOUBLE_EQ(LM.Hloc2[1].real(), columnMajor_h_lambda[1].real());
    EXPECT_DOUBLE_EQ(LM.Hloc2[1].imag(), columnMajor_h_lambda[1].imag());
    EXPECT_DOUBLE_EQ(LM.Hloc2[2].real(), columnMajor_h_lambda[2].real());
    EXPECT_DOUBLE_EQ(LM.Hloc2[2].imag(), columnMajor_h_lambda[2].imag());
    EXPECT_DOUBLE_EQ(LM.Hloc2[3].real(), columnMajor_h_lambda[3].real());
    EXPECT_DOUBLE_EQ(LM.Hloc2[3].imag(), columnMajor_h_lambda[3].imag());
}

TEST_F(ScLambdaLCAOTest, ContributeHkS2)
{
    // set paraV
    Parallel_Orbitals paraV;
    std::ofstream ofs("test.log");
    int nrow = 1;
    int ncol = 1;
    paraV.set_global2local(nrow, ncol, false, ofs);
    EXPECT_EQ(paraV.nloc, 1);
    // set LM
    LCAO_Matrix LM;
    LM.ParaV = &paraV;
    EXPECT_EQ(LM.ParaV->nloc, 1);
    LM.Sloc2 = {
        std::complex<double>{1.0, 0.0}
    };
    LM.Hloc2.resize(LM.ParaV->nloc, 0.0);
    // set sc
    SpinConstrain<std::complex<double>, psi::DEVICE_CPU>& sc
        = SpinConstrain<std::complex<double>, psi::DEVICE_CPU>::getScInstance();
    sc.set_ParaV(&paraV);
    EXPECT_EQ(sc.ParaV->nloc, 1);
    std::map<int, int> atomCounts = {
        {0, 1}
    };
    std::map<int, int> orbitalCounts = {
        {0, 1}
    };
    sc.set_atomCounts(atomCounts);
    sc.set_orbitalCounts(orbitalCounts);
    sc.set_nspin(2);
    sc.set_npol(1);
    ModuleBase::Vector3<double>* sc_lambda = new ModuleBase::Vector3<double>[1];
    sc_lambda[0][0] = 0.0;
    sc_lambda[0][1] = 0.0;
    sc_lambda[0][2] = 1.0;
    sc.set_sc_lambda(sc_lambda, 1);
    // set KS_SOLVER, which determines IS_COLUMN_MAJOR_KS_SOLVER()
    GlobalV::KS_SOLVER = "genelpa";
    EXPECT_TRUE(ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER());
    // set sc_lambda_op
    auto sc_lambda_op
        = hamilt::OperatorScLambda<hamilt::OperatorLCAO<std::complex<double>, double>>(&LM,
                                                                                                     this->kvec_d,
                                                                                                     nullptr,
                                                                                                     nullptr,
                                                                                                     isk);
    sc_lambda_op.contributeHk(0);
    std::vector<std::complex<double>> columnMajor_h_lambda = {
        std::complex<double>{-1.0,  0.0}
    };
    EXPECT_DOUBLE_EQ(LM.Hloc2[0].real(), columnMajor_h_lambda[0].real());
    EXPECT_DOUBLE_EQ(LM.Hloc2[0].imag(), columnMajor_h_lambda[0].imag());
}

TEST_F(ScLambdaLCAOTest, TemplateHelpers)
{
    auto sc_lambda_op
        = hamilt::OperatorScLambda<hamilt::OperatorLCAO<std::complex<double>, std::complex<double>>>(nullptr,
                                                                                                     this->kvec_d,
                                                                                                     nullptr,
                                                                                                     nullptr,
                                                                                                     isk);
    auto sc_lambda_op1 = hamilt::OperatorScLambda<hamilt::OperatorLCAO<std::complex<double>, double>>(nullptr,
                                                                                                      this->kvec_d,
                                                                                                      nullptr,
                                                                                                      nullptr,
                                                                                                      isk);
    auto sc_lambda_op2
        = hamilt::OperatorScLambda<hamilt::OperatorLCAO<double, double>>(nullptr, this->kvec_d, nullptr, nullptr, isk);
    EXPECT_NO_THROW(sc_lambda_op.contributeHR());
    EXPECT_NO_THROW(sc_lambda_op1.contributeHR());
    EXPECT_NO_THROW(sc_lambda_op2.contributeHR());
    EXPECT_NO_THROW(sc_lambda_op2.contributeHk(0));
}