#include "../spin_constrain.h"

#include <algorithm>
#include <string>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

/************************************************
 *  unit test of functions in class SpinConstrain
 ***********************************************/

/**
 * - Tested functions:
 *  - spinconstrain::SpinConstrain::getScInstance()
 *      get the instance of spinconstrain::SpinConstrain
 *  - spinconstrain::SpinConstrain::set_atomCounts()
 *     set the map from element index to atom number
 *  - spinconstrain::SpinConstrain::get_atomCounts()
 *     get the map from element index to atom number
 *  - spinconstrain::SpinConstrain::get_nat()
 *     get the total number of atoms
 *  - spinconstrain::SpinConstrain::get_iat()
 *     get the atom index from (itype, atom_index)
 *  - spinconstrain::SpinConstrain::set_orbitalCounts()
 *     set the map from element index to orbital number
 *  - spinconstrain::SpinConstrain::get_orbitalCounts()
 *     get the map from element index to orbital number
 *  - spinconstrain::SpinConstrain::get_nw()
 *     get the total number of orbitals
 *  - spinconstrain::SpinConstrain::set_npol()
 *     set the number of npol, which is the number of spin components
 *  - spinconstrain::SpinConstrain::get_npol()
 *     get the number of npol, which is the number of spin components
 *  - spinconstrain::SpinConstrain::get_iwt()
 *     get the index of orbital with spin component from (itype, iat, orbital_index)
 */
#include "module_cell/klist.h"
K_Vectors::K_Vectors(){}
K_Vectors::~K_Vectors(){}

template <typename T>
class SpinConstrainTest : public testing::Test
{
  protected:
    spinconstrain::SpinConstrain<T>& sc = spinconstrain::SpinConstrain<T>::getScInstance();
};

using MyTypes = ::testing::Types<double, std::complex<double>>;
TYPED_TEST_SUITE(SpinConstrainTest, MyTypes);

TYPED_TEST(SpinConstrainTest, CheckAtomCounts)
{
    // Warning 1: atomCounts is not set
    testing::internal::CaptureStdout();
    EXPECT_EXIT(this->sc.check_atomCounts(), ::testing::ExitedWithCode(0), "");
    std::string output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("atomCounts is not set"));
    // Warning 2: nat < 0
    std::map<int, int> atomCounts = {
        {0, -1},
        {1, 0 }
    };
    this->sc.set_atomCounts(atomCounts);
    testing::internal::CaptureStdout();
    EXPECT_EXIT(this->sc.check_atomCounts(), ::testing::ExitedWithCode(0), "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("nat <= 0"));
    // Warning 3: itype out of range
    std::map<int, int> atomCounts1 = {
        {1, 1},
        {2, 2}
    };
    this->sc.set_atomCounts(atomCounts1);
    testing::internal::CaptureStdout();
    EXPECT_EXIT(this->sc.check_atomCounts(), ::testing::ExitedWithCode(0), "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("itype out of range [0, ntype)"));
    // Warning 4: number of atoms <= 0 for some element
    std::map<int, int> atomCounts2 = {
        {0, 2 },
        {1, -1}
    };
    this->sc.set_atomCounts(atomCounts2);
    testing::internal::CaptureStdout();
    EXPECT_EXIT(this->sc.check_atomCounts(), ::testing::ExitedWithCode(0), "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("number of atoms <= 0 for some element"));
}

TYPED_TEST(SpinConstrainTest, AtomCounts)
{
    std::map<int, int> atomCounts = {
        {0, 5 },
        {1, 10}
    };
    this->sc.set_atomCounts(atomCounts);
    std::map<int, int> atomCounts2 = this->sc.get_atomCounts();
    int ntype = atomCounts2.size();
    EXPECT_EQ(ntype, 2);
    int nat = this->sc.get_nat();
    EXPECT_EQ(nat, 15);
    EXPECT_EQ(this->sc.get_iat(1, 4), 9); // atom_index starts from 0
    // warning 1: itype out of range
    testing::internal::CaptureStdout();
    EXPECT_EXIT(this->sc.get_iat(3, 0);, ::testing::ExitedWithCode(0), "");
    std::string output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("itype out of range [0, ntype)"));
    // warning 2: atom_index out of range
    testing::internal::CaptureStdout();
    EXPECT_EXIT(this->sc.get_iat(0, 5);, ::testing::ExitedWithCode(0), "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("atom index out of range [0, nat)"));
}

TYPED_TEST(SpinConstrainTest, NSPIN)
{
    this->sc.set_nspin(4);
    int nspin = this->sc.get_nspin();
    EXPECT_EQ(nspin, 4);
}

TYPED_TEST(SpinConstrainTest, NSPINwarning)
{
    testing::internal::CaptureStdout();
    EXPECT_EXIT(this->sc.set_nspin(1), ::testing::ExitedWithCode(0), "");
    std::string output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("nspin must be 2 or 4"));
}

TYPED_TEST(SpinConstrainTest, SetInputParameters)
{
    double sc_thr = 1e-6;
    int nsc = 100;
    int nsc_min = 2;
    double alpha_trial = 0.01;
    double sccut = 3.0;
    double sc_drop_thr = 1e-3;
    this->sc.set_input_parameters(sc_thr, nsc, nsc_min, alpha_trial, sccut, sc_drop_thr);
    EXPECT_DOUBLE_EQ(this->sc.get_sc_thr(), sc_thr);
    EXPECT_EQ(this->sc.get_nsc(), nsc);
    EXPECT_EQ(this->sc.get_nsc_min(), nsc_min);
    EXPECT_DOUBLE_EQ(this->sc.get_alpha_trial(), alpha_trial / 13.605698);
    EXPECT_DOUBLE_EQ(this->sc.get_sccut(), sccut / 13.605698);
    EXPECT_EQ(this->sc.get_sc_drop_thr(), sc_drop_thr);
}

TYPED_TEST(SpinConstrainTest, SetSolverParameters)
{
    K_Vectors kv;
    this->sc.set_nspin(4);
    this->sc.set_solver_parameters(kv, nullptr, nullptr, nullptr, "genelpa");
    EXPECT_EQ(this->sc.get_nspin(), 4);
    EXPECT_EQ(this->sc.p_hamilt, nullptr);
    EXPECT_EQ(this->sc.psi, nullptr);
    EXPECT_EQ(this->sc.pelec, nullptr);
    EXPECT_EQ(this->sc.KS_SOLVER, "genelpa");
}

TYPED_TEST(SpinConstrainTest, SetParaV)
{
    Parallel_Orbitals paraV;
    // warning 1
    paraV.nloc = 0;
    testing::internal::CaptureStdout();
    EXPECT_EXIT(this->sc.set_ParaV(&paraV), ::testing::ExitedWithCode(0), "");
    std::string output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("nloc <= 0"));
    // normal set
    int nrow = 4;
    int ncol = 4;
    std::ofstream ofs("test.log");
    paraV.set_serial(nrow, ncol);
    this->sc.set_ParaV(&paraV);
    EXPECT_EQ(this->sc.ParaV->nloc, nrow * ncol);
    remove("test.log");
}

/*
TYPED_TEST(SpinConstrainTest, PrintMi)
{
    this->sc.zero_Mi();
    testing::internal::CaptureStdout();
    this->sc.print_Mi(true);
    std::string output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("Total Magnetism (uB):"));
    EXPECT_THAT(output, testing::HasSubstr("ATOM      0         0.0000000000         0.0000000000         0.0000000000"));
    this->sc.set_nspin(2);
     testing::internal::CaptureStdout();
    this->sc.print_Mi(true);
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("Total Magnetism (uB):"));
    EXPECT_THAT(output, testing::HasSubstr("ATOM      0         0.0000000000"));
}
*/
