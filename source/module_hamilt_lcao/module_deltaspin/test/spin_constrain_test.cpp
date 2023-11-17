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
 *  - SpinConstrain::getScInstance()
 *      get the instance of SpinConstrain
 *  - SpinConstrain::set_atomCounts()
 *     set the map from element index to atom number
 *  - SpinConstrain::get_atomCounts()
 *     get the map from element index to atom number
 *  - SpinConstrain::get_nat()
 *     get the total number of atoms
 *  - SpinConstrain::get_iat()
 *     get the atom index from (itype, atom_index)
 *  - SpinConstrain::set_orbitalCounts()
 *     set the map from element index to orbital number
 *  - SpinConstrain::get_orbitalCounts()
 *     get the map from element index to orbital number
 *  - SpinConstrain::get_nw()
 *     get the total number of orbitals
 *  - SpinConstrain::set_npol()
 *     set the number of npol, which is the number of spin components
 *  - SpinConstrain::get_npol()
 *     get the number of npol, which is the number of spin components
 *  - SpinConstrain::get_iwt()
 *     get the index of orbital with spin component from (itype, iat, orbital_index)
 */

K_Vectors::K_Vectors(){}
K_Vectors::~K_Vectors(){}

template <typename T>
class SpinConstrainTest : public testing::Test
{
  protected:
    SpinConstrain<T, psi::DEVICE_CPU>& sc = SpinConstrain<T, psi::DEVICE_CPU>::getScInstance();
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

TYPED_TEST(SpinConstrainTest, OrbitalCounts)
{
	std::map<int, int> orbitalCounts = {{0,5},{1,10}};
    std::map<int, int> atomCounts = {
        {0, 1},
        {1, 2}
    };
    this->sc.set_atomCounts(atomCounts);
    this->sc.set_orbitalCounts(orbitalCounts);
    std::map<int, int> orbitalCounts2 = this->sc.get_orbitalCounts();
    int ntype = orbitalCounts2.size();
    EXPECT_EQ(ntype, 2);
    EXPECT_EQ(this->sc.get_nw(), 25);
    this->sc.set_npol(2);
    EXPECT_EQ(this->sc.get_npol(), 2);
    EXPECT_EQ(this->sc.get_nw(), 50); // npol = 2
    this->sc.set_npol(1);
    EXPECT_EQ(this->sc.get_npol(), 1);
    EXPECT_EQ(this->sc.get_iwt(1, 1, 2), 17);
    this->sc.set_npol(2);
    EXPECT_EQ(this->sc.get_iwt(1, 1, 2), 32); // npol = 2
    // warning 1: itype out of range
    testing::internal::CaptureStdout();
    EXPECT_EXIT(this->sc.get_iwt(3, 0, 0);, ::testing::ExitedWithCode(0), "");
    std::string output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("itype out of range [0, ntype)"));
    // warning 2: atom_index out of range
    testing::internal::CaptureStdout();
    EXPECT_EXIT(this->sc.get_iwt(0, 3, 0);, ::testing::ExitedWithCode(0), "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("iat out of range [0, nat)"));
    // warning 3: orbital_index out of range
    testing::internal::CaptureStdout();
    EXPECT_EXIT(this->sc.get_iwt(0, 0, 10);, ::testing::ExitedWithCode(0), "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("orbital index out of range [0, atom_nw*npol)"));
}

TYPED_TEST(SpinConstrainTest, SetScLambdaMagConstrain)
{
    this->sc.Set_ScData_From_Json("./support/sc_f1.json");
    std::map<int, int> atomCounts = {
        {0, 5 },
        {1, 10}
    };
    this->sc.set_atomCounts(atomCounts);
    int nat = this->sc.get_nat();
    this->sc.set_sc_lambda();
    this->sc.set_target_mag();
    this->sc.set_constrain();
    std::vector<ModuleBase::Vector3<double>> sc_lambda = this->sc.get_sc_lambda();
    std::vector<ModuleBase::Vector3<double>> target_mag = this->sc.get_target_mag();
    std::vector<ModuleBase::Vector3<int>> constrain = this->sc.get_constrain();
    this->sc.set_sc_lambda(sc_lambda.data(), nat);
    this->sc.set_target_mag(target_mag.data(), nat);
    this->sc.set_constrain(constrain.data(), nat);
    EXPECT_EQ(sc_lambda.size(), this->sc.get_nat());
    for (const auto& sc_elem: this->sc.get_ScData())
    {
        int itype = sc_elem.first;
        const std::vector<ScAtomData>& sc_atoms = sc_elem.second;
        for (const ScAtomData& sc_data: sc_atoms)
        {
            int index = sc_data.index;
            int iat = this->sc.get_iat(itype, index);
            EXPECT_DOUBLE_EQ(sc_data.lambda[0] * this->sc.meV_to_Ry, sc_lambda[iat].x);
            EXPECT_DOUBLE_EQ(sc_data.lambda[1] * this->sc.meV_to_Ry, sc_lambda[iat].y);
            EXPECT_DOUBLE_EQ(sc_data.lambda[2] * this->sc.meV_to_Ry, sc_lambda[iat].z);
            EXPECT_DOUBLE_EQ(sc_data.target_mag[0], target_mag[iat].x);
            EXPECT_DOUBLE_EQ(sc_data.target_mag[1], target_mag[iat].y);
            EXPECT_DOUBLE_EQ(sc_data.target_mag[2], target_mag[iat].z);
            EXPECT_EQ(sc_data.constrain[0], constrain[iat].x);
            EXPECT_EQ(sc_data.constrain[1], constrain[iat].y);
            EXPECT_EQ(sc_data.constrain[2], constrain[iat].z);
        }
    }
    for (int iat = 0; iat < this->sc.get_nat(); iat++)
    {
        if (!(iat == 1 || iat == 5 || iat == 9))
        {
            EXPECT_DOUBLE_EQ(sc_lambda[iat].x, 0.0);
            EXPECT_DOUBLE_EQ(sc_lambda[iat].y, 0.0);
            EXPECT_DOUBLE_EQ(sc_lambda[iat].z, 0.0);
            EXPECT_DOUBLE_EQ(target_mag[iat].x, 0.0);
            EXPECT_DOUBLE_EQ(target_mag[iat].y, 0.0);
            EXPECT_DOUBLE_EQ(target_mag[iat].z, 0.0);
        }
    }
    // set_sc_lambda warning
    testing::internal::CaptureStdout();
    EXPECT_EXIT(this->sc.set_sc_lambda(sc_lambda.data(), 100);, ::testing::ExitedWithCode(0), "");
    std::string output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("lambda_in size mismatch with nat"));
    // set_target_mag warning
    testing::internal::CaptureStdout();
    EXPECT_EXIT(this->sc.set_target_mag(target_mag.data(), 100);, ::testing::ExitedWithCode(0), "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("target_mag_in size mismatch with nat"));
    // set constrain warning
    testing::internal::CaptureStdout();
    EXPECT_EXIT(this->sc.set_constrain(constrain.data(), 100);, ::testing::ExitedWithCode(0), "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("constrain_in size mismatch with nat"));
}

TYPED_TEST(SpinConstrainTest, CalEscon)
{
    this->sc.zero_Mi();
    this->sc.Set_ScData_From_Json("./support/sc_f1.json");
    std::map<int, int> atomCounts = {
        {0, 5 },
        {1, 10}
    };
    this->sc.set_atomCounts(atomCounts);
    int nat = this->sc.get_nat();
    this->sc.set_sc_lambda();
    double escon = this->sc.cal_escon();
    double escon1 = this->sc.get_escon();
    EXPECT_DOUBLE_EQ(escon, escon1);
    EXPECT_DOUBLE_EQ(escon1, 0.0);
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

TYPED_TEST(SpinConstrainTest, SetScDecayGrad)
{
    std::map<int, int> atomCounts = {
        {0, 1},
        {1, 5}
    };
    this->sc.set_atomCounts(atomCounts);
    std::map<int, int> orbitalCounts = {
        {0, 1},
        {1, 1}
    };
    this->sc.set_orbitalCounts(orbitalCounts);
    this->sc.Set_ScData_From_Json("./support/sc_f2.json");
    EXPECT_DOUBLE_EQ(this->sc.get_decay_grad(1), 0.9);
    this->sc.set_decay_grad_switch(true);
    this->sc.set_decay_grad();
    EXPECT_DOUBLE_EQ(this->sc.get_decay_grad().data()[1], 0.9 * 13.605698);
    int ntype = this->sc.get_ntype();
    std::vector<double> decay_grad = this->sc.get_decay_grad();
    this->sc.set_decay_grad(decay_grad.data(), ntype);
    EXPECT_DOUBLE_EQ(this->sc.get_decay_grad().data()[0], 0.0);
    EXPECT_DOUBLE_EQ(this->sc.get_decay_grad().data()[1], 0.9 * 13.605698);
    /// warning
    testing::internal::CaptureStdout();
    EXPECT_EXIT(this->sc.set_decay_grad(decay_grad.data(), 100), ::testing::ExitedWithCode(0), "");
    std::string output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("decay_grad_in size mismatch with ntype"));
}

TYPED_TEST(SpinConstrainTest, SetTargetMagType1)
{
    std::map<int, int> atomCounts = {
        {0, 1},
        {1, 5}
    };
    this->sc.set_atomCounts(atomCounts);
    std::map<int, int> orbitalCounts = {
        {0, 1},
        {1, 1}
    };
    this->sc.set_orbitalCounts(orbitalCounts);
    this->sc.Set_ScData_From_Json("./support/sc_f2.json");
    // set target mag type 1
    this->sc.set_target_mag();
    std::vector<ModuleBase::Vector3<double>> target_mag = this->sc.get_target_mag();
    for (const auto& sc_elem: this->sc.get_ScData())
    {
        int itype = sc_elem.first;
        const std::vector<ScAtomData>& sc_atoms = sc_elem.second;
        for (const ScAtomData& sc_data: sc_atoms)
        {
            int index = sc_data.index;
            int iat = this->sc.get_iat(itype, index);
            double mag_x = sc_data.target_mag_val * std::sin(sc_data.target_mag_angle1 * M_PI / 180)
                           * std::cos(sc_data.target_mag_angle2 * M_PI / 180);
            double mag_y = sc_data.target_mag_val * std::sin(sc_data.target_mag_angle1 * M_PI / 180)
                           * std::sin(sc_data.target_mag_angle2 * M_PI / 180);
            double mag_z = sc_data.target_mag_val * std::cos(sc_data.target_mag_angle1 * M_PI / 180);
            if (std::abs(mag_x) < 1e-14)
                mag_x = 0.0;
            if (std::abs(mag_y) < 1e-14)
                mag_y = 0.0;
            if (std::abs(mag_z) < 1e-14)
                mag_z = 0.0;
            EXPECT_DOUBLE_EQ(mag_x, target_mag[iat].x);
            EXPECT_DOUBLE_EQ(mag_y, target_mag[iat].y);
            EXPECT_DOUBLE_EQ(mag_z, target_mag[iat].z);
        }
    }
    for (int iat = 0; iat < this->sc.get_nat(); iat++)
    {
        if (iat == 1)
        {
            EXPECT_DOUBLE_EQ(target_mag[iat].x, 1.0);
            EXPECT_DOUBLE_EQ(target_mag[iat].y, 0.0);
            EXPECT_DOUBLE_EQ(target_mag[iat].z, 0.0);
        }
        else if (iat == 5)
        {
            EXPECT_DOUBLE_EQ(target_mag[iat].x, 0.0);
            EXPECT_DOUBLE_EQ(target_mag[iat].y, 1.5);
            EXPECT_DOUBLE_EQ(target_mag[iat].z, 0.0);
        }
        else
        {
            EXPECT_DOUBLE_EQ(target_mag[iat].x, 0.0);
            EXPECT_DOUBLE_EQ(target_mag[iat].y, 0.0);
            EXPECT_DOUBLE_EQ(target_mag[iat].z, 0.0);
        }
    }
}

TYPED_TEST(SpinConstrainTest, SetInputParameters)
{
    double sc_thr = 1e-6;
    int nsc = 100;
    int nsc_min = 2;
    double alpha_trial = 0.01;
    double sccut = 3.0;
    bool decay_grad_switch = 1;
    this->sc.set_input_parameters(sc_thr, nsc, nsc_min, alpha_trial, sccut, decay_grad_switch);
    EXPECT_DOUBLE_EQ(this->sc.get_sc_thr(), sc_thr);
    EXPECT_EQ(this->sc.get_nsc(), nsc);
    EXPECT_EQ(this->sc.get_nsc_min(), nsc_min);
    EXPECT_DOUBLE_EQ(this->sc.get_alpha_trial(), alpha_trial / 13.605698);
    EXPECT_DOUBLE_EQ(this->sc.get_sccut(), sccut / 13.605698);
    EXPECT_EQ(this->sc.get_decay_grad_switch(), decay_grad_switch);
}

TYPED_TEST(SpinConstrainTest, SetSolverParameters)
{
    K_Vectors kv;
    this->sc.set_solver_parameters(4, kv, nullptr, nullptr, nullptr, nullptr, "genelpa", nullptr);
    EXPECT_EQ(this->sc.get_nspin(), 4);
    EXPECT_EQ(this->sc.phsol, nullptr);
    EXPECT_EQ(this->sc.p_hamilt, nullptr);
    EXPECT_EQ(this->sc.psi, nullptr);
    EXPECT_EQ(this->sc.pelec, nullptr);
    EXPECT_EQ(this->sc.KS_SOLVER, "genelpa");
    EXPECT_EQ(this->sc.LM, nullptr);
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
    paraV.set_global2local(nrow, ncol, false, ofs);
    this->sc.set_ParaV(&paraV);
    EXPECT_EQ(this->sc.ParaV->nloc, nrow * ncol);
    remove("test.log");
}

TYPED_TEST(SpinConstrainTest, PrintMi)
{
    this->sc.zero_Mi();
    testing::internal::CaptureStdout();
    this->sc.print_Mi(true);
    std::string output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("Total Magnetism on atom: 0  (0, 0, 0)"));
    this->sc.set_nspin(2);
     testing::internal::CaptureStdout();
    this->sc.print_Mi(true);
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("Total Magnetism on atom: 0  (0)"));
}