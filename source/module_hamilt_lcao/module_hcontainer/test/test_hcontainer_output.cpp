#include "../output_hcontainer.h"

#include "../hcontainer.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "module_cell/unitcell.h"

/************************************************
 *  unit test of output_hcontainer.cpp
 ***********************************************/

/**
 * - Tested Functions:
 * - write()
 * - write the matrices of all R vectors to the output stream
 * - write(int rx_in, int ry_in, int rz_in)
 * - write the matrix of a single R vector to the output stream
 * - write_single_R(int rx, int ry, int rz)
 * - write the matrix of a single R vector to the output stream
 */

class OutputHContainerTest : public testing::Test
{
  protected:
    UnitCell ucell;
    std::string output;
    void SetUp() override
    {
        ucell.ntype = 1;
        ucell.nat = 2;
        ucell.atoms = new Atom[ucell.ntype];
        ucell.iat2it = new int[ucell.nat];
        ucell.iat2ia = new int[ucell.nat];
        for (int iat = 0; iat < ucell.nat; iat++)
        {
            ucell.iat2ia[iat] = iat;
            ucell.iat2it[iat] = 0;
        }
        ucell.atoms[0].na = 2;
        ucell.atoms[0].nw = 2;
        ucell.iwt2iat = new int[4];
        ucell.iwt2iw = new int[4];
        ucell.itia2iat.create(ucell.ntype, ucell.nat);
        ucell.set_iat2iwt(1);
        ucell.itia2iat(0, 0) = 0;
        ucell.itia2iat(0, 1) = 1;
        ucell.iwt2iat[0] = 0;
        ucell.iwt2iat[1] = 0;
        ucell.iwt2iat[2] = 1;
        ucell.iwt2iat[3] = 1;
        ucell.iwt2iw[0] = 0;
        ucell.iwt2iw[1] = 1;
        ucell.iwt2iw[2] = 0;
        ucell.iwt2iw[3] = 1;
    }
    void TearDown() override
    {
        delete[] ucell.atoms;
        delete[] ucell.iat2it;
        delete[] ucell.iat2ia;
        delete[] ucell.iwt2iat;
        delete[] ucell.iwt2iw;
    }
};

TEST_F(OutputHContainerTest, Write)
{
    Parallel_Orbitals ParaV;
    ParaV.atom_begin_row.resize(3);
    ParaV.atom_begin_col.resize(3);
    for (int i = 0; i < 3; i++)
    {
        ParaV.atom_begin_row[i] = i * 2;
        ParaV.atom_begin_col[i] = i * 2;
    }
    ParaV.nrow = 4;
    ParaV.ncol = 4;
    std::ofstream ofs("output_hcontainer.log");
    ParaV.set_global2local(4, 4, false, ofs);
    // std::cout << "ParaV.global2local_row = " << ParaV.global2local_row(0) << " " << ParaV.global2local_row(1) << " "
    //           << ParaV.global2local_row(2) << " " << ParaV.global2local_row(3) << std::endl;
    // std::cout << "ParaV.global2local_col = " << ParaV.global2local_col(0) << " " << ParaV.global2local_col(1) << " "
    //           << ParaV.global2local_col(2) << " " << ParaV.global2local_col(3) << std::endl;
    ofs.close();
    remove("output_hcontainer.log");
    hamilt::HContainer<double> HR(&ParaV);
    double correct_array[16] = {1, 2, 0, 4, 5, 0, 7, 0, 3, 0, 5, 6, 7, 8, 0, 10};
    double correct_array1[16] = {1, 2, 0, 4, 5, 0, 7, 0, 3, 0, 5, 6, 7, 8, 0, 10};
    // correct_array represent a matrix of
    // 1 2 0 4
    // 5 0 7 0
    // 3 0 5 6
    // 7 8 0 10
    double test_data[8] = {0, 4, 7, 0, 5, 6, 0, 10};
    hamilt::AtomPair<double> ap1(0, 1, 0, 1, 1, &ParaV, &test_data[0]);
    hamilt::AtomPair<double> ap2(1, 1, 0, 0, 0, &ParaV, &test_data[4]);
    HR.insert_pair(ap1);
    HR.insert_pair(ap2);
    for (int ir = 0; ir < HR.size_R_loop(); ++ir)
    {
        int rx, ry, rz;
        HR.loop_R(ir, rx, ry, rz);
        HR.fix_R(rx, ry, rz);
        // std::cout << "rx = " << rx << " ry = " << ry << " rz = " << rz << std::endl;
        for (int iap = 0; iap < HR.size_atom_pairs(); ++iap)
        {
            hamilt::AtomPair<double>& tmp_ap = HR.get_atom_pair(iap);
            if (rx == 0 && ry == 1 && rz == 1)
            {
                EXPECT_DOUBLE_EQ(tmp_ap.get_value(0, 0), 0);
                EXPECT_DOUBLE_EQ(tmp_ap.get_value(0, 1), 4);
                EXPECT_DOUBLE_EQ(tmp_ap.get_value(1, 0), 7);
                EXPECT_DOUBLE_EQ(tmp_ap.get_value(1, 1), 0);
                EXPECT_DOUBLE_EQ(std::get<0>(tmp_ap.get_matrix_values())[0], 0);
                EXPECT_DOUBLE_EQ(std::get<0>(tmp_ap.get_matrix_values())[1], 2);
                EXPECT_DOUBLE_EQ(std::get<0>(tmp_ap.get_matrix_values())[2], 2);
                EXPECT_DOUBLE_EQ(std::get<0>(tmp_ap.get_matrix_values())[3], 2);
            }
            else if (rx == 0 && ry == 0 && rz == 0)
            {
                EXPECT_DOUBLE_EQ(tmp_ap.get_value(0, 0), 5);
                EXPECT_DOUBLE_EQ(tmp_ap.get_value(0, 1), 6);
                EXPECT_DOUBLE_EQ(tmp_ap.get_value(1, 0), 0);
                EXPECT_DOUBLE_EQ(tmp_ap.get_value(1, 1), 10);
                EXPECT_DOUBLE_EQ(std::get<0>(tmp_ap.get_matrix_values())[0], 2);
                EXPECT_DOUBLE_EQ(std::get<0>(tmp_ap.get_matrix_values())[1], 2);
                EXPECT_DOUBLE_EQ(std::get<0>(tmp_ap.get_matrix_values())[2], 2);
                EXPECT_DOUBLE_EQ(std::get<0>(tmp_ap.get_matrix_values())[3], 2);
            }
        }
        HR.unfix_R();
    }
    double sparse_threshold = 0.1;
    hamilt::Output_HContainer<double> output_HR(&HR, &ParaV, ucell, std::cout, sparse_threshold, 2);
    // the first R
    testing::internal::CaptureStdout();
    output_HR.write(0, 1, 1);
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("0 1 1 2"));
    EXPECT_THAT(output, testing::HasSubstr(" 4.00e+00 7.00e+00"));
    EXPECT_THAT(output, testing::HasSubstr(" 3 2"));
    EXPECT_THAT(output, testing::HasSubstr(" 0 1 2 2 2"));
    // the second R
    testing::internal::CaptureStdout();
    output_HR.write(0, 0, 0);
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("0 0 0 3"));
    EXPECT_THAT(output, testing::HasSubstr(" 5.00e+00 6.00e+00 1.00e+01"));
    EXPECT_THAT(output, testing::HasSubstr(" 2 3 3"));
    EXPECT_THAT(output, testing::HasSubstr(" 0 0 0 2 3"));
    // output all R
    testing::internal::CaptureStdout();
    output_HR.write();
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("0 1 1 2"));
    EXPECT_THAT(output, testing::HasSubstr(" 4.00e+00 7.00e+00"));
    EXPECT_THAT(output, testing::HasSubstr(" 3 2"));
    EXPECT_THAT(output, testing::HasSubstr(" 0 1 2 2 2"));
    EXPECT_THAT(output, testing::HasSubstr("0 0 0 3"));
    EXPECT_THAT(output, testing::HasSubstr(" 5.00e+00 6.00e+00 1.00e+01"));
    EXPECT_THAT(output, testing::HasSubstr(" 2 3 3"));
    EXPECT_THAT(output, testing::HasSubstr(" 0 0 0 2 3"));
}