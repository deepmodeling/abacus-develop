#include "module_basis/module_ao/parallel_orbitals.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

/***********************************************************
 *      unit test of class "Parallel_Orbitals"
 ***********************************************************/

/**
 * Tested functions:
 */

class Test_Parallel_Orbitals : public Parallel_Orbitals
{
  public:
    Test_Parallel_Orbitals() = default;
    ~Test_Parallel_Orbitals() = default;
    void set_serial(int& M_A, const int& N_A)
    {
        this->nrow = M_A;
        this->ncol = N_A;
        this->nloc = this->nrow * this->ncol;
        this->local2global_row_.resize(this->nrow);
        this->local2global_col_.resize(this->ncol);
        for (int i = 0; i < this->nrow; i++)
            this->local2global_row_[i] = i;
        for (int i = 0; i < this->ncol; i++)
            this->local2global_col_[i] = i;
    }
};

class TestParallelOrbitals : public testing::Test
{
  protected:
};

TEST_F(TestParallelOrbitals, SetAtomicTrace8x8Atom2)
{
    Test_Parallel_Orbitals paraV;
    int nloc = 8;
    // set up global2local_row_ and global2local_col_
    paraV.set_serial(nloc, nloc);
    std::ofstream ofs("test.txt");
    paraV.set_global2local(nloc, nloc, false, ofs);
    ofs.close();
    remove("test.txt");
    int nat = 2;
    // set up iat2iwt
    int iat2iwt[] = {0, 4};
    // call set_atomic_trace
    paraV.set_atomic_trace(iat2iwt, nat, nloc);
    // expected results
    std::vector<int> atom_begin_col = {0, 4};
    std::vector<int> atom_begin_row = {0, 4};
    // check results
    EXPECT_EQ(paraV.atom_begin_col, atom_begin_col);
    EXPECT_EQ(paraV.atom_begin_row, atom_begin_row);
    // check getters
    EXPECT_EQ(paraV.get_col_size(), 8);
    EXPECT_EQ(paraV.get_row_size(), 8);
    EXPECT_EQ(paraV.get_col_size(0), 4);
    EXPECT_EQ(paraV.get_row_size(0), 4);
    EXPECT_EQ(paraV.get_col_size(1), 4);
    EXPECT_EQ(paraV.get_row_size(1), 4);
}

TEST_F(TestParallelOrbitals, SetAtomicTrace10x10Atom3)
{
    Test_Parallel_Orbitals paraV;
    int nloc = 10;
    // set up global2local_row_ and global2local_col_
    paraV.set_serial(nloc, nloc);
    std::ofstream ofs("test.txt");
    paraV.set_global2local(nloc, nloc, false, ofs);
    ofs.close();
    remove("test.txt");
    int nat = 3;
    // set up iat2iwt
    int iat2iwt[] = {0, 2, 5};
    // call set_atomic_trace
    paraV.set_atomic_trace(iat2iwt, nat, nloc);
    // expected results
    std::vector<int> atom_begin_col = {0, 2, 5};
    std::vector<int> atom_begin_row = {0, 2, 5};
    // check results
    EXPECT_EQ(paraV.atom_begin_col, atom_begin_col);
    EXPECT_EQ(paraV.atom_begin_row, atom_begin_row);
    // check getters
    EXPECT_EQ(paraV.get_col_size(), 10);
    EXPECT_EQ(paraV.get_row_size(), 10);
    EXPECT_EQ(paraV.get_col_size(0), 2);
    EXPECT_EQ(paraV.get_row_size(0), 2);
    EXPECT_EQ(paraV.get_col_size(1), 3);
    EXPECT_EQ(paraV.get_row_size(1), 3);
    EXPECT_EQ(paraV.get_col_size(2), 5);
    EXPECT_EQ(paraV.get_row_size(2), 5);
}