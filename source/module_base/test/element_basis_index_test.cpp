#include "../element_basis_index.h"
#include "gtest/gtest.h"
#include "gmock/gmock.h"

/************************************************
 *  unit test of class Element_Basis_Index
 ***********************************************/

/**
 * - Tested Functions:
 *   - construct_index
 *     - make a index
 */

class IndexLNMTest : public testing::Test
{

};

TEST_F(IndexLNMTest,makeindex)
{
    ModuleBase::Element_Basis_Index::Range rangtest;
    rangtest.resize(2);
    rangtest[0].resize(3);
    rangtest[1].resize(3);
    rangtest[0][0].N=1;
    rangtest[0][0].M=2;
    rangtest[0][1].N=1;
    rangtest[0][1].M=3;
    rangtest[0][2].N=2;
    rangtest[0][2].M=1;
    rangtest[1][0].N=2;
    rangtest[1][0].M=2;
    rangtest[1][1].N=2;
    rangtest[1][1].M=3;
    rangtest[1][2].N=3;
    rangtest[1][2].M=3;
    ModuleBase::Element_Basis_Index::IndexLNM testindex;
    testindex=ModuleBase::Element_Basis_Index::construct_index(rangtest);
    EXPECT_EQ(rangtest[0][0].N,testindex[0][0].N);
    EXPECT_EQ(rangtest[0][0].M,testindex[0][0].M);
    EXPECT_EQ(rangtest[1][1].N,testindex[1][1].N);
    EXPECT_EQ(rangtest[1][1].M,testindex[1][1].M);
    EXPECT_EQ(rangtest[0][2].N,testindex[0][2].N);
    EXPECT_EQ(rangtest[0][2].M,testindex[0][2].M);
    EXPECT_EQ(testindex[0][0][0][0],0);
    EXPECT_EQ(testindex[0][0][0][1],1);
    EXPECT_EQ(testindex[0][1][0][0],2);
    EXPECT_EQ(testindex[1][1][0][0],4);
    EXPECT_EQ(testindex[1][2][0][1],11);    
    EXPECT_EQ(testindex[0].count_size,7);
    EXPECT_EQ(testindex[1].count_size,19);
}

