#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "../atom_in.h"

/************************************************
 *  unit test of functions in atom_in.h
 ***********************************************/

/**
 * - Tested functions of class atom_in:
 *   - map<string, int> atom_Z
 *     - get the atomic number
 * 
 *   - map<string, double> atom_RCS
 *     - get the atomic radius
 * 
 */

class atom_in_test : public testing::Test
{
protected:

    atom_in Atomin;
    int a = Atomin.atom_Z["H"];
    
};

TEST_F(atom_in_test, atomin)
{   
    EXPECT_EQ(atom_in_test::a, 1);
    EXPECT_EQ(Atomin.atom_RCS["H"],0.603774);
}