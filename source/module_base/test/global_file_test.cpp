#include "../global_file.h"
#include "../global_variable.h"
#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include <fstream>
#include <unistd.h>

#ifdef __MPI
#include "mpi.h"
#endif

/************************************************
 *  unit test of functions in global_file.cpp
 ***********************************************/

/**
 * - Tested Function
 * - mkdiratom
 *   - generate atom dir for each type of atom
 */

class GlobalFile : public testing::Test
{

};

TEST_F(GlobalFile,mkdiratom)
{
        GlobalV::global_out_dir = "./";
        ModuleBase::Global_File::make_dir_atom("Si");
        int a = access("./Si/",0);
        EXPECT_EQ(a , 0);
}
