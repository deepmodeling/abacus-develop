#include "../tool_title.h"
#include "../global_variable.h"
#include "gtest/gtest.h"
#include "gmock/gmock.h"

/************************************************
 *  unit test of functions in tool_title.h
 ***********************************************/

/**
 * - Tested Function
 *   - ModuleBase::TITLE
 *    - Output title for each function.
 */

class ToolTitleTest : public testing::Test
{
  protected:
	std::ifstream ifs;
    const std::string claname="ComplexMatrix";
    const std::string funname="scaled_sum()";
    const std::string cfname="ComplexMatrix::scaled_sum()";
  void SetUp()
        {
        
        }
  void TearDown()
	    {
          remove("TITLEtest2.log");
          remove("TITLEtest3.log");
	    }
};

TEST_F(ToolTitleTest, TITLE2)
{
    GlobalV::ofs_running.open("TITLEtest2.log");
    std::string output2;
    ModuleBase::TITLE(claname,funname,false);
    GlobalV::ofs_running.close();
    ifs.open("TITLEtest2.log");
    getline(ifs,output2);
    EXPECT_THAT(output2,testing::HasSubstr(" ==> ComplexMatrix::scaled_sum()"));
    ifs.close();
}

TEST_F(ToolTitleTest, TITLE3)
{
    std::ofstream oofs;  
    std::string output3a;
    std::string output3b;
    oofs.open("TITLEtest3.log");
    ModuleBase::TITLE(oofs,claname,funname,false);
    oofs.close();
    ifs.open("TITLEtest3.log");
    getline(ifs,output3a);
    EXPECT_THAT(output3a,testing::HasSubstr(" ==> ComplexMatrix::scaled_sum()"));
    ifs.close();
}
