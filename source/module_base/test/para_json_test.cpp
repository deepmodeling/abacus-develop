#include "module_base/para_json.h"
#include "gtest/gtest.h"
#ifdef __MPI
#include "mpi.h"
#endif

#include <stdlib.h>
#include "rapidjson/document.h"
/************************************************
 *  unit test of Input::ParaJson
 ***********************************************/

/**
 * - Tested Functions:
 *   - Init()
 *     - init json tree from input::Init and check if the json string is valid
 */

class ParaJsonTest : public ::testing::Test
{
  protected:
    std::string testString;
};

// check if a string is a valid JSON string
bool isValidJSON(const std::string& jsonString)
{
    rapidjson::Document document;
    document.Parse(jsonString.c_str());

    return !document.HasParseError();
}

TEST_F(ParaJsonTest, Init)
{
    //std::string input_file = "./support/INPUT";
    //Input input_tmp;
    //EXPECT_NO_THROW(input_tmp.Init(input_file));

    // int status = system("rm -r ./OUT.autotest/");
    // EXPECT_EQ(status,0);
    // Para_Json::Init_json_abacus_readinInfo();
    Para_Json::Init_json_abacus_generalInfo();
    Para_Json::Init_json_abacus();
    Para_Json::Finish_json_tree();
    rapidjson::StringBuffer buffer;
    rapidjson::Writer<rapidjson::StringBuffer> writer(buffer);
    Para_Json::doc.Accept(writer);
    std::string json = buffer.GetString();
    EXPECT_EQ(isValidJSON(json), true);
}

int main(int argc, char** argv)
{
#ifdef __MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &GlobalV::NPROC);
    MPI_Comm_rank(MPI_COMM_WORLD, &GlobalV::MY_RANK);
#endif
    testing::InitGoogleTest(&argc, argv);
    int result;
    result = RUN_ALL_TESTS();
#ifdef __MPI
    MPI_Finalize();
#endif
    return result;
}

