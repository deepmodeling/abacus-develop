#include "gtest/gtest.h"
#define private public
#define __RAPIDJSON 1
#include "../abacusjson.h"
#include "module_io/input.h"
#include "module_io/para_json.h"
#include "version.h"
TEST(AbacusJsonTest, AddJson)
{
    Json::AbacusJson::doc.SetObject();

    // add a string
    Json::AbacusJson::add_json({"key1"}, "value1",false);
    ASSERT_TRUE(Json::AbacusJson::doc.HasMember("key1"));
    ASSERT_TRUE(Json::AbacusJson::doc["key1"].IsString());
    ASSERT_STREQ(Json::AbacusJson::doc["key1"].GetString(), "value1");

    // add a string to a nested object
    Json::AbacusJson::add_json({"key2", "key3"}, "value2",false);
    ASSERT_TRUE(Json::AbacusJson::doc.HasMember("key2"));
    ASSERT_TRUE(Json::AbacusJson::doc["key2"].IsObject());
    ASSERT_TRUE(Json::AbacusJson::doc["key2"].HasMember("key3"));
    ASSERT_TRUE(Json::AbacusJson::doc["key2"]["key3"].IsString());
    ASSERT_STREQ(Json::AbacusJson::doc["key2"]["key3"].GetString(), "value2");

    // add an int
    Json::AbacusJson::add_json({"key2"}, 123,false);
    ASSERT_TRUE(Json::AbacusJson::doc.HasMember("key2"));
    ASSERT_TRUE(Json::AbacusJson::doc["key2"].IsInt());
    ASSERT_EQ(Json::AbacusJson::doc["key2"].GetInt(), 123);

    // add a bool
    Json::AbacusJson::add_json({"key3"}, true,false);
    ASSERT_TRUE(Json::AbacusJson::doc.HasMember("key3"));
    ASSERT_TRUE(Json::AbacusJson::doc["key3"].IsBool());
    ASSERT_EQ(Json::AbacusJson::doc["key3"].GetBool(), true);

    // add a double
    Json::AbacusJson::add_json({"key4"}, 1.23,false);
    ASSERT_TRUE(Json::AbacusJson::doc.HasMember("key4"));
    ASSERT_TRUE(Json::AbacusJson::doc["key4"].IsDouble());
    ASSERT_EQ(Json::AbacusJson::doc["key4"].GetDouble(), 1.23);

    // modify a value
    Json::AbacusJson::add_json({"key4"}, 4.56,false);
    ASSERT_EQ(Json::AbacusJson::doc["key4"].GetDouble(), 4.56);

    // modify a value
    Json::AbacusJson::add_json({"key4"}, 4.56,false);    


    //array test
    Json::AbacusJson::add_json({"key6","key7"}, true,true);
    Json::AbacusJson::add_json({"key6","key7"}, false,true);

    
    // add key-val to a object array
    for(int i=0;i<3;i++){
        Json::jsonValue object(JobjectType);
        object.JaddNormal("int",i);
        std::string str = std::to_string(i*100);   
        object.JaddString("string", str);
        object.JaddNormal("double", 0.01*i);    
        Json::AbacusJson::add_json({"array"}, object,true);
    }
    ASSERT_EQ(Json::AbacusJson::doc["array"][0]["int"].GetInt(), 0);
    ASSERT_STREQ(Json::AbacusJson::doc["array"][0]["string"].GetString(), "0");
    ASSERT_EQ(Json::AbacusJson::doc["array"][0]["double"].GetDouble(), 0.0);

    ASSERT_EQ(Json::AbacusJson::doc["array"][1]["int"].GetInt(), 1);
    ASSERT_STREQ(Json::AbacusJson::doc["array"][1]["string"].GetString(), "100");
    ASSERT_EQ(Json::AbacusJson::doc["array"][1]["double"].GetDouble(), 0.01);

    ASSERT_EQ(Json::AbacusJson::doc["array"][2]["int"].GetInt(), 2);
    ASSERT_STREQ(Json::AbacusJson::doc["array"][2]["string"].GetString(), "200");
    ASSERT_EQ(Json::AbacusJson::doc["array"][2]["double"].GetDouble(), 0.02);

    // add array in array
    Json::jsonValue object0(JarrayType);

    object0.JPushBack(1);
    object0.JPushBack(2);
    object0.JPushBack(3);
    
    Json::jsonValue object1(JarrayType);

    object1.JPushBack(2);
    object1.JPushBack(3);
    object1.JPushBack(4);

    Json::AbacusJson::add_json({"Darray"}, object0,true);
    Json::AbacusJson::add_json({"Darray"}, object1,true);

    ASSERT_EQ(Json::AbacusJson::doc["Darray"][0][0].GetInt(), 1);
    ASSERT_EQ(Json::AbacusJson::doc["Darray"][0][1].GetInt(), 2);
    ASSERT_EQ(Json::AbacusJson::doc["Darray"][0][2].GetInt(), 3);

    ASSERT_EQ(Json::AbacusJson::doc["Darray"][1][0].GetInt(), 2);
    ASSERT_EQ(Json::AbacusJson::doc["Darray"][1][1].GetInt(), 3);
    ASSERT_EQ(Json::AbacusJson::doc["Darray"][1][2].GetInt(), 4);
}

TEST(AbacusJsonTest, OutputJson)
{
    Json::AbacusJson::doc.SetObject();


    Json::AbacusJson::add_json({"key1"}, "value1",false);
    Json::AbacusJson::add_json({"key2", "key3"}, 1,false);
    Json::AbacusJson::add_json({"key4"}, 0.1,false);
    Json::AbacusJson::add_json({"key5"}, true,false);

    //array test
    Json::AbacusJson::add_json({"key6","key7"}, true,true);
    Json::AbacusJson::add_json({"key6","key7"}, false,true);


 
    std::string filename = "test.json";
    Json::AbacusJson::write_to_json(filename);

    std::ifstream file(filename);
    ASSERT_TRUE(file.is_open());

    std::string content((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
    ASSERT_NE(content.find("\"key1\": \"value1\","), std::string::npos);
    ASSERT_NE(content.find("\"key2\": {"), std::string::npos);
    ASSERT_NE(content.find("\"key3\": 1"), std::string::npos);
    ASSERT_NE(content.find("\"key4\": 0.1"), std::string::npos);
    ASSERT_NE(content.find("\"key5\": true"), std::string::npos);

    file.close();
}

TEST(AbacusJsonTest, GeneralInfo)
{
    std::time_t time_now = std::time(NULL);
    std::string start_time_str;
    Json::convert_time(time_now, start_time_str);
    INPUT.start_time = time_now;
    
    INPUT.device = "cpu";
    INPUT.pseudo_dir = "./abacus/test/pseudo_dir";
    INPUT.orbital_dir = "./abacus/test/orbital_dir";
    INPUT.stru_file = "./abacus/test/stru_file";
    INPUT.kpoint_file = "./abacus/test/kpoint_file";
    // output the json file
    Json::AbacusJson::doc.Parse("{}");
    Json::json_output();

    std::string filename = "abacus.json";
    std::ifstream file(filename);
    ASSERT_TRUE(file.is_open());

    std::string content((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());

    ASSERT_NE(content.find(VERSION), std::string::npos);
    ASSERT_NE(content.find("\"device\": \"cpu\","), std::string::npos);
    ASSERT_NE(content.find("\"omp_num\": 0,"), std::string::npos);
    ASSERT_NE(content.find("\"mpi_num\": 0,"), std::string::npos);
    ASSERT_NE(content.find("\"orbital_dir\": \"./abacus/test/orbital_dir\","), std::string::npos);
    ASSERT_NE(content.find("\"pseudo_dir\": \"./abacus/test/pseudo_dir\","), std::string::npos);
    ASSERT_NE(content.find("\"stru_file\": \"./abacus/test/stru_file\","), std::string::npos);
    ASSERT_NE(content.find("\"kpt_file\": \"./abacus/test/kpoint_file\","), std::string::npos);
    ASSERT_NE(content.find(start_time_str), std::string::npos);
}
