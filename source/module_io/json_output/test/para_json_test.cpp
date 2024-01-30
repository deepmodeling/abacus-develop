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
    Json::AbacusJson::add_json({"key1"}, "value1");
    ASSERT_TRUE(Json::AbacusJson::doc.HasMember("key1"));
    ASSERT_TRUE(Json::AbacusJson::doc["key1"].IsString());
    ASSERT_STREQ(Json::AbacusJson::doc["key1"].GetString(), "value1");

    // add a string to a nested object
    Json::AbacusJson::add_json({"key2", "key3"}, "value2");
    ASSERT_TRUE(Json::AbacusJson::doc.HasMember("key2"));
    ASSERT_TRUE(Json::AbacusJson::doc["key2"].IsObject());
    ASSERT_TRUE(Json::AbacusJson::doc["key2"].HasMember("key3"));
    ASSERT_TRUE(Json::AbacusJson::doc["key2"]["key3"].IsString());
    ASSERT_STREQ(Json::AbacusJson::doc["key2"]["key3"].GetString(), "value2");

    // add an int
    Json::AbacusJson::add_json({"key2"}, 123);
    ASSERT_TRUE(Json::AbacusJson::doc.HasMember("key2"));
    ASSERT_TRUE(Json::AbacusJson::doc["key2"].IsInt());
    ASSERT_EQ(Json::AbacusJson::doc["key2"].GetInt(), 123);

    // add a bool
    Json::AbacusJson::add_json({"key3"}, true);
    ASSERT_TRUE(Json::AbacusJson::doc.HasMember("key3"));
    ASSERT_TRUE(Json::AbacusJson::doc["key3"].IsBool());
    ASSERT_EQ(Json::AbacusJson::doc["key3"].GetBool(), true);

    // add a double
    Json::AbacusJson::add_json({"key4"}, 1.23);
    ASSERT_TRUE(Json::AbacusJson::doc.HasMember("key4"));
    ASSERT_TRUE(Json::AbacusJson::doc["key4"].IsDouble());
    ASSERT_EQ(Json::AbacusJson::doc["key4"].GetDouble(), 1.23);

    // modify a value
    Json::AbacusJson::add_json({"key4"}, 4.56);
    ASSERT_EQ(Json::AbacusJson::doc["key4"].GetDouble(), 4.56);
}

TEST(AbacusJsonTest, OutputJson)
{
    Json::AbacusJson::doc.SetObject();

    Json::AbacusJson::add_json({"key1"}, "value1");
    Json::AbacusJson::add_json({"key2", "key3"}, 1);
    Json::AbacusJson::add_json({"key4"}, 0.1);
    Json::AbacusJson::add_json({"key5"}, true);

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
