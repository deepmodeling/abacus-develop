#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#ifdef __RAPIDJSON
#include "rapidjson/document.h"
#include "rapidjson/prettywriter.h"
#include "rapidjson/stringbuffer.h"

//define of AbacusJson
#define JobjectType rapidjson::kObjectType
#define JarrayType rapidjson::kArrayType
#define Get_Jallocator Json::AbacusJson::allocator()
#define Set_JString(str) Json::jsonValue().SetString(str.c_str(),str.length(),Get_Jallocator)
#define JaddString(str,val) AddMember(str, Set_JString(val), Get_Jallocator)
#define JaddNormal(str,val) AddMember(str, val, Get_Jallocator)
#define JPushBack(val) PushBack(val, Get_Jallocator)


namespace Json
{

using jsonValue = rapidjson::Value;
// This class is used to construct a json value, and output some key values to a json file.
class AbacusJson
{
  public:
    // Add a key value pair to the json
    template <typename T>
    static void add_json(std::vector<std::string> keys, const T& value,bool IsArray)
    {
        if (!doc.IsObject())
        {
            doc.SetObject();
        }
        rapidjson::Value val(value);
        add_nested_member(keys.begin(), keys.end(), val, doc, doc.GetAllocator(),IsArray);
    }

    // Output the json to a file
    static void write_to_json(std::string filename);

    static rapidjson::Document::AllocatorType& allocator(){
      return doc.GetAllocator();
    }



  private:
    static rapidjson::Document doc;

    static void add_nested_member(std::vector<std::string>::iterator begin,
                                  std::vector<std::string>::iterator end,
                                  rapidjson::Value& val,
                                  rapidjson::Value& parent,
                                  rapidjson::Document::AllocatorType& allocator,
                                  bool IsArray = false
                                  );
};
template <>
void AbacusJson::add_json(std::vector<std::string> keys, const std::string& value,bool IsArray);

template <>
void AbacusJson::add_json(std::vector<std::string> keys, const rapidjson::Value& value,bool IsArray);


} // namespace Json
#endif