#ifndef ABACUS_JSON_H
#define ABACUS_JSON_H

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#ifdef __RAPIDJSON
#include "rapidjson/document.h"
#include "rapidjson/prettywriter.h"
#include "rapidjson/stringbuffer.h"

/**
* @brief Define of AbacusJson:These macro definitions simplify the complex parameters 
*          required to add objects in rapidjson, making it easy to use.
* @usage: 1. create the json value by Type: eg.  
*            Json::jsonValue object(JobjectType); 
*            Json::jsonValue array(JarrayType);
*         2 - object. add the parameter by using correct function. If the constructed object contains 
*            type std::string, select the correct macro definition function in the key 
*            or value position based on std::string.
*                 eg. key is std::string, using: object.JaddStringK(str,val)
*                 eg. val is std::string, using: object.JaddStringV (str,val)
*                 eg. both key,val is std::string, using: object.JaddStringKV(str,val)
*                 eg. none of key,val is std::string, using: object.JaddNormal(str,val)
*
*         2 - array.  using: array.JPushBack(val)
*                     or : JPushBackString(val)
*/

#define JobjectType rapidjson::kObjectType
#define JarrayType rapidjson::kArrayType
#define Get_Jallocator Json::AbacusJson::allocator()
#define Set_JString(str) Json::jsonValue().SetString(str.c_str(),str.length(),Get_Jallocator)


#define JaddStringV(str,val) AddMember(str, Set_JString(val), Get_Jallocator)
#define JaddStringK(str,val) AddMember(Set_JString(str), val, Get_Jallocator)
#define JaddStringKV(str,val) AddMember(Set_JString(str), Set_JString(val), Get_Jallocator)
#define JaddNormal(str,val) AddMember(str, val, Get_Jallocator)
#define JPushBack(val) PushBack(val, Get_Jallocator)
#define JPushBackString(val) PushBack(Set_JString(val), Get_Jallocator)
/**
* @brief The json processing module in abacus contains basic meta-operations such as 
*         adding, modifying, and checking json.
*/

namespace Json
{

using jsonValue = rapidjson::Value;
// This class is used to construct a json value, and output some key values to a json file.
class AbacusJson
{
  public:
    // The template specialization method adds a key-val object to the doc tree
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

#endif