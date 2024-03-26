#ifndef ABACUS_JSON_H
#define ABACUS_JSON_H

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "json_node.h"


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



    // Output the json to a file
    static void write_to_json(std::string filename);

    static rapidjson::Document::AllocatorType& allocator(){
      return doc.GetAllocator();
    }

    /**
     *  @brief: The template specialization method adds value to the doc tree
     *
     *  @param: 'keys' is a vector string, represents the path to be added to the json tree.
     *
     *          'value' is the value that needs to be added to the json tree, which type can be 
     *          Json::jsonValue or other common value type(such as int, double ,bool ,std::string).
     *
     *          'IsArray' is a bool value, means the whether the root node to which 'value' is added is an array.
     *
     *  @usage: 1. Add/Modify a double val to object json node (key2 is a object node): 
     *                Json::AbacusJson::add_json({"key1","key2"}, 3.1415,false);
     *
     *          2. Pushback a double val to array json node (key2 is a array node):
     *                Json::AbacusJson::add_json({"key1","key2"}, 3.1415,true);
     *
     *          3. Modify a doble val to array json node (key2 is a array node), when use this method, 
     *            The index number of the array starts at 0, if it's negative, it's going from back to front. 
     *            eg. If the index is -1, it means that the last element of the array is modified:
     *            If we have a json array: {"key":[1,2,3]}
     *                i). Json::AbacusJson::add_json({"key",0}, 4,true);  => {"key":[4,2,3]}
     *                ii). Json::AbacusJson::add_json({"key",-1}, 4,true);  => {"key":[1,2,4]}
     *                iii). Json::AbacusJson::add_json({"key",1}, 4,true);  => {"key":[1,4,3]}
     *                iv). Json::AbacusJson::add_json({"key",2}, 4,true);  => {"key":[1,2,4]}
     *                iv). Json::AbacusJson::add_json({"key",3}, 4,true);  => error!, The array element corresponding
     *                     to the index has no value.
    */
    template <typename T>
    static void add_json(std::vector<jsonKeyNode> keys, const T& value,bool IsArray)
    {
        if (!doc.IsObject())
        {
            doc.SetObject();
        }
        rapidjson::Value val(value);
        add_nested_member(keys.begin(), keys.end(), val, doc, doc.GetAllocator(),IsArray);
    }




  private:
    static rapidjson::Document doc;

    static void add_nested_member(std::vector<jsonKeyNode>::iterator begin,
                                    std::vector<jsonKeyNode>::iterator end,
                                    rapidjson::Value& val,
                                    rapidjson::Value& parent,
                                    rapidjson::Document::AllocatorType& allocator,
                                    bool IsArray
                                    );
  };

template <>
void AbacusJson::add_json(std::vector<jsonKeyNode> keys, const std::string& value,bool IsArray);

template <>
void AbacusJson::add_json(std::vector<jsonKeyNode> keys, const rapidjson::Value& value,bool IsArray);


} // namespace Json
#endif

#endif