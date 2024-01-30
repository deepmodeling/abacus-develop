#include "abacusjson.h"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace Json
{

#ifdef __RAPIDJSON
rapidjson::Document AbacusJson::doc;

void AbacusJson::add_nested_member(std::vector<std::string>::iterator begin,
                                   std::vector<std::string>::iterator end,
                                   rapidjson::Value& val,
                                   rapidjson::Value& parent,
                                   rapidjson::Document::AllocatorType& allocator)
{
    if (begin != end)
    {
        rapidjson::Value key((*begin).c_str(), allocator);
        if (begin + 1 == end)
        {
            // if key exists, then overwrite it
            if (parent.HasMember(key))
            {
                // if key is an object, then warn the user
                if (parent[key].IsObject())
                {
                    std::cout << "Warning: write to json, key " << *begin
                              << " exist and is an object, and abacus will overwrite it with a value." << std::endl;
                }
                parent[key] = val;
            }
            else
                parent.AddMember(key, val, allocator);
        }
        else
        {
            // need to check if the key exists
            if (parent.HasMember(key))
            {
                // this key should be an object
                if (!parent[key].IsObject())
                {
                    std::cout << "Warning: write to json, key " << *begin
                              << " exist and is not an object, and abacus will add it as a middle node." << std::endl;
                }
                add_nested_member(begin + 1, end, val, parent[key], allocator);
            }
            else
            {
                rapidjson::Value val1(rapidjson::kObjectType);
                add_nested_member(begin + 1, end, val, val1, allocator);
                parent.AddMember(key, val1, allocator);
            }
        }
    }
}

// Output the json to a file
void AbacusJson::write_to_json(std::string filename)
{
    rapidjson::StringBuffer buffer;
    rapidjson::PrettyWriter<rapidjson::StringBuffer> writer(buffer);
    doc.Accept(writer);

    std::ofstream ofs(filename);
    ofs << buffer.GetString();
    ofs.close();
};
template <>
void AbacusJson::add_json(std::vector<std::string> keys, const std::string& value)
{
    if (!doc.IsObject())
    {
        doc.SetObject();
    }
    rapidjson::Value val(value.c_str(), doc.GetAllocator());
    add_nested_member(keys.begin(), keys.end(), val, doc, doc.GetAllocator());
}
#endif
} // namespace Json
