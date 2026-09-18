#ifndef JSON_NODE_H
#define JSON_NODE_H

#include <string>

namespace Json
{

class jsonKeyNode
{
  public:
    jsonKeyNode(int index) : i(index), is_index(true) {}
    jsonKeyNode(const std::string& name) : key(name) {}
    jsonKeyNode(const char* name) : key(name) {}

    int i = 0;
    std::string key;
    bool is_index = false;
};

} // namespace Json

#endif
