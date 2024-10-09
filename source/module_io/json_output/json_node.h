#ifndef JSON_NODE_H
#define JSON_NODE_H

namespace Json
{

    class jsonKeyNode{
        public:
            jsonKeyNode(int i): i(i) {};
            jsonKeyNode(const std::string& s): key(s) {};

            template<size_t N>
            jsonKeyNode(const char (&s)[N]): key(s) {};
        
            int i;
            std::string key;
    };

}

#endif