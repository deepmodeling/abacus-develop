/*
To print any kind of information to an XML file
*/
#include <string>
#include <vector>
#include <memory>
/* this is the class to print information in xml format, NOT THE ONE CARRYING DATA */
class to_xml {
    public:
        to_xml();
        ~to_xml();

    private:
        std::shared_ptr<xml> root_;
};

/* this is the class to carry data, relationship between sections are record in std::vector<std::shared_ptr<>> container */
class xml {
    public:
        xml();
        ~xml();
    
    private:
        std::string name_;
        std::vector<std::string> data_;
        std::vector<std::pair<std::string, std::string>> attributes_;

        std::vector<std::shared_ptr<xml>> children_;
};