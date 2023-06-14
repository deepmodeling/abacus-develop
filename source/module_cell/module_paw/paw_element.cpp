#include "paw_element.h"
#include "module_base/tool_title.h"
#include "module_base/tool_quit.h"

void Paw_Element::read_paw_xml(std::string filename)
{
    ModuleBase::TITLE("Paw_Element","read_paw_xml");

    std::string line;
    std::string s1 = "<atom";

    std::ifstream ifs(filename.c_str());
    if(ifs.fail())
    {
        ModuleBase::WARNING_QUIT("paw_element.cpp","xml file not found!");
    }

    // Zat, core and valence electrons
    line = this->scan_file(ifs, s1);

}

std::string Paw_Element::scan_file(std::ifstream &ifs, std::string pattern)
{
    std::string line;

    while (!ifs.eof())
    {
        getline(ifs,line);
        if (line.find(pattern) != std::string::npos) return line;
    }

    ModuleBase::WARNING_QUIT("paw_element.cpp","pattern not found in xml file!");
    return 0;    
}