#include "read_input.h"

#include "module_base/tool_quit.h"
#include "module_base/tool_title.h"

#include <fstream>
#include <string.h>
namespace ModuleIO
{

void strtolower(char* sa, char* sb)
{
    char c;
    int len = strlen(sa);
    for (int i = 0; i < len; i++)
    {
        c = sa[i];
        sb[i] = tolower(c);
    }
    sb[len] = '\0';
}

void read_information(std::ifstream& ifs, std::vector<std::string>& output, const std::string& delimiters)
{
    std::string value1;
    ifs >> value1;
    output.push_back(value1);

    std::string line;
    getline(ifs, line);

    std::istringstream iss(line);
    std::string word;
    while (iss >> word)
    {
        if (delimiters.find(word[0]) != std::string::npos)
            break;
        output.push_back(word);
    }
}

ReadInput::ReadInput()
{
    this->item_general();
}

void ReadInput::ReadTxtInput(Parameter& param, const std::string& filename)
{
    ModuleBase::TITLE("Input", "Read");

    if (this->rank != 0)
        return ;

    std::ifstream ifs(filename.c_str(), std::ios::in);

    if (!ifs)
    {
        std::cout << " Can't find the INPUT file." << std::endl;
        ModuleBase::WARNING_QUIT("Input::Init", "Error during readin parameters.", 1);
    }

    ifs.clear();
    ifs.seekg(0);

    char word[80];
    char word1[80];
    int ierr = 0;

    // ifs >> std::setiosflags(ios::uppercase);
    ifs.rdstate();
    while (ifs.good())
    {
        ifs >> word;
        ifs.ignore(150, '\n');
        if (strcmp(word, "INPUT_PARAMETERS") == 0)
        {
            ierr = 1;
            break;
        }
        ifs.rdstate();
    }

    if (ierr == 0)
    {
        std::cout << " Error parameter list. "
                  << " The parameter list always starts with key word 'INPUT_PARAMETERS'. " << std::endl;
        ModuleBase::WARNING_QUIT("Input", "Bad parameter, please check the input parameters in file INPUT", 1);
    }

    ifs.rdstate();
    while (ifs.good())
    {
        ifs >> word1;
        if (ifs.eof())
            break;
        strtolower(word1, word);

        auto it = this->input_lists.find(word);
        if (it != this->input_lists.end())
        {
            Input_Item* p_item = &(it->second);
            this->readvalue_items.push_back(p_item);
            if (it->second.checkvalue)
                this->checkvalue_items.push_back(p_item);
            if (it->second.resetvalue)
                this->resetvalue_items.push_back(p_item);
            read_information(ifs, it->second.str_values, "#/!");
        }
        else
        {
            if (word[0] != '#' && word[0] != '/' && word[0] != '!')
            {
                std::cout << " THE PARAMETER NAME '" << word << "' IS NOT USED!" << std::endl;
                ModuleBase::WARNING_QUIT("Input", "Bad parameter, please check the input parameters in file INPUT", 1);
            }
            ifs.ignore(150, '\n');
        }
    }

    // 1) read the value of the parameters
    for (auto& readvalue_item: this->readvalue_items)
    {
        readvalue_item->readvalue(*readvalue_item, param);
    }

    // 2) check the value of the parameters
    for (auto& checkvalue_item: this->checkvalue_items)
    {
        checkvalue_item->checkvalue(*checkvalue_item, param);
    }

    // 3) reset the value of some parameters based on readin values
    //    e.g. if (calulation_type == "nscf") then set "init_chg" to file.
    for (auto& resetvalue_item: this->resetvalue_items)
    {
        resetvalue_item->resetvalue(*resetvalue_item, param);
    }
}

void ReadInput::add_item(const Input_Item& item)
{
    this->input_lists.insert(make_pair(item.label, item));
}

} // namespace ModuleIO