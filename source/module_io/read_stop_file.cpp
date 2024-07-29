#include "read_stop_file.h"

#include "module_io/read_input.h"

#include <fstream>
#include <iostream>
#include <sstream>

namespace ModuleIO
{

int read_stop_file(const std::string& filename, std::ofstream& ofs_running)
{
    std::ifstream ifs(filename.c_str(), std::ios::in);

    if (!ifs)
    {
        return 0;
    }

    ifs.clear();
    ifs.seekg(0);
    ifs.rdstate();
    int stop = 0;

    while (ifs.good())
    {
        std::string line;
        std::getline(ifs, line);
        if (line.empty())
        {
            continue;
        }
        std::istringstream iss(line);
        std::string word, result;
        iss >> word;
        if (iss.eof())
        {
            continue;
        }
        else
        {
            iss >> result;
        }

        bool value = convert_bool(result);

        if (word == "stop_ion" && value && stop == 0)
        {
            stop = 1;
        }
        else if (word == "stop_elec" && value && stop < 2)
        {
            stop = 2;
        }
    }
    if (stop == 1)
    {
        std::cout << "\n\n--------------------------------------------------------------------" << std::endl;
        std::cout << "--------------------------------------------------------------------" << std::endl;
        std::cout << " Read in stop_ion = true from " << filename << std::endl;
        std::cout << " The current execution stops at the ionic step " << std::endl;
        std::cout << "--------------------------------------------------------------------" << std::endl;
        std::cout << "--------------------------------------------------------------------\n" << std::endl;
        ofs_running << "\n\n--------------------------------------------------------------------" << std::endl;
        ofs_running << "--------------------------------------------------------------------" << std::endl;
        ofs_running << " Read in stop_ion = true from " << filename << std::endl;
        ofs_running << " The current execution stops at the ionic step " << std::endl;
        ofs_running << "--------------------------------------------------------------------" << std::endl;
        ofs_running << "--------------------------------------------------------------------\n" << std::endl;
    }
    else if (stop == 2)
    {
        std::cout << "\n\n--------------------------------------------------------------------" << std::endl;
        std::cout << "--------------------------------------------------------------------" << std::endl;
        std::cout << " Read in stop_elec = true from " << filename << std::endl;
        std::cout << " The current execution stops at the electronic step " << std::endl;
        std::cout << "--------------------------------------------------------------------" << std::endl;
        std::cout << "--------------------------------------------------------------------\n" << std::endl;
        ofs_running << "\n\n--------------------------------------------------------------------" << std::endl;
        ofs_running << "--------------------------------------------------------------------" << std::endl;
        ofs_running << " Read in stop_elec = true from " << filename << std::endl;
        ofs_running << " The current execution stops at the electronic step " << std::endl;
        ofs_running << "--------------------------------------------------------------------" << std::endl;
        ofs_running << "--------------------------------------------------------------------\n" << std::endl;
    }

    return stop;
}

} // namespace ModuleIO