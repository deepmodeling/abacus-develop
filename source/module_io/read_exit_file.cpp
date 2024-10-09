#include "read_exit_file.h"

#include <fstream>
#include <iostream>
#include <sstream>

#ifdef __MPI
#include <mpi.h>
#endif

namespace ModuleIO
{

int read_exit_file(const int& my_rank, const std::string& filename, std::ofstream& ofs_running)
{
    auto str2bool = [](std::string str) {
        for (auto& i: str)
        {
            i = tolower(i);
        }
        if (str == "true" || str == "t" || str == "1")
        {
            return true;
        }
        else
        {
            return false;
        }
    };
    int stop = 0;
    if (my_rank == 0)
    {
        std::ifstream ifs(filename.c_str(), std::ios::in);

        if (ifs)
        {
            ifs.clear();
            ifs.seekg(0);
            ifs.rdstate();

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

                if (word == "stop_ion" && str2bool(result) && stop < 1)
                {
                    stop = 1;
                }
                else if (word == "stop_elec" && str2bool(result) && stop < 2)
                {
                    stop = 2;
                }
            }
        }
    }

#ifdef __MPI
    MPI_Bcast(&stop, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

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