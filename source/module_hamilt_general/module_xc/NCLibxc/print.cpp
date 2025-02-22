#include "NCLibxc.h"
#include <iostream>
#include <fstream>
#include <iomanip>

namespace NCXC{

///////////////////////////////////////////////////////////////////////////////////
// print message and citation of the program
// Function to print message and citation
void NCLibxc::print_NCLibxc()
{
    std::ofstream log_file("NCLibxc.log", std::ios::out | std::ios::app);
    if (log_file.is_open())
    {
        log_file << "You are using the multi-collinear approach implemented by Xiaoyu Zhang. Please cite:\n";
        log_file.close();
    }
    else
    {
        std::cerr << "Unable to open log file." << std::endl;
    }
}

//////////////////////////////////////////////////////////////////////////////////
// print xc torque
void NCLibxc::print_torque(const std::vector<double>& torque)
{
    std::ofstream log_file("xc_torque.log", std::ios::out | std::ios::app);
    if (log_file.is_open())
    {
        log_file << "Torque values from gga_torque:\n";
        size_t num_points = torque.size() / 3;
        for (size_t i = 0; i < num_points; ++i)
        {
            log_file << "Point " << i << ": "
                     << "Tx = " << std::setw(12) << torque[3 * i] << ", "
                     << "Ty = " << std::setw(12) << torque[3 * i + 1] << ", "
                     << "Tz = " << std::setw(12) << torque[3 * i + 2] << "\n";
        }
        log_file.close();
    }
    else
    {
        std::cerr << "Unable to open xc_torque.log for writing." << std::endl;
    }

    
}

}