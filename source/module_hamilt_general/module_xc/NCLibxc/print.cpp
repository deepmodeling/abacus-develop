#include "NCLibxc.h"
#include <iostream>
#include <fstream>
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