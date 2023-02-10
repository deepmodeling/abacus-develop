#ifndef RHO_IO_H
#define RHO_IO_H

#include<string>

namespace ModuleIO
{
    bool read_rho(const int &is, const std::string &fn, double* rho, int &prenspin);//mohan add 2007-10-17
    void write_rho(const double* rho_save, const int &is, const int &iter, const std::string &fn, 
		const int &precision = 11, const bool for_plot = false);//mohan add 2007-10-17
}

#endif
