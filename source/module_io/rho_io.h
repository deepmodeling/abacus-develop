#ifndef RHO_IO_H
#define RHO_IO_H
#include <string>
#include "module_cell/unitcell.h"

namespace ModuleIO
{
    bool read_rho(const int &is,
		    const int &nspin,
		    const std::string &fn,
		    double* rho,
		    int& nx,
		    int& ny,
		    int& nz,
		    double& ef,
		    UnitCell& ucell,
		    int &prenspin);//mohan add 2007-10-17
    void write_rho(const double* rho_save,
		    const int &is,
		    const int &iter,
		    const std::string &fn,
		    const int &precision = 11,
		    const bool for_plot = false);//mohan add 2007-10-17
}

#endif
