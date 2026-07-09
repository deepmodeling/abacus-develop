#ifndef ABACUS_SOURCE_RELAX_SOCKET_DRIVER_H
#define ABACUS_SOURCE_RELAX_SOCKET_DRIVER_H

#include "source_cell/unitcell.h"
#include "source_esolver/esolver.h"
#include "source_io/module_parameter/input_parameter.h"

#include <fstream>

class Socket_Driver
{
  public:
    Socket_Driver() = default;
    ~Socket_Driver() = default;

    void socket_driver(ModuleESolver::ESolver* p_esolver,
                       UnitCell& ucell,
                       const Input_para& inp,
                       std::ofstream& ofs_running);
};

#endif
