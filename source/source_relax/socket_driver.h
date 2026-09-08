#ifndef ABACUS_SOURCE_RELAX_SOCKET_DRIVER_H
#define ABACUS_SOURCE_RELAX_SOCKET_DRIVER_H

#include <iosfwd>

class UnitCell;
struct Input_para;

namespace ModuleESolver
{
class ESolver;
}

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
