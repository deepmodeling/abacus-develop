#ifndef ABACUS_SOURCE_RELAX_SOCKET_DRIVER_H
#define ABACUS_SOURCE_RELAX_SOCKET_DRIVER_H

#include <iosfwd>

class UnitCell;
class Parameter;
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

    bool socket_driver(ModuleESolver::ESolver* solver, UnitCell& cell,
                       Parameter& parameters, std::ofstream& log);

    bool socket_driver(ModuleESolver::ESolver* p_esolver,
                       UnitCell& ucell,
                       const Input_para& inp,
                       std::ofstream& ofs_running);
};

#endif
