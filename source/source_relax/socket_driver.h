#ifndef ABACUS_SOURCE_RELAX_SOCKET_DRIVER_H
#define ABACUS_SOURCE_RELAX_SOCKET_DRIVER_H

#include <iosfwd>

class UnitCell;
class Parameter;
struct Input_para;
struct System_para;

namespace ModuleESolver
{
class ESolver;
}

namespace elecstate
{
struct ElectronicStep;
}

namespace SocketCP
{
// Implemented only by socket KS adapters. The outer driver owns each request
// and binds this rank-local, non-owning slot for one evaluation.
struct ElectronicSlot
{
    virtual ~ElectronicSlot() = default;
    elecstate::ElectronicStep *electronic_step = nullptr;
};

ModuleESolver::ESolver *create_pw_socket_solver();
#ifdef __LCAO
ModuleESolver::ESolver *create_lcao_socket_solver(bool gamma_only);
#endif
} // namespace SocketCP

class Socket_Driver
{
  public:
    Socket_Driver() = default;
    ~Socket_Driver() = default;

    bool socket_driver(
        ModuleESolver::ESolver *solver,
        UnitCell &cell,
        Parameter &parameters,
        const Input_para &input,
        const System_para &system,
        std::ofstream &log);

    bool socket_driver(
        ModuleESolver::ESolver *p_esolver,
        UnitCell &ucell,
        const Input_para &inp,
        std::ofstream &ofs_running);
};

#endif
