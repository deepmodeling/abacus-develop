#include "socket_driver.h"

#include "source_esolver/esolver.h"
#include "source_relax/socket_handlers.h"

bool Socket_Driver::socket_driver(
    ModuleESolver::ESolver *solver,
    UnitCell &cell,
    const Input_para &input,
    std::ofstream &log)
{
    return SocketHandlers::run_socket_driver(solver, cell, input, nullptr, nullptr, nullptr, log);
}

bool Socket_Driver::socket_driver(
    ModuleESolver::ESolver *solver,
    UnitCell &cell,
    Parameter &parameters,
    const Input_para &input,
    const System_para &system,
    std::ofstream &log)
{
    SocketCP::ElectronicSlot *slot = dynamic_cast<SocketCP::ElectronicSlot *>(solver);
    elecstate::ElectronicStep **electronic_slot = nullptr;
    if (slot != nullptr)
    {
        electronic_slot = &slot->electronic_step;
    }
    return SocketHandlers::run_socket_driver(solver, cell, input, &parameters, &system, electronic_slot, log);
}
