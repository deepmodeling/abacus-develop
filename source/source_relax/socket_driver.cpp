#include "socket_driver.h"

#include "source_esolver/esolver_fp.h"
#include "source_relax/socket_handlers.h"

bool Socket_Driver::socket_driver(ModuleESolver::ESolver* solver, UnitCell& cell,
                                 const Input_para& input, std::ofstream& log)
{
    return SocketHandlers::run_socket_driver(solver, cell, input, nullptr, nullptr, nullptr, log);
}

bool Socket_Driver::socket_driver(ModuleESolver::ESolver* solver, UnitCell& cell,
                                 Parameter& parameters, const Input_para& input,
                                 const System_para& system, std::ofstream& log)
{
    auto* fp_solver = dynamic_cast<ModuleESolver::ESolver_FP*>(solver);
    elecstate::ElectronicStep** electronic_slot = nullptr;
    if (fp_solver != nullptr)
    {
        electronic_slot = &fp_solver->electronic_step;
    }
    return SocketHandlers::run_socket_driver(solver, cell, input, &parameters, &system, electronic_slot, log);
}
