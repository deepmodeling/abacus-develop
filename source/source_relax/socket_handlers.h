#ifndef ABACUS_SOURCE_RELAX_SOCKET_HANDLERS_H
#define ABACUS_SOURCE_RELAX_SOCKET_HANDLERS_H

#include "source_relax/socket_utils.h"
#include "source_relax/socket_ipi.h"

#include <fstream>
#include <string>
#include <vector>

class UnitCell;
class Parameter;

namespace ModuleESolver
{
class ESolver;
}
class Input_para;

namespace SocketHandlers
{
struct DriverContext
{
    ModuleESolver::ESolver* esolver = nullptr;
    UnitCell* ucell = nullptr;
    const Input_para* inp = nullptr;
    Parameter* parameters = nullptr;
    SocketUtils::DriverState state = SocketUtils::DriverState::NeedInit;
    int istep = 0;
    int nat_return = 0;
    SocketUtils::ComputedFrame published;
    std::vector<double> reference_cell;
    bool checked_initial_positions = false;
    SocketCP::Configuration cp;
    ModuleESolver::ElectronicControl pending_control;
    bool charge_pending = false;
    std::vector<double> cp_wire_cell;
};

// Reads the next i-PI header on the root rank and broadcasts it to all ranks.
// Returns an empty string when the peer closed the connection while no frame
// is pending. Calls WARNING_QUIT on unrecoverable I/O failure.
std::string read_header_bcast(IpiSocket& socket, const DriverContext& context);

void handle_status(IpiSocket& socket, const SocketUtils::DriverState state);

void handle_init(IpiSocket& socket,
                 DriverContext& context,
                 std::ofstream& ofs_running);

void handle_chgdata(IpiSocket& socket, DriverContext& context, std::ofstream& ofs_running);

void handle_posdata(IpiSocket& socket,
                    DriverContext& context,
                    std::ofstream& ofs_running);

void handle_getforce(IpiSocket& socket, DriverContext& context, std::ofstream& ofs_running);

void handle_exit(const DriverContext& context, std::ofstream& ofs_running);
} // namespace SocketHandlers

#endif
