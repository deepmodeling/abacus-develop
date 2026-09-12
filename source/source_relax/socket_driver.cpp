#include "socket_driver.h"

#include "source_relax/socket_handlers.h"
#include "source_relax/socket_utils.h"
#include "source_base/timer.h"
#include "source_cell/unitcell.h"
#include "source_esolver/esolver.h"
#include "source_io/module_parameter/input_parameter.h"

#include <exception>
#include <fstream>
#include <string>

using SocketUtils::ComputedFrame;
using SocketUtils::DriverState;
using SocketUtils::ipi_cell_bohr_from_unitcell;
using SocketUtils::is_root;
using SocketUtils::quit_if_root_failed;
using SocketUtils::address;
using SocketHandlers::DriverContext;

namespace
{
void connect_on_root(IpiSocket& socket, std::ofstream& ofs_running)
{
    int io_failed = 0;
    std::string io_message;
    if (is_root())
    {
        try
        {
            const std::string endpoint = address();
            ofs_running << " ABACUS socket driver connecting to i-PI endpoint " << endpoint << std::endl;
            socket.connect(endpoint);
        }
        catch (const std::exception& exc)
        {
            io_failed = 1;
            io_message = exc.what();
        }
    }
    quit_if_root_failed(io_failed, io_message);
}

void log_peer_closed(std::ofstream& ofs_running)
{
    if (is_root())
    {
        ofs_running << " ABACUS socket driver exiting after peer closed connection" << std::endl;
    }
}
} // namespace

void Socket_Driver::socket_driver(ModuleESolver::ESolver* p_esolver,
                                  UnitCell& ucell,
                                  const Input_para& inp,
                                  std::ofstream& ofs_running)
{
    ModuleBase::TITLE("Socket_Driver", "socket_driver");
    ModuleBase::timer::start("Socket_Driver", "socket_driver");

    if (p_esolver == nullptr)
    {
        ModuleBase::WARNING_QUIT("ABACUS socket", "socket driver requires a valid ESolver.");
    }
    IpiSocket socket;

    try
    {
        connect_on_root(socket, ofs_running);

        DriverContext context;
        context.esolver = p_esolver;
        context.ucell = &ucell;
        context.inp = &inp;
        context.state = DriverState::NeedInit;
        context.nat_return = ucell.nat;
        context.reference_cell = ipi_cell_bohr_from_unitcell(ucell);

        while (true)
        {
            const std::string header
                = SocketHandlers::read_header_bcast(socket, context.state);
            if (header.empty())
            {
                log_peer_closed(ofs_running);
                break;
            }
            else if (header == "STATUS")
            {
                SocketHandlers::handle_status(socket, context.state);
            }
            else if (header == "INIT")
            {
                SocketHandlers::handle_init(socket, context, ofs_running);
            }
            else if (header == "POSDATA")
            {
                SocketHandlers::handle_posdata(socket, context, ofs_running);
            }
            else if (header == "GETFORCE")
            {
                SocketHandlers::handle_getforce(socket, context);
            }
            else if (header == "EXIT")
            {
                SocketHandlers::handle_exit(ofs_running);
                break;
            }
            else
            {
                quit_if_root_failed(is_root() ? 1 : 0,
                                       is_root() ? "unknown i-PI header: " + header : "");
            }
        }
    }
    catch (const std::exception& exc)
    {
        ModuleBase::WARNING_QUIT("ABACUS socket", exc.what());
    }

    if (is_root())
    {
        socket.close();
    }

    ModuleBase::timer::end("Socket_Driver", "socket_driver");
}
