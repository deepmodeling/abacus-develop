#include "socket_handlers.h"

#include "source_relax/socket_frame.h"
#include "source_base/global_function.h"
#include "source_cell/unitcell.h"
#include "source_esolver/esolver.h"
#include "source_io/module_parameter/input_parameter.h"
#ifdef __JSON
#include "source_io/module_json/output_info.h"
#endif

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace SocketHandlers
{
using SocketUtils::ComputedFrame;
using SocketUtils::DriverState;
using SocketUtils::kInverseAbsoluteTolerance;
using SocketUtils::kInverseRelativeTolerance;
using SocketUtils::kMaxCellCondition;
using SocketUtils::kMaxInitBytes;
using SocketUtils::kRyToHartree;
using SocketUtils::kStressAbsoluteTolerance;
using SocketUtils::kStressRelativeTolerance;
using SocketUtils::all_ranks_converged;
using SocketUtils::bcast_double_vector;
using SocketUtils::bcast_header;
using SocketUtils::bcast_int32;
using SocketUtils::bcast_string;
using SocketUtils::fail_during_collective_stage;
using SocketUtils::flatten_forces_hartree_per_bohr;
using SocketUtils::ipi_cell_bohr_from_unitcell;
using SocketUtils::is_root;
using SocketUtils::matrix9_from_stress;
using SocketUtils::max_abs_delta;
using SocketUtils::max_wrapped_direct_delta_from_unitcell;
using SocketUtils::properties_extra;
using SocketUtils::quit_if_root_failed;
using SocketUtils::set_positions_from_ipi_bohr;
using SocketUtils::throw_if_any_rank_failed;
using SocketUtils::unchanged_cell_tolerance;
using SocketUtils::vector_from_matrix9;

std::string read_header_bcast(IpiSocket& socket, const DriverContext& context)
{
    std::string header;
    int io_failed = 0;
    std::string io_message;
    if (is_root())
    {
        try
        {
            header = socket.read_header();
        }
        catch (const IpiSocketClosed&)
        {
            if (context.state == DriverState::HasData || context.charge_pending)
            {
                io_failed = 1;
                io_message = "i-PI peer closed while a computed frame was pending";
            }
        }
        catch (const std::exception& exc)
        {
            io_failed = 1;
            io_message = exc.what();
        }
    }
    quit_if_root_failed(io_failed, io_message);
    return bcast_header(header);
}

void handle_status(IpiSocket& socket, const DriverState state)
{
    int io_failed = 0;
    std::string io_message;
    if (is_root())
    {
        try
        {
            if (state == DriverState::HasData)
            {
                socket.write_header("HAVEDATA");
            }
            else if (state == DriverState::Ready)
            {
                socket.write_header("READY");
            }
            else
            {
                socket.write_header("NEEDINIT");
            }
        }
        catch (const std::exception& exc)
        {
            io_failed = 1;
            io_message = exc.what();
        }
    }
    quit_if_root_failed(io_failed, io_message);
}

void handle_init(IpiSocket& socket,
                 DriverContext& context,
                 std::ofstream& ofs_running)
{
    std::int32_t rid = 0;
    std::int32_t nbytes = 0;
    std::string params;
    int io_failed = 0;
    std::string io_message;
    if (is_root())
    {
        if (context.state != DriverState::NeedInit)
        {
            io_failed = 1;
            io_message = "INIT requires NEEDINIT state";
        }
        else
        {
            try
            {
                rid = socket.read_int32();
                nbytes = socket.read_int32();
                if (nbytes < 0)
                {
                    io_failed = 1;
                    io_message = "negative INIT payload length from i-PI socket";
                }
                else if (nbytes > kMaxInitBytes)
                {
                    io_failed = 1;
                    io_message = "INIT payload exceeds the 1 MiB socket limit";
                }
                else if (nbytes > 0)
                {
                    params = socket.read_string(static_cast<std::size_t>(nbytes));
                }
            }
            catch (const std::exception& exc)
            {
                io_failed = 1;
                io_message = exc.what();
            }
        }
    }
    quit_if_root_failed(io_failed, io_message);
    bcast_int32(rid);
    bcast_int32(nbytes);
    if (nbytes > 0 && is_root())
    {
        ofs_running << " ABACUS socket INIT params bytes " << nbytes << std::endl;
    }
    bcast_string(params);
    int local_failed = 0;
    std::string local_message;
    try
    {
        context.cp = SocketCP::parse_init(params);
        if (context.cp.enabled && (context.parameters == nullptr || !context.inp->cal_force
            || context.inp->cal_stress || (context.inp->symmetry != "0" && context.inp->symmetry != "-1")))
        {
            throw std::runtime_error("constant-potential socket requires runtime parameters, "
                                     "cal_force=1, cal_stress=0 and symmetry=0/-1");
        }
    }
    catch (const std::exception& exc)
    {
        local_failed = 1;
        local_message = exc.what();
    }
    throw_if_any_rank_failed(local_failed, local_message);
    context.state = DriverState::Ready;
    if (is_root())
    {
        ofs_running << " ABACUS socket INIT replica " << rid << std::endl;
    }
}

void handle_chgdata(IpiSocket& socket, DriverContext& context, std::ofstream& ofs_running)
{
    std::vector<double> charge(1, 0.0);
    std::int32_t flags = 0;
    int io_failed = 0;
    std::string io_message;
    if (is_root())
    {
        try
        {
            if (context.state != DriverState::Ready || context.charge_pending)
            {
                throw std::runtime_error("CHGDATA requires READY with no pending charge or force result");
            }
            charge[0] = socket.read_double();
            flags = socket.read_int32();
        }
        catch (const std::exception& exc)
        {
            io_failed = 1;
            io_message = exc.what();
        }
    }
    quit_if_root_failed(io_failed, io_message);
    bcast_double_vector(charge);
    bcast_int32(flags);
    int local_failed = 0;
    std::string local_message;
    try
    {
        context.pending_control = SocketCP::make_control(context.cp, charge[0], flags);
    }
    catch (const std::exception& exc)
    {
        local_failed = 1;
        local_message = exc.what();
    }
    throw_if_any_rank_failed(local_failed, local_message);
    context.charge_pending = true;
    if (is_root())
    {
        ofs_running << " @IPI DRIVER: CHGDATA nelect=" << std::setprecision(17)
                    << charge[0] << " flag=" << flags << std::endl;
    }
}

namespace
{
struct PosdataPayload
{
    SocketFrame::Matrix9 cell = {{0.0}};
    SocketFrame::Matrix9 inv_cell = {{0.0}};
    std::int32_t nat_socket = 0;
    std::vector<double> positions;
};

// Root rank reads and validates the POSDATA frame; the results are then
// broadcast to all ranks. Calls WARNING_QUIT on protocol/validation failure.
// The READY-state check must stay inside this function so that every rank
// passes through the single quit_if_root_failed collective sequence below;
// a root-only early quit outside would deadlock if any collective were
// added ahead of it.
PosdataPayload read_posdata(IpiSocket& socket, const DriverContext& context)
{
    PosdataPayload payload;
    int io_failed = 0;
    std::string io_message;
    if (is_root())
    {
        if (context.state != DriverState::Ready)
        {
            io_failed = 1;
            io_message = "POSDATA requires READY state";
        }
        else if (context.cp.enabled && !context.charge_pending)
        {
            io_failed = 1;
            io_message = "constant-potential POSDATA requires one fresh CHGDATA";
        }
        else
        {
            try
            {
                const std::vector<double> cell_values = socket.read_doubles(9);
                const std::vector<double> inverse_values = socket.read_doubles(9);
                std::copy(cell_values.begin(), cell_values.end(), payload.cell.begin());
                std::copy(inverse_values.begin(), inverse_values.end(), payload.inv_cell.begin());
                payload.nat_socket = socket.read_int32();
                const SocketFrame::CellValidation validation
                    = SocketFrame::validate_ipi_cell(payload.cell,
                                                     payload.inv_cell,
                                                     kMaxCellCondition,
                                                     kInverseAbsoluteTolerance,
                                                     kInverseRelativeTolerance);
                if (!validation.ok)
                {
                    io_failed = 1;
                    io_message = "invalid POSDATA cell: " + validation.message;
                }
                std::size_t coordinate_count = 0;
                if (io_failed == 0
                    && !SocketFrame::checked_position_count(payload.nat_socket,
                                                            context.ucell->nat,
                                                            coordinate_count,
                                                            io_message))
                {
                    io_failed = 1;
                }
                if (io_failed == 0)
                {
                    payload.positions = socket.read_doubles(coordinate_count);
                    if (!SocketFrame::validate_positions(payload.positions,
                                                         coordinate_count,
                                                         io_message))
                    {
                        io_failed = 1;
                    }
                }
            }
            catch (const std::exception& exc)
            {
                io_failed = 1;
                io_message = exc.what();
            }
        }
    }
    quit_if_root_failed(io_failed, io_message);
    return payload;
}

void bcast_posdata(PosdataPayload& payload)
{
    bcast_int32(payload.nat_socket);
    std::vector<double> cell_values(payload.cell.begin(), payload.cell.end());
    std::vector<double> inverse_values(payload.inv_cell.begin(), payload.inv_cell.end());
    bcast_double_vector(cell_values);
    bcast_double_vector(inverse_values);
    if (!is_root())
    {
        payload.cell = {{0.0}};
        payload.inv_cell = {{0.0}};
        std::copy(cell_values.begin(), cell_values.end(), payload.cell.begin());
        std::copy(inverse_values.begin(), inverse_values.end(), payload.inv_cell.begin());
        if (payload.nat_socket >= 0)
        {
            payload.positions.assign(static_cast<std::size_t>(3 * payload.nat_socket), 0.0);
        }
    }
    bcast_double_vector(payload.positions);
}

void check_posdata_geometry(DriverContext& context, const PosdataPayload& payload)
{
    const std::vector<double> wire_cell(payload.cell.begin(), payload.cell.end());
    const auto& expected_cell = context.cp.enabled && !context.cp_wire_cell.empty()
                                    ? context.cp_wire_cell : context.reference_cell;
    double cell_tolerance = unchanged_cell_tolerance(payload.cell);
    if (context.cp.enabled && context.cp_wire_cell.empty())
    {
        // Allow the two applications' rounded length constants once; lock later wire cells.
        double scale = 1.0;
        for (const double value : context.reference_cell) { scale = std::max(scale, std::abs(value)); }
        cell_tolerance = 1.0e-6 * scale;
    }
    if (max_abs_delta(wire_cell, expected_cell) > cell_tolerance)
    {
        ModuleBase::WARNING_QUIT("ABACUS socket",
                                 "variable-cell socket updates are not supported yet.");
    }
    if (context.cp.enabled && context.cp_wire_cell.empty()) { context.cp_wire_cell = wire_cell; }
    if (!context.checked_initial_positions)
    {
        context.checked_initial_positions = true;
        if (max_wrapped_direct_delta_from_unitcell(*context.ucell, payload.positions) > 1.0e-5
            && is_root())
        {
            ModuleBase::WARNING(
                "ABACUS socket",
                "first POSDATA positions are not PBC-equivalent to STRU atom order; "
                "i-PI POSDATA carries no species, so the client atoms should use the same atom order "
                "as STRU.");
        }
    }
}

void run_esolver_for_positions(DriverContext& context, const std::vector<double>& positions)
{
    try
    {
        set_positions_from_ipi_bohr(*context.ucell, positions);
    }
    catch (const std::exception& exc)
    {
        fail_during_collective_stage("set_positions", exc.what());
    }
    catch (...)
    {
        fail_during_collective_stage("set_positions", "unknown socket position update failure");
    }
    try
    {
        if (context.cp.enabled)
        {
            context.charge_pending = false;
            if (!context.esolver->apply_electronic_control(*context.parameters, context.pending_control))
            {
                throw std::runtime_error("selected ESolver does not support constant-potential electron control");
            }
        }
        context.esolver->runner(*context.ucell, context.istep);
    }
    catch (const std::exception& exc)
    {
        fail_during_collective_stage("runner", exc.what());
    }
    catch (...)
    {
        fail_during_collective_stage("runner", "unknown socket runner failure");
    }
}

void compute_energy_hartree(ModuleESolver::ESolver* esolver,
                            ComputedFrame& computed,
                            std::ofstream& ofs_running)
{
    double energy_ry = 0.0;
    try
    {
        energy_ry = esolver->cal_energy();
    }
    catch (const std::exception& exc)
    {
        fail_during_collective_stage("cal_energy", exc.what());
    }
    catch (...)
    {
        fail_during_collective_stage("cal_energy", "unknown socket energy failure");
    }
    const int local_failed = std::isfinite(energy_ry) ? 0 : 1;
    throw_if_any_rank_failed(local_failed,
                             local_failed == 0 ? "" : "socket energy is not finite");
    if (!std::isfinite(energy_ry))
    {
        ModuleBase::WARNING_QUIT("ABACUS socket", "socket energy is not finite.");
    }
    computed.energy_hartree = energy_ry * kRyToHartree;
    if (is_root())
    {
        ofs_running << " ABACUS socket return energy "
                    << energy_ry << " Ry, "
                    << energy_ry * ModuleBase::Ry_to_eV << " eV, "
                    << computed.energy_hartree << " Ha" << std::endl;
    }
}

void compute_forces(UnitCell& ucell,
                    ModuleESolver::ESolver* esolver,
                    ComputedFrame& computed)
{
    ModuleBase::matrix force;
    try
    {
        esolver->cal_force(ucell, force);
    }
    catch (const std::exception& exc)
    {
        fail_during_collective_stage("cal_force", exc.what());
    }
    catch (...)
    {
        fail_during_collective_stage("cal_force", "unknown socket force failure");
    }
    int local_failed = 0;
    std::string local_message;
    try
    {
        computed.forces_hartree_per_bohr = flatten_forces_hartree_per_bohr(force, ucell.nat);
    }
    catch (const std::exception& exc)
    {
        local_failed = 1;
        local_message = exc.what();
    }
    catch (...)
    {
        local_failed = 1;
        local_message = "unknown socket force validation failure";
    }
    throw_if_any_rank_failed(local_failed, local_message);
    computed.forces_present = true;
}

void compute_stress(UnitCell& ucell,
                    ModuleESolver::ESolver* esolver,
                    ComputedFrame& computed)
{
    ModuleBase::matrix stress;
    try
    {
        esolver->cal_stress(ucell, stress);
    }
    catch (const std::exception& exc)
    {
        fail_during_collective_stage("cal_stress", exc.what());
    }
    catch (...)
    {
        fail_during_collective_stage("cal_stress", "unknown socket stress failure");
    }
    int local_failed = 0;
    std::string local_message;
    try
    {
        const SocketFrame::VirialConversion virial
            = SocketFrame::make_ipi_virial(matrix9_from_stress(stress),
                                           ucell.omega,
                                           kStressAbsoluteTolerance,
                                           kStressRelativeTolerance);
        if (!virial.ok)
        {
            throw std::runtime_error(virial.message);
        }
        computed.virial_wire_hartree = virial.wire_virial_hartree;
    }
    catch (const std::exception& exc)
    {
        local_failed = 1;
        local_message = exc.what();
    }
    catch (...)
    {
        local_failed = 1;
        local_message = "unknown socket stress validation failure";
    }
    throw_if_any_rank_failed(local_failed, local_message);
    computed.stress_present = true;
}
} // namespace

void handle_posdata(IpiSocket& socket,
                    DriverContext& context,
                    std::ofstream& ofs_running)
{
    PosdataPayload payload = read_posdata(socket, context);
    bcast_posdata(payload);
    check_posdata_geometry(context, payload);
#ifdef __JSON
    // Each POSDATA frame is one electronic-structure output step.
    Json::init_output_array_obj();
#endif
    run_esolver_for_positions(context, payload.positions);

    ComputedFrame computed;
    computed.scf_converged = all_ranks_converged(context.esolver->conv_esolver);
    if (context.cp.enabled && !computed.scf_converged)
    {
        throw std::runtime_error("SCF did not converge; refusing to return constant-potential results");
    }
    if (!computed.scf_converged && is_root())
    {
        ModuleBase::WARNING(
            "ABACUS socket",
            "SCF did not converge; returning the available frame and marking it in i-PI extras.");
    }
    compute_energy_hartree(context.esolver, computed, ofs_running);
    if (context.inp->cal_force)
    {
        compute_forces(*context.ucell, context.esolver, computed);
    }
    if (context.inp->cal_stress)
    {
        compute_stress(*context.ucell, context.esolver, computed);
    }
    if (context.cp.enabled)
    {
        try
        {
            computed.electronic = SocketCP::collect_observables(context.cp, *context.esolver,
                                                                *context.ucell, *context.parameters);
        }
        catch (const std::exception& exc)
        {
            fail_during_collective_stage("electronic_observables", exc.what());
        }
        catch (...)
        {
            fail_during_collective_stage("electronic_observables", "unknown electronic observables failure");
        }
        computed.electronic_present = true;
    }
    computed.valid = true;
    context.published = computed;
    ++context.istep;
    context.state = DriverState::HasData;
}

void handle_getforce(IpiSocket& socket, DriverContext& context, std::ofstream& ofs_running)
{
    int io_failed = 0;
    std::string io_message;
    if (is_root())
    {
        try
        {
            if (context.state != DriverState::HasData || !context.published.valid)
            {
                throw std::runtime_error("GETFORCE requires HAVEDATA state and a valid frame");
            }
            socket.write_header("FORCEREADY");
            socket.write_double(context.published.energy_hartree);
            socket.write_int32(static_cast<std::int32_t>(context.nat_return));
            const std::vector<double> forces
                = context.published.forces_present
                      ? context.published.forces_hartree_per_bohr
                      : std::vector<double>(static_cast<std::size_t>(3 * context.nat_return), 0.0);
            socket.write_doubles(forces);
            socket.write_doubles(vector_from_matrix9(context.published.virial_wire_hartree));
            const std::string extra = properties_extra(context.published);
            if (extra.size() > static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max()))
            {
                throw std::overflow_error("i-PI extras payload is larger than int32");
            }
            socket.write_int32(static_cast<std::int32_t>(extra.size()));
            socket.write_string(extra);
            ofs_running << " @IPI DRIVER: Sending extras JSON = " << extra << std::endl;
        }
        catch (const std::exception& exc)
        {
            io_failed = 1;
            io_message = exc.what();
        }
    }
    quit_if_root_failed(io_failed, io_message);
    context.published = ComputedFrame();
    context.state = DriverState::Ready;
}

void handle_exit(const DriverContext& context, std::ofstream& ofs_running)
{
    const bool unfinished = context.cp.enabled
                            && (context.state == DriverState::HasData || context.charge_pending);
    throw_if_any_rank_failed(unfinished ? 1 : 0, "EXIT with an unfinished request/result");
    if (is_root())
    {
        ofs_running << " ABACUS socket driver received i-PI EXIT" << std::endl;
    }
}
} // namespace SocketHandlers
