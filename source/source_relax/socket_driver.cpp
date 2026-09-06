#include "socket_driver.h"

#include "source_relax/socket_ipi.h"
#include "source_relax/socket_frame.h"
#include "source_base/global_function.h"
#include "source_base/mathzone.h"
#include "source_base/parallel_common.h"
#include "source_base/timer.h"
#include "source_cell/unitcell.h"
#include "source_cell/update_cell.h"
#include "source_esolver/esolver.h"
#include "source_io/module_parameter/input_parameter.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
constexpr double RY_TO_HARTREE = 0.5;
constexpr int IPI_RANK_ROOT = 0;
constexpr double MAX_CELL_CONDITION = 1.0e12;
constexpr double INVERSE_ABSOLUTE_TOLERANCE
    = 64.0 * std::numeric_limits<double>::epsilon();
constexpr double INVERSE_RELATIVE_TOLERANCE = 64.0;
constexpr double STRESS_ABSOLUTE_TOLERANCE = 1.0e-10;
constexpr double STRESS_RELATIVE_TOLERANCE = 1.0e-8;
constexpr std::int32_t MAX_INIT_BYTES = INT32_C(1048576);

enum class DriverState
{
    NeedInit,
    Ready,
    HasData
};

struct ComputedFrame
{
    bool valid = false;
    bool forces_present = false;
    bool stress_present = false;
    bool scf_converged = true;
    double energy_hartree = 0.0;
    std::vector<double> forces_hartree_per_bohr;
    SocketFrame::Matrix9 virial_wire_hartree = {{0.0}};
};

bool all_ranks_converged(const bool local_converged)
{
    int converged = local_converged ? 1 : 0;
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, &converged, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
#endif
    return converged != 0;
}

void throw_if_any_rank_failed(int local_failed, std::string local_message)
{
    int any_failed = local_failed;
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, &any_failed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
#endif
    if (any_failed != 0)
    {
        if (local_message.empty())
        {
            local_message = "socket frame validation failed on another MPI rank";
        }
        throw std::runtime_error(local_message);
    }
}

[[noreturn]] void fail_during_collective_stage(const char* stage,
                                               const std::string& message)
{
#ifdef __MPI
    int rank = -1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::fprintf(stderr,
                 "ABACUS_SOCKET_MPI_FATAL stage=%s rank=%d message=%s\n",
                 stage,
                 rank,
                 message.c_str());
    std::fflush(stderr);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    std::abort();
#else
    (void)stage;
    throw std::runtime_error(message);
#endif
}

std::string properties_extra(const ComputedFrame& frame)
{
    std::ostringstream extra;
    extra << "{\"schema\":\"abacus.socket.properties.v1\",\"present\":[\"energy\"";
    if (frame.forces_present)
    {
        extra << ",\"forces\"";
    }
    if (frame.stress_present)
    {
        extra << ",\"stress\"";
    }
    extra << "],\"scf_converged\":"
          << (frame.scf_converged ? "true" : "false") << "}";
    return extra.str();
}

bool is_root()
{
#ifdef __MPI
    int rank = IPI_RANK_ROOT;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return rank == IPI_RANK_ROOT;
#else
    return true;
#endif
}

void bcast_double_vector(std::vector<double>& values)
{
#ifdef __MPI
    if (!values.empty())
    {
        Parallel_Common::bcast_double(values.data(), static_cast<int>(values.size()));
    }
#else
    (void)values;
#endif
}

void bcast_socket_int(int& value)
{
#ifdef __MPI
    Parallel_Common::bcast_int(value);
#else
    (void)value;
#endif
}

void bcast_socket_int32(std::int32_t& value)
{
#ifdef __MPI
    MPI_Bcast(&value, 1, MPI_INT32_T, IPI_RANK_ROOT, MPI_COMM_WORLD);
#else
    (void)value;
#endif
}

void bcast_socket_chars(char* value, const int size)
{
#ifdef __MPI
    Parallel_Common::bcast_char(value, size);
#else
    (void)value;
    (void)size;
#endif
}

void bcast_socket_string(std::string& value)
{
    int size = static_cast<int>(value.size());
    bcast_socket_int(size);
    if (!is_root())
    {
        value.resize(static_cast<std::size_t>(size));
    }
    if (size > 0)
    {
        bcast_socket_chars(&value[0], size);
    }
}

void quit_if_root_io_failed(int root_failed, std::string root_message)
{
    bcast_socket_int(root_failed);
    bcast_socket_string(root_message);
    if (root_failed != 0)
    {
        ModuleBase::WARNING_QUIT("ABACUS socket", root_message.empty() ? "i-PI socket I/O failed" : root_message);
    }
}

std::string bcast_header(std::string header)
{
    bcast_socket_string(header);
    return header;
}

std::string socket_address()
{
    const char* env = std::getenv("ABACUS_SOCKET_ADDRESS");
    if (env == nullptr || std::string(env).empty())
    {
        return "localhost:31415";
    }
    return std::string(env);
}

std::vector<double> ipi_cell_bohr_from_unitcell(const UnitCell& ucell)
{
    const double lat0 = ucell.lat0;
    // ASE/i-PI sends POSDATA cell as cell.T in C order. ABACUS stores
    // lattice vectors as rows in latvec, so use the transposed order here.
    return {
        ucell.latvec.e11 * lat0, ucell.latvec.e21 * lat0, ucell.latvec.e31 * lat0,
        ucell.latvec.e12 * lat0, ucell.latvec.e22 * lat0, ucell.latvec.e32 * lat0,
        ucell.latvec.e13 * lat0, ucell.latvec.e23 * lat0, ucell.latvec.e33 * lat0,
    };
}

double max_wrapped_direct_delta_from_unitcell(const UnitCell& ucell, const std::vector<double>& positions_bohr)
{
    if (positions_bohr.size() != static_cast<std::size_t>(3 * ucell.nat))
    {
        return 1.0e99;
    }

    double out = 0.0;
    int iat = 0;
    for (int it = 0; it < ucell.ntype; ++it)
    {
        const Atom* atom = &ucell.atoms[it];
        for (int ia = 0; ia < atom->na; ++ia)
        {
            const double tau_x = positions_bohr[3 * iat + 0] / ucell.lat0;
            const double tau_y = positions_bohr[3 * iat + 1] / ucell.lat0;
            const double tau_z = positions_bohr[3 * iat + 2] / ucell.lat0;

            double dx = 0.0;
            double dy = 0.0;
            double dz = 0.0;
            ModuleBase::Mathzone::Cartesian_to_Direct(tau_x,
                                                      tau_y,
                                                      tau_z,
                                                      ucell.latvec.e11,
                                                      ucell.latvec.e12,
                                                      ucell.latvec.e13,
                                                      ucell.latvec.e21,
                                                      ucell.latvec.e22,
                                                      ucell.latvec.e23,
                                                      ucell.latvec.e31,
                                                      ucell.latvec.e32,
                                                      ucell.latvec.e33,
                                                      dx,
                                                      dy,
                                                      dz);

            double ddx = dx - atom->taud[ia].x;
            double ddy = dy - atom->taud[ia].y;
            double ddz = dz - atom->taud[ia].z;
            ddx -= std::round(ddx);
            ddy -= std::round(ddy);
            ddz -= std::round(ddz);
            out = std::max(out, std::abs(ddx));
            out = std::max(out, std::abs(ddy));
            out = std::max(out, std::abs(ddz));
            ++iat;
        }
    }
    return out;
}

double max_abs_delta(const std::vector<double>& a, const std::vector<double>& b)
{
    if (a.size() != b.size())
    {
        return 1.0e99;
    }
    double out = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i)
    {
        out = std::max(out, std::abs(a[i] - b[i]));
    }
    return out;
}

double unchanged_cell_tolerance(const SocketFrame::Matrix9& cell)
{
    double maximum = 0.0;
    for (std::size_t index = 0; index < cell.size(); ++index)
    {
        maximum = std::max(maximum, std::fabs(cell[index]));
    }
    return 32.0 * std::numeric_limits<double>::epsilon() * std::max(1.0, maximum);
}

void set_positions_from_ipi_bohr(UnitCell& ucell, const std::vector<double>& positions_bohr)
{
    if (positions_bohr.size() != static_cast<std::size_t>(3 * ucell.nat))
    {
        ModuleBase::WARNING_QUIT("ABACUS socket", "POSDATA atom count does not match STRU.");
    }

    int iat = 0;
    for (int it = 0; it < ucell.ntype; ++it)
    {
        Atom* atom = &ucell.atoms[it];
        for (int ia = 0; ia < atom->na; ++ia)
        {
            const double tau_x = positions_bohr[3 * iat + 0] / ucell.lat0;
            const double tau_y = positions_bohr[3 * iat + 1] / ucell.lat0;
            const double tau_z = positions_bohr[3 * iat + 2] / ucell.lat0;

            double dx = 0.0;
            double dy = 0.0;
            double dz = 0.0;
            ModuleBase::Mathzone::Cartesian_to_Direct(tau_x,
                                                      tau_y,
                                                      tau_z,
                                                      ucell.latvec.e11,
                                                      ucell.latvec.e12,
                                                      ucell.latvec.e13,
                                                      ucell.latvec.e21,
                                                      ucell.latvec.e22,
                                                      ucell.latvec.e23,
                                                      ucell.latvec.e31,
                                                      ucell.latvec.e32,
                                                      ucell.latvec.e33,
                                                      dx,
                                                      dy,
                                                      dz);

            atom->dis[ia].x = dx - atom->taud[ia].x;
            atom->dis[ia].y = dy - atom->taud[ia].y;
            atom->dis[ia].z = dz - atom->taud[ia].z;
            atom->taud[ia].x = dx;
            atom->taud[ia].y = dy;
            atom->taud[ia].z = dz;
            atom->tau[ia].x = tau_x;
            atom->tau[ia].y = tau_y;
            atom->tau[ia].z = tau_z;
            ++iat;
        }
    }
    unitcell::periodic_boundary_adjustment(ucell.atoms, ucell.latvec, ucell.ntype);
    ucell.ionic_position_updated = true;
    ucell.cell_parameter_updated = false;
}

std::vector<double> flatten_forces_hartree_per_bohr(const ModuleBase::matrix& force, const int nat)
{
    if (nat < 0 || force.nr != nat || force.nc != 3)
    {
        throw std::runtime_error("force matrix must have nat rows and three columns");
    }
    std::vector<double> out(static_cast<std::size_t>(force.nr * force.nc));
    for (int iat = 0; iat < force.nr; ++iat)
    {
        for (int idir = 0; idir < force.nc; ++idir)
        {
            const double value = force(iat, idir);
            if (!std::isfinite(value))
            {
                throw std::runtime_error("force entries must be finite");
            }
            out[static_cast<std::size_t>(3 * iat + idir)] = value * RY_TO_HARTREE;
        }
    }
    return out;
}

SocketFrame::Matrix9 matrix9_from_stress(const ModuleBase::matrix& stress)
{
    if (stress.nr != 3 || stress.nc != 3)
    {
        throw std::runtime_error("stress matrix must have three rows and three columns");
    }
    SocketFrame::Matrix9 values;
    for (int row = 0; row < 3; ++row)
    {
        for (int column = 0; column < 3; ++column)
        {
            values[3 * row + column] = stress(row, column);
        }
    }
    return values;
}

std::vector<double> vector_from_matrix9(const SocketFrame::Matrix9& values)
{
    return std::vector<double>(values.begin(), values.end());
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
        int io_failed = 0;
        std::string io_message;
        if (is_root())
        {
            try
            {
                const std::string address = socket_address();
                ofs_running << " ABACUS socket driver connecting to i-PI endpoint " << address << std::endl;
                socket.connect(address);
            }
            catch (const std::exception& exc)
            {
                io_failed = 1;
                io_message = exc.what();
            }
        }
        quit_if_root_io_failed(io_failed, io_message);

        DriverState state = DriverState::NeedInit;
        int istep = 0;
        const int nat_return = ucell.nat;
        ComputedFrame published;

        const std::vector<double> reference_cell = ipi_cell_bohr_from_unitcell(ucell);
        bool checked_initial_positions = false;

        while (true)
        {
            std::string header;
            io_failed = 0;
            io_message.clear();
            if (is_root())
            {
                try
                {
                    header = socket.read_header();
                }
                catch (const IpiSocketClosed&)
                {
                    if (state == DriverState::HasData)
                    {
                        io_failed = 1;
                        io_message = "i-PI peer closed while a computed frame was pending";
                    }
                    else
                    {
                        header.clear();
                    }
                }
                catch (const std::exception& exc)
                {
                    io_failed = 1;
                    io_message = exc.what();
                }
            }
            quit_if_root_io_failed(io_failed, io_message);
            header = bcast_header(header);

            if (header.empty())
            {
                if (is_root())
                {
                    ofs_running << " ABACUS socket driver exiting after peer closed connection" << std::endl;
                }
                break;
            }
            else if (header == "STATUS")
            {
                io_failed = 0;
                io_message.clear();
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
                quit_if_root_io_failed(io_failed, io_message);
            }
            else if (header == "INIT")
            {
                std::int32_t rid = 0;
                std::int32_t nbytes = 0;
                std::string params;
                io_failed = 0;
                io_message.clear();
                if (is_root())
                {
                    if (state != DriverState::NeedInit)
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
                            else if (nbytes > MAX_INIT_BYTES)
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
                quit_if_root_io_failed(io_failed, io_message);
                bcast_socket_int32(rid);
                bcast_socket_int32(nbytes);
                if (nbytes > 0 && is_root())
                {
                    ofs_running << " ABACUS socket INIT params bytes " << nbytes << std::endl;
                }
                state = DriverState::Ready;
                if (is_root())
                {
                    ofs_running << " ABACUS socket INIT replica " << rid << std::endl;
                }
            }
            else if (header == "POSDATA")
            {
                SocketFrame::Matrix9 cell = {{0.0}};
                SocketFrame::Matrix9 inv_cell = {{0.0}};
                std::int32_t nat_socket = 0;
                std::vector<double> positions;
                io_failed = 0;
                io_message.clear();
                if (is_root())
                {
                    if (state != DriverState::Ready)
                    {
                        io_failed = 1;
                        io_message = "POSDATA requires READY state";
                    }
                    else
                    {
                        try
                        {
                            const std::vector<double> cell_values = socket.read_doubles(9);
                            const std::vector<double> inverse_values = socket.read_doubles(9);
                            std::copy(cell_values.begin(), cell_values.end(), cell.begin());
                            std::copy(inverse_values.begin(), inverse_values.end(), inv_cell.begin());
                            nat_socket = socket.read_int32();
                            SocketFrame::CellValidation validation
                                = SocketFrame::validate_ipi_cell(cell,
                                                                 inv_cell,
                                                                 MAX_CELL_CONDITION,
                                                                 INVERSE_ABSOLUTE_TOLERANCE,
                                                                 INVERSE_RELATIVE_TOLERANCE);
                            if (!validation.ok)
                            {
                                io_failed = 1;
                                io_message = "invalid POSDATA cell: " + validation.message;
                            }
                            std::size_t coordinate_count = 0;
                            if (io_failed == 0
                                && !SocketFrame::checked_position_count(nat_socket,
                                                                        ucell.nat,
                                                                        coordinate_count,
                                                                        io_message))
                            {
                                io_failed = 1;
                            }
                            if (io_failed == 0)
                            {
                                positions = socket.read_doubles(coordinate_count);
                                if (!SocketFrame::validate_positions(positions,
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
                quit_if_root_io_failed(io_failed, io_message);
                bcast_socket_int32(nat_socket);
                std::vector<double> cell_values(cell.begin(), cell.end());
                std::vector<double> inverse_values(inv_cell.begin(), inv_cell.end());
                bcast_double_vector(cell_values);
                bcast_double_vector(inverse_values);
                if (!is_root())
                {
                    cell = {{0.0}};
                    inv_cell = {{0.0}};
                    std::copy(cell_values.begin(), cell_values.end(), cell.begin());
                    std::copy(inverse_values.begin(), inverse_values.end(), inv_cell.begin());
                    if (nat_socket >= 0)
                    {
                        positions.assign(static_cast<std::size_t>(3 * nat_socket), 0.0);
                    }
                }
                bcast_double_vector(positions);

                const double max_cell_delta_bohr = max_abs_delta(std::vector<double>(cell.begin(), cell.end()), reference_cell);
                if (max_cell_delta_bohr > unchanged_cell_tolerance(cell))
                {
                    ModuleBase::WARNING_QUIT("ABACUS socket", "variable-cell socket updates are not supported yet.");
                }
                if (!checked_initial_positions)
                {
                    checked_initial_positions = true;
                    if (max_wrapped_direct_delta_from_unitcell(ucell, positions) > 1.0e-5 && is_root())
                    {
                        ModuleBase::WARNING(
                            "ABACUS socket",
                            "first POSDATA positions are not PBC-equivalent to STRU atom order; "
                            "i-PI POSDATA carries no species, so the client atoms should use the same atom order as STRU.");
                    }
                }

                try
                {
                    set_positions_from_ipi_bohr(ucell, positions);
                }
                catch (const std::exception& exc)
                {
                    fail_during_collective_stage("set_positions", exc.what());
                }
                catch (...)
                {
                    fail_during_collective_stage("set_positions",
                                                 "unknown socket position update failure");
                }
                try
                {
                    p_esolver->runner(ucell, istep);
                }
                catch (const std::exception& exc)
                {
                    fail_during_collective_stage("runner", exc.what());
                }
                catch (...)
                {
                    fail_during_collective_stage("runner",
                                                 "unknown socket runner failure");
                }
                ComputedFrame computed;
                computed.scf_converged = all_ranks_converged(p_esolver->conv_esolver);
                if (!computed.scf_converged && is_root())
                {
                    ModuleBase::WARNING(
                        "ABACUS socket",
                        "SCF did not converge; returning the available frame and marking it in i-PI extras.");
                }
                double energy_ry = 0.0;
                try
                {
                    energy_ry = p_esolver->cal_energy();
                }
                catch (const std::exception& exc)
                {
                    fail_during_collective_stage("cal_energy", exc.what());
                }
                catch (...)
                {
                    fail_during_collective_stage("cal_energy",
                                                 "unknown socket energy failure");
                }
                int local_failed = std::isfinite(energy_ry) ? 0 : 1;
                throw_if_any_rank_failed(local_failed,
                                         local_failed == 0 ? "" : "socket energy is not finite");
                if (!std::isfinite(energy_ry))
                {
                    ModuleBase::WARNING_QUIT("ABACUS socket", "socket energy is not finite.");
                }
                computed.energy_hartree = energy_ry * RY_TO_HARTREE;
                if (is_root())
                {
                    ofs_running << " ABACUS socket return energy "
                                << energy_ry << " Ry, "
                                << energy_ry * ModuleBase::Ry_to_eV << " eV, "
                                << computed.energy_hartree << " Ha" << std::endl;
                }
                ModuleBase::matrix force;
                if (inp.cal_force)
                {
                    try
                    {
                        p_esolver->cal_force(ucell, force);
                    }
                    catch (const std::exception& exc)
                    {
                        fail_during_collective_stage("cal_force", exc.what());
                    }
                    catch (...)
                    {
                        fail_during_collective_stage("cal_force",
                                                     "unknown socket force failure");
                    }
                    local_failed = 0;
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
                if (inp.cal_stress)
                {
                    ModuleBase::matrix stress;
                    try
                    {
                        p_esolver->cal_stress(ucell, stress);
                    }
                    catch (const std::exception& exc)
                    {
                        fail_during_collective_stage("cal_stress", exc.what());
                    }
                    catch (...)
                    {
                        fail_during_collective_stage("cal_stress",
                                                     "unknown socket stress failure");
                    }
                    local_failed = 0;
                    std::string local_message;
                    try
                    {
                        const SocketFrame::VirialConversion virial
                            = SocketFrame::make_ipi_virial(matrix9_from_stress(stress),
                                                           ucell.omega,
                                                           STRESS_ABSOLUTE_TOLERANCE,
                                                           STRESS_RELATIVE_TOLERANCE);
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
                computed.valid = true;
                published = computed;
                ++istep;
                state = DriverState::HasData;
            }
            else if (header == "GETFORCE")
            {
                io_failed = 0;
                io_message.clear();
                if (is_root())
                {
                    try
                    {
                        if (state != DriverState::HasData || !published.valid)
                        {
                            throw std::runtime_error("GETFORCE requires HAVEDATA state and a valid frame");
                        }
                        socket.write_header("FORCEREADY");
                        socket.write_double(published.energy_hartree);
                        socket.write_int32(static_cast<std::int32_t>(nat_return));
                        const std::vector<double> forces
                            = published.forces_present
                                  ? published.forces_hartree_per_bohr
                                  : std::vector<double>(static_cast<std::size_t>(3 * nat_return), 0.0);
                        socket.write_doubles(forces);
                        socket.write_doubles(vector_from_matrix9(published.virial_wire_hartree));
                        const std::string extra = properties_extra(published);
                        if (extra.size() > static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max()))
                        {
                            throw std::overflow_error("i-PI extras payload is larger than int32");
                        }
                        socket.write_int32(static_cast<std::int32_t>(extra.size()));
                        socket.write_string(extra);
                    }
                    catch (const std::exception& exc)
                    {
                        io_failed = 1;
                        io_message = exc.what();
                    }
                }
                quit_if_root_io_failed(io_failed, io_message);
                published = ComputedFrame();
                state = DriverState::Ready;
            }
            else if (header == "EXIT")
            {
                if (is_root())
                {
                    ofs_running << " ABACUS socket driver received i-PI EXIT" << std::endl;
                }
                break;
            }
            else
            {
                if (is_root())
                {
                    io_failed = 1;
                    io_message = "unknown i-PI header: " + header;
                }
                quit_if_root_io_failed(io_failed, io_message);
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
