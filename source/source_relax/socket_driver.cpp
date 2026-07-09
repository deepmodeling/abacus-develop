#include "socket_driver.h"

#include "source_relax/ipi_socket.h"
#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "source_base/mathzone.h"
#include "source_base/parallel_common.h"
#include "source_base/timer.h"
#include "source_cell/update_cell.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <string>
#include <vector>

namespace
{
constexpr double RY_TO_HARTREE = 0.5;
constexpr int IPI_RANK_ROOT = 0;

bool is_root()
{
    return GlobalV::MY_RANK == IPI_RANK_ROOT;
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

std::string ipi_address()
{
    const char* env = std::getenv("ABACUS_IPI_ADDRESS");
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

std::vector<double> flatten_forces_hartree_per_bohr(const ModuleBase::matrix& force)
{
    std::vector<double> out(static_cast<std::size_t>(force.nr * force.nc));
    for (int iat = 0; iat < force.nr; ++iat)
    {
        for (int idir = 0; idir < force.nc; ++idir)
        {
            out[static_cast<std::size_t>(3 * iat + idir)] = force(iat, idir) * RY_TO_HARTREE;
        }
    }
    return out;
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
    if (!inp.cal_force)
    {
        ModuleBase::WARNING_QUIT("ABACUS socket", "socket calculation requires cal_force=1 for i-PI GETFORCE.");
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
                const std::string address = ipi_address();
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

        bool isinit = false;
        bool hasdata = false;
        int istep = 0;
        const int nat_return = ucell.nat;
        double energy_hartree = 0.0;
        std::vector<double> forces_hartree_bohr(static_cast<std::size_t>(3 * ucell.nat), 0.0);
        std::vector<double> virial_hartree(9, 0.0);

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
                    header.clear();
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
                        if (hasdata)
                        {
                            socket.write_header("HAVEDATA");
                        }
                        else if (isinit)
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
                int rid = 0;
                int nbytes = 0;
                std::string params;
                io_failed = 0;
                io_message.clear();
                if (is_root())
                {
                    try
                    {
                        rid = socket.read_int();
                        nbytes = socket.read_int();
                        if (nbytes < 0)
                        {
                            io_failed = 1;
                            io_message = "negative INIT payload length from i-PI socket";
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
                quit_if_root_io_failed(io_failed, io_message);
                bcast_socket_int(rid);
                bcast_socket_int(nbytes);
                if (nbytes > 0 && is_root())
                {
                    ofs_running << " ABACUS socket INIT params bytes " << nbytes << std::endl;
                }
                isinit = true;
                if (is_root())
                {
                    ofs_running << " ABACUS socket INIT replica " << rid << std::endl;
                }
            }
            else if (header == "POSDATA")
            {
                std::vector<double> cell(9, 0.0);
                std::vector<double> inv_cell(9, 0.0);
                int nat_socket = 0;
                std::vector<double> positions;
                io_failed = 0;
                io_message.clear();
                if (is_root())
                {
                    try
                    {
                        cell = socket.read_doubles(9);
                        inv_cell = socket.read_doubles(9);
                        nat_socket = socket.read_int();
                        if (nat_socket < 0)
                        {
                            io_failed = 1;
                            io_message = "negative POSDATA atom count from i-PI socket";
                        }
                        else
                        {
                            positions = socket.read_doubles(static_cast<std::size_t>(3 * nat_socket));
                        }
                    }
                    catch (const std::exception& exc)
                    {
                        io_failed = 1;
                        io_message = exc.what();
                    }
                }
                quit_if_root_io_failed(io_failed, io_message);
                bcast_double_vector(cell);
                bcast_double_vector(inv_cell);
                bcast_socket_int(nat_socket);
                if (!is_root())
                {
                    positions.assign(static_cast<std::size_t>(3 * nat_socket), 0.0);
                }
                bcast_double_vector(positions);

                if (nat_socket != ucell.nat)
                {
                    ModuleBase::WARNING_QUIT("ABACUS socket", "POSDATA atom count does not match STRU.");
                }
                const double max_cell_delta_bohr = max_abs_delta(cell, reference_cell);
                if (max_cell_delta_bohr > 1.0e-6)
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

                set_positions_from_ipi_bohr(ucell, positions);
                p_esolver->runner(ucell, istep);
                const double energy_ry = p_esolver->cal_energy();
                energy_hartree = energy_ry * RY_TO_HARTREE;
                if (is_root())
                {
                    ofs_running << " ABACUS socket return energy "
                                << energy_ry << " Ry, "
                                << energy_ry * ModuleBase::Ry_to_eV << " eV, "
                                << energy_hartree << " Ha" << std::endl;
                }
                ModuleBase::matrix force;
                if (inp.cal_force)
                {
                    p_esolver->cal_force(ucell, force);
                    forces_hartree_bohr = flatten_forces_hartree_per_bohr(force);
                }
                ++istep;
                hasdata = true;
            }
            else if (header == "GETFORCE")
            {
                io_failed = 0;
                io_message.clear();
                if (is_root())
                {
                    try
                    {
                        socket.write_header("FORCEREADY");
                        socket.write_double(energy_hartree);
                        socket.write_int(nat_return);
                        socket.write_doubles(forces_hartree_bohr);
                        socket.write_doubles(virial_hartree);
                        socket.write_int(0);
                    }
                    catch (const std::exception& exc)
                    {
                        io_failed = 1;
                        io_message = exc.what();
                    }
                }
                quit_if_root_io_failed(io_failed, io_message);
                isinit = false;
                hasdata = false;
            }
            else
            {
                if (is_root())
                {
                    ofs_running << " ABACUS socket driver exiting on header " << header << std::endl;
                }
                break;
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
