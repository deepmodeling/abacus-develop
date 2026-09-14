#include "socket_utils.h"

#include "source_base/global_function.h"
#include "source_base/mathzone.h"
#include "source_base/parallel_common.h"
#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include "source_cell/unitcell.h"
#include "source_cell/update_cell.h"

#ifdef __MPI
#include <mpi.h>
#endif

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <sstream>
#include <stdexcept>

namespace SocketUtils
{
bool all_ranks_converged(const bool local_converged)
{
    int converged = local_converged ? 1 : 0;
    Parallel_Reduce::reduce_min(converged);
    return converged != 0;
}

void throw_if_any_rank_failed(int local_failed, std::string local_message)
{
    int any_failed = local_failed;
    Parallel_Reduce::reduce_max(any_failed);
    if (any_failed != 0)
    {
        if (local_message.empty())
        {
            local_message = "socket frame validation failed on another MPI rank";
        }
        throw std::runtime_error(local_message);
    }
}

/// @brief Terminate the socket run after a failure inside a collective stage.
///
/// Under MPI the failure may deadlock other ranks waiting in a collective,
/// so print a diagnostic and abort the whole job. In serial builds there is
/// no collective: the exception propagates to Socket_Driver::socket_driver,
/// which reports a clean WARNING_QUIT into the log. Do not replace the
/// serial throw with abort(); serial builds are user-facing.
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
    int rank = kIpiRankRoot;
#ifdef __MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
    return rank == kIpiRankRoot;
}

void bcast_double_vector(std::vector<double>& values)
{
    if (!values.empty())
    {
        Parallel_Common::bcast_double(values.data(), static_cast<int>(values.size()));
    }
}

void bcast_int(int& value)
{
    Parallel_Common::bcast_int(value);
}

void bcast_int32(std::int32_t& value)
{
    int tmp = static_cast<int>(value);
    Parallel_Common::bcast_int(tmp);
    value = static_cast<std::int32_t>(tmp);
}

void bcast_chars(char* value, const int size)
{
    Parallel_Common::bcast_char(value, size);
}

void bcast_string(std::string& value)
{
    int size = static_cast<int>(value.size());
    bcast_int(size);
    if (!is_root())
    {
        value.resize(static_cast<std::size_t>(size));
    }
    if (size > 0)
    {
        bcast_chars(&value[0], size);
    }
}

void quit_if_root_failed(int root_failed, std::string root_message)
{
    bcast_int(root_failed);
    bcast_string(root_message);
    if (root_failed != 0)
    {
        ModuleBase::WARNING_QUIT("ABACUS socket", root_message.empty() ? "i-PI socket I/O failed" : root_message);
    }
}

std::string bcast_header(std::string header)
{
    bcast_string(header);
    return header;
}

std::string address()
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
            out[static_cast<std::size_t>(3 * iat + idir)] = value * kRyToHartree;
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
} // namespace SocketUtils
