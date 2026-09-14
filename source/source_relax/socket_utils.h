#ifndef ABACUS_SOURCE_RELAX_SOCKET_UTILS_H
#define ABACUS_SOURCE_RELAX_SOCKET_UTILS_H

#include "source_relax/socket_frame.h"
#include "source_base/matrix.h"

#include <cstdint>
#include <limits>
#include <string>
#include <vector>

class UnitCell;

namespace SocketUtils
{
constexpr double kRyToHartree = 0.5;
constexpr int kIpiRankRoot = 0;
constexpr double kMaxCellCondition = 1.0e12;
constexpr double kInverseAbsoluteTolerance
    = 64.0 * std::numeric_limits<double>::epsilon();
constexpr double kInverseRelativeTolerance = 64.0;
constexpr double kStressAbsoluteTolerance = 1.0e-10;
constexpr double kStressRelativeTolerance = 1.0e-8;
constexpr std::int32_t kMaxInitBytes = static_cast<std::int32_t>(1048576);

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

bool all_ranks_converged(bool local_converged);
void throw_if_any_rank_failed(int local_failed, std::string local_message);
[[noreturn]] void fail_during_collective_stage(const char* stage,
                                               const std::string& message);
std::string properties_extra(const ComputedFrame& frame);
bool is_root();
void bcast_double_vector(std::vector<double>& values);
void bcast_int(int& value);
void bcast_int32(std::int32_t& value);
void bcast_chars(char* value, int size);
void bcast_string(std::string& value);
void quit_if_root_failed(int root_failed, std::string root_message);
std::string bcast_header(std::string header);
std::string address();
std::vector<double> ipi_cell_bohr_from_unitcell(const UnitCell& ucell);
double max_wrapped_direct_delta_from_unitcell(const UnitCell& ucell,
                                              const std::vector<double>& positions_bohr);
double max_abs_delta(const std::vector<double>& a, const std::vector<double>& b);
double unchanged_cell_tolerance(const SocketFrame::Matrix9& cell);
void set_positions_from_ipi_bohr(UnitCell& ucell,
                                 const std::vector<double>& positions_bohr);
std::vector<double> flatten_forces_hartree_per_bohr(const ModuleBase::matrix& force,
                                                    int nat);
SocketFrame::Matrix9 matrix9_from_stress(const ModuleBase::matrix& stress);
std::vector<double> vector_from_matrix9(const SocketFrame::Matrix9& values);
} // namespace SocketUtils

#endif
