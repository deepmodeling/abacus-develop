#ifndef SOURCE_RELAX_SOCKET_FRAME_H
#define SOURCE_RELAX_SOCKET_FRAME_H

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace SocketFrame
{
using Matrix9 = std::array<double, 9>;

struct CellValidation
{
    bool ok;
    std::string message;
    double determinant_bohr3;
    double condition_number_2;
    double inverse_residual;
    Matrix9 computed_inverse_wire_bohr_inv;
};

struct VirialConversion
{
    bool ok;
    std::string message;
    Matrix9 wire_virial_hartree;
    double max_antisymmetric_component;
};

Matrix9 transpose_matrix9(const Matrix9& values);
CellValidation validate_ipi_cell(const Matrix9& cell_wire,
                                 const Matrix9& inverse_wire,
                                 double max_condition_number,
                                 double inverse_absolute_tolerance,
                                 double inverse_relative_tolerance);
bool validate_positions(const std::vector<double>& positions_bohr,
                        std::size_t coordinate_count,
                        std::string& message);
bool checked_position_count(std::int32_t nat_socket,
                            int nat_expected,
                            std::size_t& coordinate_count,
                            std::string& message);
VirialConversion make_ipi_virial(const Matrix9& stress_ry_per_bohr3,
                                  double volume_bohr3,
                                  double antisymmetric_absolute_tolerance,
                                  double antisymmetric_relative_tolerance);
} // namespace SocketFrame

#endif
