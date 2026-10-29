#ifndef LIBRPA_STRU_SYMMETRY_H
#define LIBRPA_STRU_SYMMETRY_H

#include <array>
#include <iosfwd>
#include <string>
#include <vector>

namespace ModuleBase
{
class Matrix3;
template <typename T>
class Vector3;
}

namespace RpaLriDetail
{
struct LibRpaSymmetryOperation
{
    double rotation[9];
    double translation[3];
};

struct LibRpaSpinSymmetryOperation
{
    int antiunitary = 0;
    std::array<double, 8> spin_u{};
};

int checked_near_int(double value, const std::string& context);

LibRpaSymmetryOperation make_librpa_symmetry_operation(const ModuleBase::Matrix3& rotation,
                                                       const ModuleBase::Vector3<double>& translation);

std::vector<LibRpaSpinSymmetryOperation> make_librpa_spin_symmetry_operations(
    const ModuleBase::Matrix3& latvec,
    const std::vector<ModuleBase::Matrix3>& unitary_rotations,
    const std::vector<ModuleBase::Matrix3>& antiunitary_rotations);

void write_librpa_symmetry_rows(std::ostream& output,
                                const std::vector<LibRpaSymmetryOperation>& unitary,
                                const std::vector<LibRpaSymmetryOperation>& antiunitary);

void write_librpa_spin_symmetry(std::ostream& output,
                                int grey_group,
                                int spin_source,
                                const std::vector<LibRpaSpinSymmetryOperation>& operations);
} // namespace RpaLriDetail

#endif
