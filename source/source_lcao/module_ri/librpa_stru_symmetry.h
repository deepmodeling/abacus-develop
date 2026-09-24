#ifndef LIBRPA_STRU_SYMMETRY_H
#define LIBRPA_STRU_SYMMETRY_H

#include "source_base/matrix3.h"
#include "source_base/vector3.h"

#include <cmath>
#include <iomanip>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace RpaLriDetail
{
struct LibRpaSymmetryOperation
{
    double rotation[9];
    double translation[3];
};

inline int checked_near_int(const double value, const std::string& context)
{
    const double rounded = std::round(value);
    if (std::abs(value - rounded) > 1e-8)
    {
        throw std::runtime_error(context + " is not close to an integer.");
    }
    return static_cast<int>(rounded);
}

inline LibRpaSymmetryOperation make_librpa_symmetry_operation(
    const ModuleBase::Matrix3& rotation,
    const ModuleBase::Vector3<double>& translation)
{
    LibRpaSymmetryOperation operation = {{rotation.e11,
                                          rotation.e12,
                                          rotation.e13,
                                          rotation.e21,
                                          rotation.e22,
                                          rotation.e23,
                                          rotation.e31,
                                          rotation.e32,
                                          rotation.e33},
                                         {translation.x, translation.y, translation.z}};
    return operation;
}

inline void write_librpa_symmetry_rows(std::ostream& output,
                                       const std::vector<LibRpaSymmetryOperation>& unitary,
                                       const std::vector<LibRpaSymmetryOperation>& antiunitary)
{
    const auto write_operation = [&output](const LibRpaSymmetryOperation& operation) {
        output << std::setw(4) << checked_near_int(operation.rotation[0], "symmetry rotation e11")
               << std::setw(4) << checked_near_int(operation.rotation[1], "symmetry rotation e12")
               << std::setw(4) << checked_near_int(operation.rotation[2], "symmetry rotation e13")
               << std::setw(4) << checked_near_int(operation.rotation[3], "symmetry rotation e21")
               << std::setw(4) << checked_near_int(operation.rotation[4], "symmetry rotation e22")
               << std::setw(4) << checked_near_int(operation.rotation[5], "symmetry rotation e23")
               << std::setw(4) << checked_near_int(operation.rotation[6], "symmetry rotation e31")
               << std::setw(4) << checked_near_int(operation.rotation[7], "symmetry rotation e32")
               << std::setw(4) << checked_near_int(operation.rotation[8], "symmetry rotation e33")
               << std::setw(24) << std::scientific << std::setprecision(15) << operation.translation[0]
               << std::setw(24) << std::scientific << std::setprecision(15) << operation.translation[1]
               << std::setw(24) << std::scientific << std::setprecision(15) << operation.translation[2]
               << std::endl;
    };

    output << (unitary.size() + antiunitary.size()) << " row" << std::endl;
    for (const auto& operation: unitary)
    {
        write_operation(operation);
    }
    for (const auto& operation: antiunitary)
    {
        write_operation(operation);
    }
}
} // namespace RpaLriDetail

#endif
