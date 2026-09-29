#ifndef PCC_BOUNDARY_H
#define PCC_BOUNDARY_H

#include <string>

namespace ModuleSccs
{

// Electrostatic boundary of the cell: fully periodic, or periodic with the
// point-counter-charge (PCC) open-boundary correction.
enum class Boundary
{
    Periodic,
    Pcc0d,
    Pcc2d
};

// Parse periodic, pcc_0d or pcc_2d (case-insensitive; '-' may replace '_').
Boundary parse_boundary(const std::string& value);

} // namespace ModuleSccs

#endif
