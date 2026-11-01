#ifndef SOL_GRID_FIELD_H
#define SOL_GRID_FIELD_H

#include <string>
#include <vector>

namespace elecstate
{
/// A real-space field of an implicit-solvent model on this rank's density
/// grid, written by out_sol as sol_{name}.cube.
struct SolventGridField
{
    std::string name;
    std::vector<double> values;
};
} // namespace elecstate

#endif
