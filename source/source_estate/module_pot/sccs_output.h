#ifndef SCCS_OUTPUT_H
#define SCCS_OUTPUT_H

#include "source_estate/charge_moments.h"
#include "source_hamilt/module_sccs/sccs_coulomb.h"
#include <iosfwd>
#include <vector>

class UnitCell;
namespace ModulePW { class PW_Basis; }
namespace unitcell { struct AtomData; }
namespace ModuleSccs
{
struct SccsConfig;
struct SccsResponse;
struct FunctionalResult;
}
namespace elecstate
{
// Scalar reporting snapshot only: no retained grid buffers or model controls.
struct PccOutput
{
    bool slab = false;
    int axis = 2;
    double coordinate = 0.0;
    ModuleBase::Vector3<double> origin;
    ModuleBase::Vector3<double> normal;
    ChargeMoments moments[4]; // smooth, point, polarization, screened
    double smooth_energy = 0.0; // Hartree
};
struct SccsOutput
{
    bool valid = false;
    bool reused_fixed_sources = false;
    int iterations = 0;
    bool warm_started = false;
    double residual_rms = 0.0;
    double residual_max = 0.0;
    double reaction_energy = 0.0;
    double volume = 0.0;
    double surface = 0.0;
    bool filled = false; // solvent-aware cavity
    double filled_volume = 0.0;
    double far_field_charge = 0.0;
    double expected_charge = 0.0;
    ModuleSccs::CoulombTransformCounts transforms;
    PccOutput pcc;
};
// All ranks in the PW pool call together.
void collect_output(const UnitCell& cell,
                     const ModulePW::PW_Basis& basis,
                     const std::vector<unitcell::AtomData>& atoms,
                     const std::vector<ModuleBase::Vector3<double>>& grid_positions,
                     const std::vector<double>& ions,
                     const std::vector<double>& charge,
                     const ModuleSccs::SccsConfig& config,
                     const ModuleSccs::SccsResponse& response,
                     SccsOutput& result);
// sccs_debug lines of one SCF iteration, tagged and column-aligned.
void write_sccs_out(std::ostream& output, const SccsOutput& result, int level,
                    double solvation_energy, bool pcc, double pcc_energy);
// The line of an SCF iteration while a delayed SCCS waits.
void write_sccs_wait(std::ostream& output, double residual, double start_drho, int start_nmax);
} // namespace elecstate
#endif
