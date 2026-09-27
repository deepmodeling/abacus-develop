#ifndef ESOLVER_DFTB_NATIVE_INPUT_H
#define ESOLVER_DFTB_NATIVE_INPUT_H

#include "source_dftb/periodic_scc.h"

#include <map>
#include <string>
#include <vector>

// Internal declarations shared by native DFTB solver translation units.
class UnitCell;

namespace ModuleESolver
{
namespace NativeDftb
{
struct NativeDftbConfig
{
    std::string skf_directory;
    std::string band_path_file;
    std::map<std::string, double> hubbard_derivatives;
    double temperature_kelvin = 0.0;
    double scc_tolerance = 1.0e-6;
    int maximum_scc_iterations = 200;
    double mixing_parameter = 0.2;
    std::string mixing_method = "linear";
    int mixing_history = 6;
    double broyden_inverse_jacobi_weight = 0.01;
    double broyden_minimal_weight = 1.0;
    double broyden_maximal_weight = 1.0e5;
    double broyden_weight_factor = 1.0e-2;
    int output_precision = 12;
    bool third_order = false;
};

NativeDftbConfig read_native_config(const std::string& filename);
std::vector<ModuleDFTB::DftbWeightedKPoint> read_native_kpoints(const std::string& filename, const UnitCell& ucell);
std::vector<ModuleDFTB::DftbBandKPoint> read_native_band_path(const std::string& filename);
} // namespace NativeDftb
} // namespace ModuleESolver

#endif // ESOLVER_DFTB_NATIVE_INPUT_H
