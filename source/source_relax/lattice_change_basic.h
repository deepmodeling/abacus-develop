#ifndef LATTICE_CHANGE_BASIC_H
#define LATTICE_CHANGE_BASIC_H

#include <fstream>
#include <vector>
#include "relax_data.h"
#include "source_base/matrix.h"
#include "source_cell/unitcell.h"

/**
 * @namespace Lattice_Change_Basic
 * @brief Basic utilities and shared state for lattice relaxation algorithms.
 *
 * This namespace provides common functions and parameters used by lattice
 * optimization methods. It shares core state variables through references to
 * Relax_Data, ensuring consistent data across different optimization algorithms.
 *
 */
namespace Lattice_Change_Basic
{
// Shared state variables (referenced from Relax_Data for unified data sharing)
static int& dim = Relax_Data::dim_lattice;      ///< Dimension of free variables (9 for full lattice)
static double& largest_grad = Relax_Data::largest_grad; ///< Largest gradient component

// Lattice-specific parameters (not shared with ion movement)
extern int update_iter;          ///< Number of successfully updated iterations
extern int stress_step;          ///< Index of stress optimization step
extern double lattice_change_ini; ///< Initial trust radius for lattice change (default: 0.01)
extern std::string fixed_axes;    ///< Fixed axes constraint ("None", "shape", "volume", or specific axes)

/**
 * @brief Setup gradient from stress tensor.
 * @param ucell Unit cell containing lattice information
 * @param lat Output lattice vector array (9 elements: e11, e12, e13, e21, e22, e23, e31, e32, e33)
 * @param grad Output gradient array (9 elements)
 * @param stress Stress tensor (3x3 matrix)
 */
void setup_gradient(const UnitCell &ucell, double *lat, double *grad, ModuleBase::matrix &stress);

/**
 * @brief Update lattice vectors according to displacement.
 * @param ucell Unit cell to update
 * @param move Displacement vector for lattice change (9 elements)
 * @param lat Current lattice vectors (9 elements)
 */
void change_lattice(UnitCell &ucell, double *move, double *lat);

/**
 * @brief Check convergence based on stress threshold.
 * @param ucell Unit cell containing lattice constraints
 * @param stress Stress tensor (3x3 matrix)
 * @param grad Gradient array (9 elements)
 * @param ofs Output stream for logging
 * @return true if converged, false otherwise
 */
bool check_converged(const UnitCell &ucell, ModuleBase::matrix &stress, double *grad, std::ofstream& ofs);

/**
 * @brief Terminate lattice optimization and output results.
 * @param converged Convergence flag
 * @param ofs Output stream for logging
 */
void terminate(const bool converged, std::ofstream& ofs);

/**
 * @brief Update energy values and compute energy difference.
 * @param energy_in Input energy value
 * @param judgement Flag for SD method (true) or BFGS (false)
 * @param etot_info Vector containing [etot, etot_p]
 */
void setup_etot(const double &energy_in, std::vector<double>& etot_info);
} // namespace Lattice_Change_Basic
#endif
