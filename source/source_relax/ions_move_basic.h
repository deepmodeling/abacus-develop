#ifndef IONS_MOVE_BASIC_H
#define IONS_MOVE_BASIC_H

#include <fstream>
#include <iostream>
#include "relax_data.h"
#include "source_base/matrix.h"
#include "source_cell/unitcell.h"

/**
 * @namespace Ions_Move_Basic
 * @brief Basic utilities and shared state for ionic relaxation algorithms.
 * 
 * This namespace provides common functions and parameters used by all
 * ion movement methods (BFGS, CG, SD, etc.). It shares core state variables
 * through references to Relax_Data, ensuring consistent data across different
 * optimization algorithms.
 */
namespace Ions_Move_Basic
{
static int& dim = Relax_Data::dim;              ///< Dimension of free variables (3 * number of atoms)
static double& largest_grad = Relax_Data::largest_grad; ///< Largest gradient component

// Ions-specific parameters (not shared with lattice change)
extern double trust_radius;          ///< Current trust radius
extern double trust_radius_old;      ///< Previous trust radius
extern double relax_bfgs_rmax;       ///< Maximum trust radius (default: 0.8 Bohr)
extern double relax_bfgs_rmin;       ///< Minimum trust radius (default: 1e-5 Bohr)
extern double relax_bfgs_init;       ///< Initial trust radius (default: 0.5 Bohr)
extern double best_xxx;              ///< Last step length from CG, used as BFGS initial guess

/**
 * @brief Setup gradient from atomic forces.
 * @param ucell Unit cell containing atomic information
 * @param force Force matrix (nat x 3)
 * @param pos Output position array (dimension: dim)
 * @param grad Output gradient array (dimension: dim)
 * @param ofs Output stream for logging
 */
void setup_gradient(const UnitCell &ucell, const ModuleBase::matrix &force, double *pos, double *grad, std::ofstream& ofs);

/**
 * @brief Move atoms according to displacement vector.
 * @param ucell Unit cell to update
 * @param move Displacement vector (dimension: dim)
 * @param pos Current position array (dimension: dim)
 * @param ofs Output stream for logging
 */
void move_atoms(UnitCell &ucell, double *move, double *pos, std::ofstream& ofs);

/**
 * @brief Check convergence based on gradient threshold.
 * @param ucell Unit cell containing lattice information
 * @param grad Gradient array (dimension: dim)
 * @param update_iter Number of successfully updated iterations (will be incremented if converged)
 * @param ofs Output stream for logging
 * @param etot_info Energy information array [etot, etot_p, ediff]
 * @return true if converged, false otherwise
 */
bool check_converged(const UnitCell &ucell, const double *grad, int& update_iter, std::ofstream& ofs, std::vector<double>& etot_info);

/**
 * @brief Terminate geometry optimization and output results.
 * @param converged Convergence flag
 * @param update_iter Number of successfully updated iterations
 * @param ucell Unit cell to output
 * @param istep Current ionic step index
 * @param ofs Output stream for logging
 */
void terminate(const bool converged, const int update_iter, const UnitCell &ucell, const int istep, std::ofstream& ofs);

/**
 * @brief Update energy values and compute energy difference.
 * @param energy_in Input energy value
 * @param judgement Flag for SD method (true) or BFGS (false)
 * @param istep Current ionic step index
 * @param etot_info Energy information array [etot, etot_p, ediff]
 */
void setup_etot(const double &energy_in, const bool judgement, const int istep, std::ofstream& ofs, std::vector<double>& etot_info);

/**
 * @brief Compute dot product of two vectors.
 * @param a First vector
 * @param b Second vector
 * @param dim_in Dimension of vectors
 * @return Dot product value
 */
double dot_func(const double *a, const double *b, const int &dim_in);

/**
 * @brief Third-order polynomial interpolation for line search.
 */
void third_order();

} // namespace Ions_Move_Basic
#endif
