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
// Shared state variables (referenced from Relax_Data for unified data sharing)
static int& dim = Relax_Data::dim;              ///< Dimension of free variables (3 * number of atoms)
static bool& converged = Relax_Data::converged; ///< Convergence flag
static double& largest_grad = Relax_Data::largest_grad; ///< Largest gradient component
static double& ediff = Relax_Data::ediff;       ///< Energy difference from previous step
static double& etot = Relax_Data::etot;         ///< Total energy of current step
static double& etot_p = Relax_Data::etot_p;     ///< Total energy of previous step

// Ions-specific parameters (not shared with lattice change)
extern int update_iter;              ///< Number of successfully updated iterations
extern double trust_radius;          ///< Current trust radius
extern double trust_radius_old;      ///< Previous trust radius
extern double relax_bfgs_rmax;       ///< Maximum trust radius (default: 0.8 Bohr)
extern double relax_bfgs_rmin;       ///< Minimum trust radius (default: 1e-5 Bohr)
extern double relax_bfgs_init;       ///< Initial trust radius (default: 0.5 Bohr)
extern double best_xxx;              ///< Last step length from CG, used as BFGS initial guess
extern std::vector<std::string> relax_method; ///< Relaxation method settings
extern int out_stru;                 ///< Structure output flag

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
 * @param ofs Output stream for logging
 */
void check_converged(const UnitCell &ucell, const double *grad, std::ofstream& ofs);

/**
 * @brief Terminate geometry optimization and output results.
 * @param ucell Unit cell to output
 * @param istep Current ionic step index
 * @param ofs Output stream for logging
 */
void terminate(const UnitCell &ucell, const int istep, std::ofstream& ofs);

/**
 * @brief Update energy values and compute energy difference.
 * @param energy_in Input energy value
 * @param judgement Flag for SD method (true) or BFGS (false)
 * @param istep Current ionic step index
 */
void setup_etot(const double &energy_in, const bool judgement, const int istep);

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