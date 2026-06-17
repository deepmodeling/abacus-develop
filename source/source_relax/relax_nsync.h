#ifndef ION_CELL_OPTIMIZER_H
#define ION_CELL_OPTIMIZER_H

#include "ions_move_methods.h"
#include "lattice_change_methods.h"
#include "source_cell/unitcell.h"
#include <fstream>

/**
 * @class IonCellOptimizer
 * @brief Optimizer for ionic position relaxation and cell parameter optimization.
 * 
 * This class handles both atomic relaxation and cell relaxation calculations.
 * It manages the geometry optimization process by coordinating between
 * Ions_Move_Methods (for atomic position updates) and Lattice_Change_Methods
 * (for cell parameter updates).
 * 
 * The optimization follows a sequential approach:
1. First perform atomic relaxation until convergence
2. Then perform cell relaxation (only in cell-relax mode)
 */
class IonCellOptimizer
{
  public:
    /**
     * @brief Initialize relaxation algorithms.
     * @param natom Number of atoms in the system
     */
    void init_relax(const int& natom);

    /**
     * @brief Perform one step of relaxation (atomic and/or cell).
     * 
     * This is the main interface for relaxation. Depending on the calculation type
     * ("relax" or "cell-relax"), it will perform atomic relaxation, cell relaxation,
     * or both in sequence.
     * 
     * @param istep Current total iteration step
     * @param energy Total energy of the system
     * @param ucell Unit cell containing atomic positions and lattice vectors
     * @param force Ionic forces matrix
     * @param stress Stress tensor matrix
     * @param force_step Current step counter for force-based relaxation
     * @param stress_step Current step counter for stress-based relaxation
     * @param ofs_running Output stream for running log
     * @return true if relaxation is converged, false otherwise
     */
    bool relax_step(const int& istep,
                    const double& energy,
                    UnitCell& ucell,
                    ModuleBase::matrix force,
                    ModuleBase::matrix stress,
                    int& force_step,
                    int& stress_step,
                    std::ofstream& ofs_running);

  private:
    Ions_Move_Methods IMM;       ///< Ionic movement methods for atom relaxation
    Lattice_Change_Methods LCM;  ///< Lattice change methods for cell relaxation
};

#endif