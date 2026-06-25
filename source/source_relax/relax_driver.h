#ifndef RELAX_DRIVER_H
#define RELAX_DRIVER_H

#include "source_cell/unitcell.h"
#include "source_esolver/esolver.h"
#include "relax_sync.h"
#include "relax_nsync.h"
#include "source_io/module_parameter/input_parameter.h"
#include "source_base/matrix.h"
#include <vector>
#include <fstream>

/**
 * @brief Driver class for geometry relaxation calculations.
 *
 * This class manages the main relaxation loop, including energy solving,
 * force/stress calculation, and structure optimization steps.
 */
class Relax_Driver
{

  public:
    Relax_Driver(){};
    ~Relax_Driver(){};

    /**
     * @brief Main driver function for relaxation calculations.
     *
     * This function executes the main iteration loop for relaxation,
     * calling energy solver and structure optimization steps until
     * convergence or maximum steps reached.
     *
     * @param p_esolver Pointer to the energy solver.
     * @param ucell Reference to the unit cell to be relaxed.
     * @param inp Input parameters for the calculation.
     * @param ofs_running Output stream for running log.
     */
    void relax_driver(ModuleESolver::ESolver* p_esolver,
            UnitCell& ucell,
            const Input_para& inp,
            std::ofstream& ofs_running);

  private:
    /// New relaxation optimizer (Relax class)
    Relax rl;
    /// Old relaxation optimizer (IonCellOptimizer class)
    IonCellOptimizer rl_old;

    /**
     * @brief Initialize the relaxation optimizer.
     *
     * @param nat Number of atoms in the unit cell.
     * @param inp Input parameters for the calculation.
     */
    void init_relax(const int nat, const Input_para& inp);

    /**
     * @brief Print iteration information to screen.
     *
     * @param steps Vector containing step counters: steps[0]=istep, steps[1]=force_step, steps[2]=stress_step.
     * @param inp Input parameters for the calculation.
     */
    void iter_info(const std::vector<int>& steps, const Input_para& inp);

    /**
     * @brief Perform energy solving for the current step.
     *
     * @param istep Current iteration step.
     * @param p_esolver Pointer to the energy solver.
     * @param ucell Reference to the unit cell.
     * @param inp Input parameters for the calculation.
     * @param force Output matrix for calculated forces.
     * @param stress Output matrix for calculated stress.
     * @param etot Output total energy.
     */
    void esolve(const int istep, ModuleESolver::ESolver* p_esolver, UnitCell& ucell,
            const Input_para& inp, ModuleBase::matrix& force, ModuleBase::matrix& stress, double& etot);

    /**
     * @brief Perform one relaxation step.
     *
     * @param steps Vector containing step counters (modified in place): steps[0]=istep, steps[1]=force_step, steps[2]=stress_step.
     * @param p_esolver Pointer to the energy solver.
     * @param ucell Reference to the unit cell.
     * @param inp Input parameters for the calculation.
     * @param force Matrix of calculated forces.
     * @param stress Matrix of calculated stress.
     * @param etot Total energy.
     * @param ofs_running Output stream for running log.
     * @return True if relaxation converged, false otherwise.
     */
    bool relax_step(std::vector<int>& steps, ModuleESolver::ESolver* p_esolver, UnitCell& ucell,
            const Input_para& inp, const ModuleBase::matrix& force, const ModuleBase::matrix& stress,
            const double etot, std::ofstream& ofs_running);

    /**
     * @brief Output structure files after relaxation step.
     *
     * @param istep Current iteration step.
     * @param ucell Reference to the unit cell.
     * @param inp Input parameters for the calculation.
     */
    void stru_out(const int istep, UnitCell& ucell, const Input_para& inp);

    /**
     * @brief Output JSON format results.
     *
     * @param p_esolver Pointer to the energy solver.
     * @param ucell Reference to the unit cell.
     * @param inp Input parameters for the calculation.
     * @param force Matrix of calculated forces.
     * @param stress Matrix of calculated stress.
     */
    void json_out(ModuleESolver::ESolver* p_esolver, UnitCell& ucell, const Input_para& inp,
            const ModuleBase::matrix& force, const ModuleBase::matrix& stress);

    /**
     * @brief Output final results after relaxation.
     *
     * @param istep Final iteration step.
     * @param ucell Reference to the unit cell.
     * @param inp Input parameters for the calculation.
     */
    void final_out(const int istep, UnitCell& ucell, const Input_para& inp);
};

#endif
