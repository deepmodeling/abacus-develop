#ifndef LBFGS_H
#define LBFGS_H

#include <vector>
#include <tuple> 
#include <algorithm>
#include <cmath>
#include "matrix_methods.h"
//#include "line_search.h"
#include "source_base/matrix.h"
#include "source_base/matrix3.h"
#include "source_cell/unitcell.h"
#include "source_esolver/esolver.h"
#include "source_esolver/esolver_ks.h"

/**
 * @class LBFGS
 * @brief Implements L-BFGS optimization algorithm for structural relaxation
 */
class LBFGS
{
public:
    /**
     * @brief Initialize L-BFGS parameters
     * @param _size Number of atoms in system
     */
    void allocate(const int _size);

    /**
     * @brief Perform one L-BFGS relaxation step
     * @param _force Current force
     * @param ucell Unit cell to optimize
     * @param etot Current total energy
     * @param p_esolver  Structure solver
     */
    void relax_step(const ModuleBase::matrix _force,
                    UnitCell& ucell,
                    const double &etot);

private:
    //LineSearch l_search;
    double alpha;                           ///< Initial Hessian diagonal element
    double maxstep;                         ///< Maximum allowed step length
    int size;                               ///< Number of atoms in system
    int memory;                             ///< Number of previous steps to store
    double H0;                              ///< Initial inverse Hessian approximation
    int iteration;                          ///< Current iteration count
    double energy;                          ///< Current system energy
    double alpha_k;                         ///< Step size parameter

    ModuleESolver::ESolver* solver;         ///< Structure solver
    std::vector<double> steplength;//the length of atoms displacement 
    std::vector<std::vector<double>> H;//Hessian matrix
    std::vector<double> force0;//force in previous step
    std::vector<ModuleBase::Vector3<double>> force;
    std::vector<double> pos0;//atom pos in previous step(cartesian coordinates)
    std::vector<ModuleBase::Vector3<double>> pos;
    std::vector<double> pos_taud0;//atom pos in previous step(relative coordinates)
    std::vector<ModuleBase::Vector3<double>> pos_taud;
    std::vector<ModuleBase::Vector3<double>> dpos;
    std::vector<std::vector<double>> s;     ///< Position difference vectors
    std::vector<std::vector<double>> y;     ///< Force difference vectors
    std::vector<double> rho;                ///< Scalar products for L-BFGS update

    /**
     * @brief Prepare optimization step parameters
     */
    void prepare_step(std::vector<ModuleBase::Vector3<double>>& force,
                      std::vector<ModuleBase::Vector3<double>>& pos,
                      std::vector<std::vector<double>>& H,
                      std::vector<double>& pos0,
                      std::vector<double>& force0,
                      std::vector<ModuleBase::Vector3<double>>& dpos,
                      UnitCell& ucell,
                      const double &etot);

    /**
     * @brief Judge if the cell is restrain
     * @param dpos Position displacements to constrain
     */
    void is_restrain();

    /**
     * @brief Calculate maximum gradient component
     * @param _force Current force matrix
     * @param ucell Unit cell being optimized
     */
    void calculate_largest_grad(const ModuleBase::matrix& _force,
                                UnitCell& ucell);

    /**
     * @brief Extract atomic positions from unit cell
     * @param ucell Unit cell to read
     * @param pos Output position vector
     */
    void get_pos(UnitCell& ucell,
                 std::vector<ModuleBase::Vector3<double>>& pos);

    /**
     * @brief Get fractional positions from unit cell
     * @param ucell Unit cell to read
     * @param pos_taud Output fractional positions
     */
    void get_pos_taud(UnitCell& ucell,
                      std::vector<ModuleBase::Vector3<double>>& pos_taud);

    /**
     * @brief Update L-BFGS history buffers
     * @param pos_taud Current fractional positions
     * @param pos_taud0 Previous fractional positions
     * @param force Current forces
     * @param force0 Previous forces
     * @param ucell Unit cell being optimized
     * @param iteration Current step number
     * @param memory History buffer size
     * @param s Position differences buffer
     * @param y Force differences buffer
     * @param rho Scalar products buffer
     */
    void update(std::vector<ModuleBase::Vector3<double>>& pos_taud, 
                std::vector<double>& pos_taud0, 
                std::vector<double>& force,
                std::vector<double>& force0, 
                UnitCell& ucell,
                int iteration,
                int memory,
                std::vector<std::vector<double>>& s,
                std::vector<std::vector<double>>& y,
                std::vector<double>& rho);

    /**
     * @brief Determine optimal step lengths
     * @param steplength Output step lengths
     * @param dpos Position displacements
     * @param maxstep Maximum allowed step length
     */
    void determine_step(std::vector<double>& steplength,
                       std::vector<ModuleBase::Vector3<double>>& dpos,
                       double& maxstep);

    /**
     * @brief Update atomic positions in unit cell
     * @param ucell Unit cell to update
     */
    void update_pos(UnitCell& ucell);  
};

#endif