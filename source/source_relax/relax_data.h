#ifndef RELAX_DATA_H
#define RELAX_DATA_H

#include <vector>

/**
 * @brief Unified data structure for geometry optimization algorithms.
 * 
 * This class provides a common data container for all relaxation methods
 * (BFGS, CG, SD, etc.) to share and reuse key variables like positions,
 * gradients, and move vectors. It serves as the single source of truth
 * for optimization state across different algorithms.
 */
class Relax_Data {
public:
    /**
     * @brief Default constructor.
     */
    Relax_Data() = default;

    /**
     * @brief Default destructor.
     */
    ~Relax_Data() = default;

    /**
     * @brief Allocate memory for data vectors.
     * @param dim_in Dimension of the optimization problem (3 * number of atoms).
     */
    void allocate(const int dim_in);

    // Static members - shared global state across all relaxation instances
    static int dim;              ///< Dimension of free variables (3 * number of atoms) for ion movement
    static int dim_lattice;      ///< Dimension of free variables (9) for lattice change

    static double largest_grad;  ///< Largest gradient component (force) in current step

    // Instance members - per-iteration data
    std::vector<double> pos;     ///< Current atomic positions in Bohr
    std::vector<double> grad;    ///< Current gradient (negative force) in Ry/Bohr
    std::vector<double> move;    ///< Displacement vector for next step

    std::vector<double> pos_p;   ///< Previous atomic positions
    std::vector<double> grad_p;  ///< Previous gradient
    std::vector<double> move_p;  ///< Previous displacement
};

#endif