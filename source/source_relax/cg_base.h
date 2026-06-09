#ifndef CG_BASE_H
#define CG_BASE_H

/**
 * @class CG_Base
 * @brief Base class for conjugate gradient optimization algorithms.
 * 
 * This class provides common CG utility functions that can be shared between
 * different optimization implementations (ions movement and lattice change).
 * The functions are designed to be independent of specific data sources.
 */
class CG_Base {
public:
    /**
     * @brief Setup conjugate gradient direction.
     * @param dim Dimension of the optimization problem
     * @param grad Current gradient
     * @param grad0 Previous gradient
     * @param cg_grad Output CG direction
     * @param cg_grad0 Previous CG direction
     * @param ncggrad CG iteration counter
     * @param flag Restart flag
     */
    void setup_cg_grad(const int dim, double *grad, const double *grad0,
                       double *cg_grad, const double *cg_grad0,
                       const int &ncggrad, int &flag);

    /**
     * @brief Brent's method for one-dimensional minimization.
     * @param fa Function value at xa
     * @param fb Function value at xb
     * @param fc Function value at xc
     * @param xa Left bound
     * @param xb Middle point
     * @param xc Right bound
     * @param best_x Output best step length
     * @param xpt Previous best point
     */
    void Brent(double &fa, double &fb, double &fc,
               double &xa, double &xb, double &xc,
               double &best_x, double &xpt);

    /**
     * @brief Third-order polynomial interpolation for line search.
     * @param e0 Energy at previous step
     * @param e1 Energy at current step
     * @param fa Gradient projection at previous step
     * @param fb Gradient projection at current step
     * @param x Current step length
     * @param best_x Output optimal step length
     */
    void third_order(const double &e0, const double &e1,
                     const double &fa, const double &fb,
                     const double x, double &best_x);

    /**
     * @brief Calculate projection of gradient onto search direction.
     * @param dim Dimension of vectors
     * @param g0 First vector (usually search direction)
     * @param g1 Second vector (usually gradient)
     * @param f_value Output projection value
     */
    void f_cal(const int dim, const double *g0, const double *g1, double &f_value);

    /**
     * @brief Calculate movement vector from CG direction.
     * @param dim Dimension of the optimization problem
     * @param move Output movement vector
     * @param cg_gradn Normalized CG direction
     * @param trust_radius Step size limit
     */
    void setup_move(const int dim, double *move, double *cg_gradn, const double &trust_radius);

protected:

    /**
     * @brief Normalize vector.
     * @param dim Dimension of vector
     * @param cg_gradn Output normalized vector
     * @param cg_grad Input vector
     */
    void normalize(const int dim, double *cg_gradn, const double *cg_grad);
};

#endif