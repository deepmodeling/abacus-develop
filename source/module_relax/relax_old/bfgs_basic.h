#ifndef BFGS_BASIC
#define BFGS_BASIC

#include "module_base/matrix.h"

// references
// 1) Roger Fletcher, Practical Methods of Optimization, John Wiley and
// Sons, Chichester, 2nd edn, 1987.
// 2) Salomon R. Billeter, Alexander J. Turner, Walter Thiel,
// Phys. Chem. Chem. Phys. 2, 2177 (2000).
// 3) Salomon R. Billeter, Alessandro Curioni, Wanda Andreoni,
// Comput. Mat. Science 27, 437, (2003).
// 4) Ren Weiqing, PhD Thesis: Numerical Methods for the Study of Energy
// Landscapes and Rare Events.

class BFGS_Basic
{

  public:
    BFGS_Basic();
    ~BFGS_Basic();

  protected:
    void allocate_basic(void);
    void new_step(const double& lat0);
    void reset_hessian(void);
    void save_bfgs(void);

    double* pos;  // std::vector containing 3N coordinates of the system ( x )
    double* grad; // std::vector containing 3N components of ( grad( V(x) ) )
    double* move; // pos = pos_p + move.

    double* pos_p;  // p: previous
    double* grad_p; // p: previous
    double* move_p;

  public:                        // mohan update 2011-06-12
    static double relax_bfgs_w1; // fixed: parameters for Wolfe conditions.
    static double relax_bfgs_w2; // fixed: parameters for Wolfe conditions.

  protected:
    bool save_flag;
    bool tr_min_hit; //.TRUE. if the trust_radius has already been set
                     // to the minimum value at the previous step

    // mohan add 2010-07-27
    double check_move(const double& lat0, const double& pos, const double& pos_p);

  private:
    bool wolfe_flag;
    ModuleBase::matrix inv_hess;

    int bfgs_ndim;

    void update_inverse_hessian(const double& lat0);
    void check_wolfe_conditions(void);
    void compute_trust_radius(void);
};

#endif
