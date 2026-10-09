#ifndef BFGS_BASIC
#define BFGS_BASIC

#include <fstream>
#include <iostream>
#include "source_base/matrix.h"
#include <vector>

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
    ~BFGS_Basic() = default;

    //====================================================================
    // Test seam.
    //
    // The BFGS state below is protected and the update machinery private,
    // because only Ions_Move_BFGS drives them. The unit tests seed that state
    // and step the algorithm one stage at a time, so each piece they touch is
    // reachable through the accessors and wrappers here rather than by
    // reinterpreting the access specifiers.
    //
    // Production code must keep using the protected/private names directly;
    // nothing outside the tests should call the *_for_testing() wrappers.
    //====================================================================

    /// @brief 3N coordinates of the system ( x )
    std::vector<double>& get_pos()
    {
        return pos;
    }
    /// @brief 3N components of ( grad( V(x) ) )
    std::vector<double>& get_grad()
    {
        return grad;
    }
    /// @brief the step taken, pos = pos_p + move
    std::vector<double>& get_move()
    {
        return move;
    }
    /// @brief coordinates of the previous step
    std::vector<double>& get_pos_p()
    {
        return pos_p;
    }
    /// @brief gradient of the previous step
    std::vector<double>& get_grad_p()
    {
        return grad_p;
    }
    /// @brief step taken at the previous step
    std::vector<double>& get_move_p()
    {
        return move_p;
    }
    /// @brief whether a bfgs state has been saved
    bool& get_save_flag()
    {
        return save_flag;
    }
    /// @brief whether the trust radius already hit its minimum last step
    bool& get_tr_min_hit()
    {
        return tr_min_hit;
    }
    /// @brief whether the Wolfe conditions were satisfied
    bool& get_wolfe_flag()
    {
        return wolfe_flag;
    }
    /// @brief the inverse Hessian of the BFGS update
    ModuleBase::matrix& get_inv_hess()
    {
        return inv_hess;
    }
    /// @brief number of previous steps kept by the BFGS update
    int& get_bfgs_ndim()
    {
        return bfgs_ndim;
    }

    void allocate_basic_for_testing()
    {
        allocate_basic();
    }
    void new_step_for_testing(const double& lat0,
                              int& update_iter,
                              std::ofstream& ofs,
                              std::vector<double>& etot_info,
                              const int test_relax_method)
    {
        new_step(lat0, update_iter, ofs, etot_info, test_relax_method);
    }
    void reset_hessian_for_testing()
    {
        reset_hessian();
    }
    void save_bfgs_for_testing()
    {
        save_bfgs();
    }
    void update_inverse_hessian_for_testing(const double& lat0, std::ofstream& ofs)
    {
        update_inverse_hessian(lat0, ofs);
    }
    void check_wolfe_conditions_for_testing(std::ofstream& ofs, std::vector<double>& etot_info)
    {
        check_wolfe_conditions(ofs, etot_info);
    }
    void compute_trust_radius_for_testing(std::ofstream& ofs,
                                          std::vector<double>& etot_info,
                                          const int test_relax_method)
    {
        compute_trust_radius(ofs, etot_info, test_relax_method);
    }

  protected:
    void allocate_basic(void);
    void new_step(const double& lat0, int& update_iter, std::ofstream& ofs, std::vector<double>& etot_info, const int test_relax_method);
    void reset_hessian(void);
    void save_bfgs(void);

    std::vector<double> pos;  // std::vector containing 3N coordinates of the system ( x )
    std::vector<double> grad; // std::vector containing 3N components of ( grad( V(x) ) )
    std::vector<double> move; // pos = pos_p + move.

    std::vector<double> pos_p;  // p: previous
    std::vector<double> grad_p; // p: previous
    std::vector<double> move_p;

  public:
    static double relax_bfgs_w1; // fixed: parameters for Wolfe conditions.
    static double relax_bfgs_w2; // fixed: parameters for Wolfe conditions.

  protected:
    bool save_flag=false;
    bool tr_min_hit=false; //.TRUE. if the trust_radius has already been set
                     // to the minimum value at the previous step

    // mohan add 2010-07-27
    double check_move(const double& lat0, const double& pos, const double& pos_p);

  private:
    bool wolfe_flag=false;
    ModuleBase::matrix inv_hess;

    int bfgs_ndim;

    void update_inverse_hessian(const double& lat0, std::ofstream& ofs);
    void check_wolfe_conditions(std::ofstream& ofs, std::vector<double>& etot_info);
    void compute_trust_radius(std::ofstream& ofs, std::vector<double>& etot_info, const int test_relax_method);
};

#endif
