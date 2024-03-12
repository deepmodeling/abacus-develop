#ifndef FIRE_H
#define FIRE_H

#include "md_base.h"

/**
 * @brief FIRE method
 *
 * A MD-based relaxation algorithm, named fast inertial relaxation engine. [Phys. Rev. Lett. 97, 170201 (2006)]
 * It is based on conventional molecular dynamics with additional velocity modifications and adaptive time steps.
 * The MD trajectory will descend to an energy-minimum.
 */
class FIRE : public MD_base
{
  public:
    FIRE(MD_para& MD_para_in, UnitCell& unit_in);
    ~FIRE();

  private:
    void setup(ModuleESolver::ESolver* p_esolver, const std::string& global_readin_dir);
    void first_half(std::ofstream& ofs);
    void second_half();
    void print_md(std::ofstream& ofs, const bool& cal_stress);
    void restart(const std::string& global_readin_dir);
    void write_restart(const std::string& global_out_dir);

    /**
     * @brief check the atomic forces converged or not
     */
    void check_force();

    /**
     * @brief update related parameters
     */
    void check_fire();

    double max;         ///< max force
    double alpha_start; ///< alpha_start begin
    double alpha;       ///< alpha begin
    double finc;        ///< finc begin
    double fdec;        ///< fdec begin
    double f_alpha;     ///< f_alpha
    int n_min;          ///< n_min
    double dt_max;      ///< dt_max
    int negative_count; ///< Negative count
};

#endif