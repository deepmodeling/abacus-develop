#ifndef VERLET_H
#define VERLET_H

#include "md_base.h"

/**
 * @brief the md methods based on the velocity-Verlet equation
 *
 */
class Verlet : public MD_base
{
  public:
    Verlet(MD_para& MD_para_in, UnitCell& unit_in);
    ~Verlet();

  private:
    void setup(ModuleESolver::ESolver* p_esolver, const std::string& global_readin_dir);
    void first_half(std::ofstream& ofs);
    void second_half();
    void restart(const std::string& global_readin_dir);
    void print_md(std::ofstream& ofs, const bool& cal_stress);
    void write_restart(const std::string& global_out_dir);

    /**
     * @brief apply specifical thermostats according to the input para
     *
     */
    void apply_thermostat();

    /**
     * @brief rescale atomic velocities
     *
     * @param nraise a parameter related to thermostats
     * @param current_temp the current temperature
     * @param target_temp the target temperature
     */
    void thermalize(const int& nraise, const double& current_temp, const double& target_temp);
};

#endif