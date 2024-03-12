#ifndef MSST_H
#define MSST_H

#include "md_base.h"

/**
 * @brief MSST method
 *
 * Multi-Scale Shock Technique (MSST) integration [Phys. Rev. Lett. 90, 235503 (2003)]
 * update positions and velocities each timestep to mimic a compressive shock wave passing over the system.
 * The MSST varies the cell volume and temperature to restrain the system to the shock Hugoniot and the Rayleigh line.
 * These restraints correspond to the macroscopic conservation laws dictated by a shock front.
 */
class MSST : public MD_base
{
  public:
    MSST(MD_para& MD_para_in, UnitCell& unit_in);
    ~MSST();

  private:
    void setup(ModuleESolver::ESolver* p_esolver, const std::string& global_readin_dir);
    void first_half(std::ofstream& ofs);
    void second_half();
    void print_md(std::ofstream& ofs, const bool& cal_stress);
    void write_restart(const std::string& global_out_dir);
    void restart(const std::string& global_readin_dir);

    /**
     * @brief get the sum of square of velocities
     *
     * @return the sum of square of velocities
     */
    double vel_sum();

    /**
     * @brief rescale the lattice and velocities
     *
     * @param ofs determine the output files
     * @param volume the current cell volume
     */
    void rescale(std::ofstream& ofs, const double& volume);

    /**
     * @brief propagate atomic velocities
     *
     */
    void propagate_vel();

    /**
     * @brief propagate the volume change rate
     *
     */
    void propagate_voldot();

    ModuleBase::Vector3<double>* old_v;   ///< old atomic velocities
    ModuleBase::Vector3<double> dilation; ///< dilation scale
    ModuleBase::Vector3<double> omega;    ///< time derivative of volume
    double p0;                            ///< initial pressure
    double v0;                            ///< initial volume
    double e0;                            ///< initial energy
    double totmass;                       ///< total mass of the cell
    double lag_pos;                       ///< Lagrangian location of cell
    double vsum;                          ///< sum over v^2
};

#endif