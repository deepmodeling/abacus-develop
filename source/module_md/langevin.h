#ifndef LANGEVIN_H
#define LANGEVIN_H

#include "md_base.h"

/**
 * @brief Langevin method
 *
 * Assume the atoms are embedded in a sea of much smaller fictional particles.
 * The solvent influences the dynamics of the solute(typically nanoparticles) via random collisions,
 * and by imposing a frictional drag force on the motion of the nanoparticle in the solvent.
 * The damping factor and the random force combine to give the correct NVT ensemble.
 */
class Langevin : public MD_base
{
  public:
    Langevin(MD_para& MD_para_in, UnitCell& unit_in);
    ~Langevin();

  private:
    void setup(ModuleESolver::ESolver* p_esolver, const std::string& global_readin_dir);
    void first_half(std::ofstream& ofs);
    void second_half();
    void print_md(std::ofstream& ofs, const bool& cal_stress);
    void write_restart(const std::string& global_out_dir);
    void restart(const std::string& global_readin_dir);

    /**
     * @brief calculate fictitious forces
     *
     */
    void post_force();

    ModuleBase::Vector3<double>* total_force; ///< total force = true force + Langevin fictitious_force
};

#endif