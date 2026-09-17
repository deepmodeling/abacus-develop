#ifndef MD_BASE_H
#define MD_BASE_H

#include "source_cell/mdcell.h"
#include "source_esolver/esolver.h"
#include "source_io/module_parameter/md_parameter.h"

#include <cstdint>

class DomainDecomposition;

/**
 * @brief base class of md
 *
 * This class implements the velocity-Verlet method.
 * The system is assumed to be isolated in the sense that it cannot exchange
 * energy or particles with its environment, so that the energy of the system
 * does not change with time.
 */
class MD_base
{
  public:
    /**
     * @brief construct the integrator from the values it actually uses
     * @param mdp_in the md input parameters; the reference is kept, so it must
     *               outlive the integrator
     * @param cal_stress_in whether stress is calculated
     * @param init_vel whether initial velocities are read from STRU
     * @param my_rank_in MPI rank of the processor; only consulted in serial
     *                   builds, where MDCell cannot supply it
     * @param mdcell_in mdcell information
     */
    MD_base(const MD_para& mdp_in,
            const bool cal_stress_in,
            const bool init_vel,
            const int my_rank_in,
            MDCell& mdcell_in);
    virtual ~MD_base();

    /**
     * @brief init before running md, calculate energy, force, and stress of the
     * initial configuration.
     * @param p_esolver the energy solver used in md
     * @param global_readin_dir directory of files for reading
     */
    virtual void setup(ModuleESolver::ESolver* p_esolver, const std::string& global_readin_dir, DomainDecomposition& decomp);

    /**
     * @brief the first half of equation of motion, update velocities and
     * positions
     * @param ofs determine the output files
     */
    virtual void first_half(std::ofstream& ofs);

    /**
     * @brief the second half of equation of motion, update velocities
     */
    virtual void second_half();

    /**
     * @brief output MD information such as energy, temperature, and pressure
     * @param ofs determine the output files
     * @param cal_stress whether calculate and output stress
     */
    virtual void print_md(std::ofstream& ofs, const bool& cal_stress);

    /**
     * @brief write the information into files used for MD restarting
     * @param global_out_dir directory of output files
     */
    virtual void write_restart(const std::string& global_out_dir);

    /**
     * @brief restart MD when md_restart is true
     *
     * Public counterpart of write_restart(): setup() calls it internally when
     * md_restart is set, and a caller that has just written a restart file may
     * read it back through here.
     *
     * @param global_readin_dir directory of files for reading
     */
    virtual void restart(const std::string& global_readin_dir);

  protected:
    /**
     * @brief perform one step update of pos due to atomic velocity
     */
    virtual void update_pos();

    /**
     * @brief perform half-step update of vel due to atomic force
     * @param force atomic forces
     */
    virtual void update_vel();

  public:
    bool stop;                          ///< MD stop or not
    double t_current;                   ///< current temperature
    int step_;                          ///< the MD step finished in current calculation
    int step_rst_;                      ///< the MD step finished in previous calculations
    std::int64_t frozen_freedom_;       ///< the fixed freedom of the system
    ModuleBase::matrix virial;          ///< virial for this lattice
    ModuleBase::matrix stress;          ///< stress for this lattice
    double potential=0.0;               ///< potential energy
    double kinetic;                     ///< kinetic energy

  protected:
    const MD_para& mdp; ///< input parameters used in md
    MDCell& mdcell;     ///< mdcell information
    double energy_=0.0; ///< total energy of the system

    bool cal_stress;  ///< whether calculate stress
    int my_rank;      ///< MPI rank of the processor
    double md_dt;     ///< Time increment (hbar/E_hartree)
    double md_tfirst; ///< Temperature (in Hartree, 1 Hartree ~ 3E5 K)
    double md_tlast;  ///< Target temperature
};

#endif // MD_BASE_H
