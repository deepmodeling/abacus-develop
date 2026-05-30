#ifndef MD_BASE_H
#define MD_BASE_H

#include "source_esolver/esolver.h"
#include "source_io/module_parameter/parameter.h"

#include <array>
#include <memory>
#include <string>
#include <vector>

#include "potential/potential_base.h"

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
    MD_base(const Parameter& param_in, UnitCell& unit_in);
    virtual ~MD_base();

    /**
     * @brief initialize an optional machine-learning potential.
     *
     * If this function is not called, MD keeps using the original ESolver
     * force/virial calculation path.
     *
     * @param potential_type backend name, e.g. "nep", "dpmd", or "deepmd"
     * @param model_file path to the model file, e.g. nep.txt or frozen_model.pb
     */
    void init_ml_potential(const std::string& potential_type,
                           const std::string& model_file);

    /**
     * @brief update potential energy, force, and virial.
     *
     * This is the single force-evaluation entrance used by both setup()
     * and the production MD loop.
     */
    void update_force_virial(ModuleESolver::ESolver* p_esolver);

    /**
     * @brief init before running md, calculate energy, force, and stress of the
     * initial configuration.
     * @param p_esolver the energy solver used in md
     * @param global_readin_dir directory of files for reading
     */
    virtual void setup(ModuleESolver::ESolver* p_esolver, const std::string& global_readin_dir);

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

  protected:
    /**
     * @brief restart MD when md_restart is true
     * @param global_readin_dir directory of files for reading
     */
    virtual void restart(const std::string& global_readin_dir);

    /**
     * @brief perform one step update of pos due to atomic velocity
     */
    virtual void update_pos();

    /**
     * @brief perform half-step update of vel due to atomic force
     * @param force atomic forces
     */
    virtual void update_vel(const ModuleBase::Vector3<double>* force);

  public:
    bool stop;                          ///< MD stop or not
    double t_current;                   ///< current temperature
    int step_;                          ///< the MD step finished in current calculation
    int step_rst_;                      ///< the MD step finished in previous calculations
    int frozen_freedom_;                ///< the fixed freedom of the system
    double* allmass = nullptr;                    ///< atom mass
    ModuleBase::Vector3<double>* pos;   ///< atom displacements  liuyu modify 2023-03-22
    ModuleBase::Vector3<double>* vel;   ///< atom velocity
    ModuleBase::Vector3<int>* ionmbl;   ///< atom is frozen or not
    ModuleBase::Vector3<double>* force; ///< force of each atom
    ModuleBase::matrix virial;          ///< virial for this lattice
    ModuleBase::matrix stress;          ///< stress for this lattice
    double potential=0.0;               ///< potential energy
    double kinetic;                     ///< kinetic energy

  protected:
    /**
     * @brief evaluate the optional ML potential.
     *
     * The ML potential adapter should return all quantities in ABACUS internal
     * MD units: energy in Hartree, force in Hartree/Bohr, and virial in Hartree.
     */
    void update_force_virial_ml();

    std::vector<std::array<double, 3>> pack_positions_for_ml() const;
    std::vector<int> pack_atom_types_for_ml() const;
    std::array<double, 9> pack_cell_for_ml() const;

    bool use_ml_potential_ = false;
    std::unique_ptr<ModuleMD::PotentialBase> ml_potential_;

  protected:
    const MD_para& mdp; ///< input parameters used in md
    UnitCell& ucell;    ///< unitcell information
    double energy_=0.0; ///< total energy of the system

    bool cal_stress;  ///< whether calculate stress
    int my_rank;      ///< MPI rank of the processor
    double md_dt;     ///< Time increment (hbar/E_hartree)
    double md_tfirst; ///< Temperature (in Hartree, 1 Hartree ~ 3E5 K)
    double md_tlast;  ///< Target temperature
};

#endif // MD_BASE_H