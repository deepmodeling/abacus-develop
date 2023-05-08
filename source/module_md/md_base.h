#ifndef MDRUN_H
#define MDRUN_H

#include "md_para.h"
#include "module_esolver/esolver.h"

class MD_base
{
  public:
    MD_base(MD_parameters& MD_para_in, UnitCell& unit_in);
    virtual ~MD_base();

    /**
     * @brief init before running md, calculate energy, force, and stress of the initial configuration.
     * @param p_esolver the energy solver used in md
     * @param my_rank MPI rank of the processer
     * @param global_readin_dir directory of files for reading
     */
    virtual void setup(ModuleESolver::ESolver* p_esolver, const int& my_rank, const std::string& global_readin_dir);

    /**
     * @brief the first half of equation of motion, update velocities and positions
     * @param my_rank MPI rank of the processer
     * @param ofs determine the output files
     */
    virtual void first_half(const int& my_rank, std::ofstream& ofs);

    /**
     * @brief the second half of equation of motion, update velocities
     * @param my_rank MPI rank of the processer
     */
    virtual void second_half(const int& my_rank);

    /**
     * @brief output MD information such as energy, temperature, and pressure
     * @param ofs determine the output files
     * @param cal_stress whether calculate and output stress
     * @param my_rank MPI rank of the processer
     */
    virtual void outputMD(std::ofstream& ofs, const bool& cal_stress, const int& my_rank);

    /**
     * @brief write the information into files used for MD restarting
     * @param my_rank MPI rank of the processer
     * @param global_out_dir directory of output files
     */
    virtual void write_restart(const int& my_rank, const std::string& global_out_dir);

  protected:
    /**
     * @brief restart MD when md_restart is true
     * @param my_rank MPI rank of the processer
     * @param global_readin_dir directory of files for reading
     */
    virtual void restart(const int& my_rank, const std::string& global_readin_dir);

    /**
     * @brief perform one step update of pos due to atomic velocity
     * @param my_rank MPI rank of the processer
     */
    virtual void update_pos(const int& my_rank);

    /**
     * @brief perform half-step update of vel due to atomic force
     * @param force atomic forces
     * @param my_rank MPI rank of the processer
     */
    virtual void update_vel(const ModuleBase::Vector3<double>* force, const int& my_rank);

    // All parameters are in a.u. unit.
  public:
    bool stop;                          // MD stop or not
    double t_current;                   // current temperature
    int step_;                          // the MD step finished in current calculation
    int step_rst_;                      // the MD step finished in previous calculations
    int frozen_freedom_;                // the fixed freedom of the system
    double* allmass;                    // atom mass
    ModuleBase::Vector3<double>* pos;   // atom displacements  liuyu modify 2023-03-22
    ModuleBase::Vector3<double>* vel;   // atom velocity
    ModuleBase::Vector3<int>* ionmbl;   // atom is frozen or not
    ModuleBase::Vector3<double>* force; // force of each atom
    ModuleBase::matrix virial;          // virial for this lattice
    ModuleBase::matrix stress;          // stress for this lattice
    double potential;                   // potential energy
    double kinetic;                     // kinetic energy

  protected:
    MD_parameters& mdp;
    UnitCell& ucell;
    double energy_; // total energy of the system
};

#endif