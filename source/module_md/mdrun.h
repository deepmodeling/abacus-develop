#ifndef MDRUN_H
#define MDRUN_H

#include "MD_parameters.h"
#include "../module_cell/unitcell.h"
#include "../module_base/matrix.h"
#include "module_esolver/esolver.h"

class MDrun
{
public:
    MDrun(MD_parameters& MD_para_in, UnitCell &unit_in);
    virtual ~MDrun();

    /**
     * @brief init before running md, calculate energy, force, and stress of the initial configuration.
     * @param p_esolver the energy solver used in md
     */
    virtual void setup(ModuleESolver::ESolver *p_esolver);

    /**
     * @brief the first half of equation of motion, update velocities and positions
     */
    virtual void first_half();

    /**
     * @brief the second half of equation of motion, update velocities
     */
    virtual void second_half();

    /**
     * @brief output MD information such as energy, temperature, and pressure
     * @param ofs determine the output files
     * @param cal_stress whether calculate and output stress
     */
    virtual void outputMD(std::ofstream &ofs, bool cal_stress);

    /**
     * @brief write the information into files used for MD restarting
     */
    virtual void write_restart();

    /**
     * @brief restart MD when md_restart is true
     */
    virtual void restart();

    MD_parameters &mdp;
	UnitCell &ucell;
    bool stop;                           // MD stop or not

    // All parameters are in a.u. unit.
    double t_current;                    // current temperature
    int step_;                           // the MD step finished in current calculation
    int step_rst_;                       // the MD step finished in previous calculations
    double energy_;                      // total energy of the system
	int frozen_freedom_;                 // the fixed freedom of the system

    double *allmass;                     // atom mass 
    ModuleBase::Vector3<double> *pos;    // atom position
    ModuleBase::Vector3<double> *vel;    // atom velocity
    ModuleBase::Vector3<int> *ionmbl;    // atom is frozen or not
    ModuleBase::Vector3<double> *force;  // force of each atom
    ModuleBase::matrix virial;           // virial for this lattice
	ModuleBase::matrix stress;           // stress for this lattice
    double potential;                    // potential energy
    double kinetic;                      // kinetic energy

};

#endif