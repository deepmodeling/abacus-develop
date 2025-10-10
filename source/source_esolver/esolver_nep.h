#ifndef ESOLVER_NEP_H
#define ESOLVER_NEP_H

#include "esolver.h"
#ifdef __NEP
#include "nep.h"
#endif
#include <vector>
#include <string>

namespace ModuleESolver
{

class ESolver_NEP : public ESolver
{
  public:
#ifdef __NEP
    ESolver_NEP(const std::string& pot_file): nep(pot_file)
  {
      classname = "ESolver_NEP";
      nep_file = pot_file;
  }
#else
    ESolver_NEP(const std::string& pot_file)
  {
      classname = "ESolver_NEP";
      nep_file = pot_file;
  }
#endif

    /**
     * @brief Initialize the NEP solver with given input parameters and unit cell
     *
     * @param inp input parameters
     * @param cell unitcell information
     */
    void before_all_runners(UnitCell& ucell, const Input_para& inp) override;
  
    /**
     * @brief Run the NEP solver for a given ion/md step and unit cell
     *
     * @param istep the current ion/md step
     * @param cell unitcell information
     */
    void runner(UnitCell& ucell, const int istep) override;

    /**
     * @brief get the total energy without ion kinetic energy
     *
     * @param etot the computed energy
     * @return total energy without ion kinetic energy
     */
    double cal_energy() override;

    /**
     * @brief get the computed atomic forces
     *
     * @param force the computed atomic forces
     */
    void cal_force(UnitCell& ucell, ModuleBase::matrix& force) override;

    /**
     * @brief get the computed lattice virials
     *
     * @param stress the computed lattice virials
     */
    void cal_stress(UnitCell& ucell, ModuleBase::matrix& stress) override;

    /**
     * @brief Prints the final total energy of the NEP model to the output file
     *
     * This function prints the final total energy of the NEP model in eV to the output file along with some formatting.
     */
    void after_all_runners(UnitCell& ucell) override;

  private:
    /**
     * @brief determine the type map of NEP model
     *
     * @param ucell unitcell information
     */
    void type_map(const UnitCell& ucell);

    /**
     * @brief NEP related variables for ESolver_NEP class
     *
     * These variables are related to the NEP method and are used in the ESolver_NEP class to compute the potential
     * energy and forces.
     *
     * @note These variables are only defined if the __NEP preprocessor macro is defined.
     */
#ifdef __NEP
    NEP3 nep; ///< NEP3 object for NEP calculations
#endif

    std::string nep_file;                ///< directory of NEP model file
    std::vector<int> atype = {};         ///< atom type mapping for NEP model
    double nep_potential;                ///< computed potential energy
    ModuleBase::matrix nep_force;        ///< computed atomic forces
    ModuleBase::matrix nep_virial;       ///< computed lattice virials
    std::vector<double> _e;              ///< temporary storage for energy computation
    std::vector<double> _f;              ///< temporary storage for force computation
    std::vector<double> _v;              ///< temporary storage for virial computation
};

} // namespace ModuleESolver

#endif