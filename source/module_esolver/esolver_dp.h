#ifndef ESOLVER_DP_H
#define ESOLVER_DP_H

#include "esolver.h"
#ifdef __DPMD
#ifdef __DPMDC
#include "deepmd/deepmd.hpp"
#else
#include "deepmd/DeepPot.h"
#endif
#endif

namespace ModuleESolver
{

class ESolver_DP : public ESolver
{
  public:
#ifdef __DPMD
    ESolver_DP(const std::string& pot_file) : dp(pot_file)
    {
        classname = "ESolver_DP";
        dp_file = pot_file;
    }
#else
    ESolver_DP(const std::string& pot_file)
    {
        classname = "ESolver_DP";
        dp_file = pot_file;
    }
#endif

    /**
     * @brief Initialize the DP solver with given input parameters and unit cell
     *
     * @param inp input parameters
     * @param cell unitcell information
     */
    void before_all_runners(const Input_para& inp, UnitCell& cell) override;

    /**
     * @brief Run the DP solver for a given ion/md step and unit cell
     *
     * @param istep the current ion/md step
     * @param cell unitcell information
     */
    void runner(const int istep, UnitCell& cell) override;

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
    void cal_force(ModuleBase::matrix& force) override;

    /**
     * @brief get the computed lattice virials
     *
     * @param stress the computed lattice virials
     */
    void cal_stress(ModuleBase::matrix& stress) override;

    /**
     * @brief Prints the final total energy of the DP model to the output file
     *
     * This function prints the final total energy of the DP model in eV to the output file along with some formatting.
     */
    void after_all_runners() override;

  private:
    /**
     * @brief determine the type map of DP model
     *
     * @param ucell unitcell information
     */
    void type_map(const UnitCell& ucell);

    /**
     * @brief DeePMD related variables for ESolver_DP class
     *
     * These variables are related to the DeePMD method and are used in the ESolver_DP class to compute the potential
     * energy and forces.
     *
     * @note These variables are only defined if the __DPMD preprocessor macro is defined.
     */
#ifdef __DPMD
#ifdef __DPMDC
    deepmd::hpp::DeepPot dp; ///< C interface
#else
    deepmd::DeepPot dp; ///< C++ interface
#endif
#endif

    /**
     * @brief Variables for storing simulation data in ESolver_DP class
     *
     * These variables are used in the ESolver_DP class to store simulation data such as atomic positions, types, and
     * the potential energy and forces.
     *
     * @param dp_file the directory of DP model file
     * @param cell the lattice vectors
     * @param atype the atom type corresponding to DP model
     * @param coord the atomic positions
     * @param fparam The frame parameter for dp potential. The array can be of size:
     *               dim_fparam. Then all frames are assumed to be provided with the same fparam.
     * @param aparam The atomic parameterfor dp potential. The array can be of size:
     *               natoms x dim_aparam. Then all frames are assumed to be provided with the same aparam.
     *               dim_aparam. Then all frames and atoms are assumed to be provided with the same aparam.
     * @param dp_potential the computed potential energy
     * @param dp_force the computed atomic forces
     * @param dp_virial the computed lattice virials
     * @param ucell_ pointer to the unitcell information
     */
    std::string dp_file;
    std::vector<double> cell = {};
    std::vector<int> atype = {};
    std::vector<double> coord = {};
    std::vector<double> fparam = {};
    std::vector<double> aparam = {};
    double dp_potential;
    ModuleBase::matrix dp_force;
    ModuleBase::matrix dp_virial;
    UnitCell* ucell_;
};

} // namespace ModuleESolver

#endif
