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
    void Init(Input& inp, UnitCell& cell) override;

    /**
     * @brief Run the DP solver for a given ion/md step and unit cell
     *
     * @param istep the current ion/md step
     * @param cell unitcell information
     */
    void Run(const int istep, UnitCell& cell) override;

    /**
     * @brief get the total energy without ion kinetic energy
     *
     * @param etot the computed energy
     * @return total energy without ion kinetic energy
     */
    double cal_Energy() override;

    /**
     * @brief get the computed atomic forces
     *
     * @param force the computed atomic forces
     */
    void cal_Force(ModuleBase::matrix& force) override;

    /**
     * @brief get the computed lattice virials
     *
     * @param stress the computed lattice virials
     */
    void cal_Stress(ModuleBase::matrix& stress) override;

    /**
     * @brief Prints the final total energy of the DP model to the output file
     *
     * This function prints the final total energy of the DP model in eV to the output file along with some formatting.
     */
    void postprocess() override;

  private:
    /**
     * @brief determine the type map of DP model
     *
     * @param ucell unitcell information
     * @return true if find keyword "type_map" in DP model
     * @return false if not find keyword "type_map" in DP model
     */
    bool type_map(const UnitCell& ucell);

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
     */
    std::string dp_file;          ///< the directory of DP model file
    std::vector<int> dp_type;     ///< convert atom type to dp type if find type_map
    std::vector<double> cell;     ///< the lattice vectors
    std::vector<int> atype;       ///< the atom type corresponding to DP model
    std::vector<double> coord;    ///< the atomic positions
    double dp_potential;          ///< the computed potential energy
    ModuleBase::matrix dp_force;  ///< the computed atomic forces
    ModuleBase::matrix dp_virial; ///< the computed lattice virials
    UnitCell* ucell_;             ///< pointer to the unitcell information
};

} // namespace ModuleESolver

#endif
