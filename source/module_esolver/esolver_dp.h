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

/**
 * @brief esolver dp
 *
 * Take DP machine learning model as energy esolver.
 */
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
     * @brief initialize related variables
     *
     * @param inp input parameters
     * @param cell unitcell information
     */
    void Init(Input& inp, UnitCell& cell) override;

    /**
     * @brief run the energy esolver
     *
     * Using DP model, energy, forces, and virials are obtained
     *
     * @param istep the current ion/md step
     * @param cell unitcell information
     */
    void Run(const int istep, UnitCell& cell) override;

    /**
     * @brief get the total energy without ion kinetic energy
     *
     * @param etot the computed energy
     */
    void cal_Energy(double& etot) override;

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
     * @brief the postprocess of esolver dp
     *
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

    /// the DP model
#ifdef __DPMD
#ifdef __DPMDC
    deepmd::hpp::DeepPot dp;
#else
    deepmd::DeepPot dp;
#endif
#endif

    std::string dp_file;      ///< the directory of DP model file
    std::vector<int> dp_type; ///< convert atom type to dp type if find type_map
    std::vector<double> cell; ///< the lattice vectors
    std::vector<int> atype;   ///< the atom type corresponding to DP model
    std::vector<double> coord;    ///< the atomic positions
    double dp_potential;          ///< the computed potential energy
    ModuleBase::matrix dp_force;  ///< the computed atomic forces
    ModuleBase::matrix dp_virial; ///< the computed lattice virials
};

} // namespace ModuleESolver

#endif
