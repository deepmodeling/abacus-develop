#ifndef ESOLVER_NEP_H
#define ESOLVER_NEP_H

#include "esolver.h"
#include "esolver_nep_postprocess.h"
#ifdef __NEP
#include "nep.h"
#endif
#include <string>
#include <vector>

namespace ModuleESolver
{

class ESolver_NEP : public ESolver
{
  public:
#ifdef __NEP
    ESolver_NEP(const std::string& pot_file) : nep(pot_file)
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
    void before_all_runners(BaseCell& basecell, const Input_para& inp) override;

    /**
     * @brief Run the NEP solver for a given ion/md step and unit cell
     *
     * @param istep the current ion/md step
     * @param cell unitcell information
     */
    void runner(BaseCell& basecell, const int istep) override;

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
    void cal_force(BaseCell& basecell, ModuleBase::matrix& force) override;

    /**
     * @brief get the computed lattice virials
     *
     * @param stress the computed lattice virials
     */
    void cal_stress(BaseCell& basecell, ModuleBase::matrix& stress) override;

    /**
     * @brief Prints the final total energy of the NEP model to the output file
     *
     * This function prints the final total energy of the NEP model in eV to the output file along with some formatting.
     */
    void after_all_runners(BaseCell& basecell) override;

  private:
    void prepare_input_buffers(const UnitCell& ucell);
    void postprocess_outputs(const UnitCell& ucell);
    void type_map(const UnitCell& ucell);

#ifdef __NEP
    NEP nep;
#endif

    std::string nep_file;          ///< directory of NEP model file
    std::vector<int> atype = {};   ///< atom type mapping for NEP model
    double nep_potential;          ///< computed potential energy
    ModuleBase::matrix nep_force;  ///< computed atomic forces
    ModuleBase::matrix nep_virial; ///< computed lattice virials
    std::vector<double> _e;        ///< temporary storage for energy computation
    std::vector<double> _f;        ///< temporary storage for force computation
    std::vector<double> _v;        ///< temporary storage for virial computation
    std::vector<double> cell;
    std::vector<double> coord;
    bool use_gpu_ = false;
#ifdef __CUDA
    NepCudaPostprocessWorkspace cuda_postprocess_workspace;

    // Neighbor list buffers for GPU compute path
    static constexpr int NEP_GPU_MN = 1000;
    int num_cells[3];
    double ebox[18];
    std::vector<int> g_NN_radial;
    std::vector<int> g_NL_radial;
    std::vector<int> g_NN_angular;
    std::vector<int> g_NL_angular;
    std::vector<double> r12;
#endif
};

} // namespace ModuleESolver

#endif
