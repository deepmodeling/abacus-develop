#ifndef ESOLVER_NEP_H
#define ESOLVER_NEP_H

#include "esolver.h"
#include "esolver_nep_postprocess.h"
#ifdef __NEP
#include "nep.h"
#endif
#ifdef __CUDA
#include "neighbor_nep.h"
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

    void before_all_runners(UnitCell& ucell, const Input_para& inp) override;
    void runner(UnitCell& ucell, const int istep) override;
    double cal_energy() override;
    void cal_force(UnitCell& ucell, ModuleBase::matrix& force) override;
    void cal_stress(UnitCell& ucell, ModuleBase::matrix& stress) override;
    void after_all_runners(UnitCell& ucell) override;

  private:
    void prepare_input_buffers(const UnitCell& ucell);
    void postprocess_outputs(const UnitCell& ucell);
    void type_map(const UnitCell& ucell);

#ifdef __NEP
    NEP nep;
#endif

    std::string nep_file;
    std::vector<int> atype = {};
    double nep_potential;
    ModuleBase::matrix nep_force;
    ModuleBase::matrix nep_virial;
    std::vector<double> _e;
    std::vector<double> _f;
    std::vector<double> _v;
    std::vector<double> cell;
    std::vector<double> coord;
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
