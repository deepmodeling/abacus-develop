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
    ESolver_NEP(const std::string& pot_file);

    void before_all_runners(UnitCell& ucell, const Input_para& inp) override;
    void runner(UnitCell& ucell, const int istep) override;
    double cal_energy() override;
    void cal_force(UnitCell& ucell, ModuleBase::matrix& force) override;
    void cal_stress(UnitCell& ucell, ModuleBase::matrix& stress) override;
    void after_all_runners(UnitCell& ucell) override;

  private:
    void type_map(const UnitCell& ucell);
#ifdef __NEP
    NEP3 nep;
#endif
    std::string nep_file;
    std::vector<int> atype = {};
    double nep_potential;
    ModuleBase::matrix nep_force;
    ModuleBase::matrix nep_virial;
    std::vector<double> _e;
    std::vector<double> _f;
    std::vector<double> _v;
};

} // namespace ModuleESolver

#endif