#ifndef MD_FUNC_H
#define MD_FUNC_H

#include "module_esolver/esolver.h"

namespace MD_func
{
double gaussrand();

void InitVel(const UnitCell& unit_in,
             const double& temperature,
             double* allmass,
             int& frozen_freedom,
             ModuleBase::Vector3<int>* ionmbl,
             ModuleBase::Vector3<double>* vel);

void ReadVel(const UnitCell& unit_in, ModuleBase::Vector3<double>* vel);

void RandomVel(const int& numIon,
               const double& temperature,
               const double* allmass,
               const int& frozen_freedom,
               const ModuleBase::Vector3<int> frozen,
               const ModuleBase::Vector3<int>* ionmbl,
               ModuleBase::Vector3<double>* vel);

void force_virial(ModuleESolver::ESolver* p_esolver,
                  const int& istep,
                  UnitCell& unit_in,
                  double& potential,
                  ModuleBase::Vector3<double>* force,
                  ModuleBase::matrix& virial);

double GetAtomKE(const int& numIon, const ModuleBase::Vector3<double>* vel, const double* allmass);

void compute_stress(const UnitCell& unit_in,
                    const ModuleBase::Vector3<double>* vel,
                    const double* allmass,
                    const ModuleBase::matrix& virial,
                    ModuleBase::matrix& stress);

void outStress(const ModuleBase::matrix& virial, const ModuleBase::matrix& stress);

void MDdump(const int& step,
            const UnitCell& unit_in,
            const MD_parameters& mdp,
            const ModuleBase::matrix& virial,
            const ModuleBase::Vector3<double>* force,
            const ModuleBase::Vector3<double>* vel);

void getMassMbl(const UnitCell& unit_in,
                double* allmass,
                ModuleBase::Vector3<int>& frozen,
                ModuleBase::Vector3<int>* ionmbl);

double target_temp(const int& istep, const int& nstep, const double& tfirst, const double& tlast);

double current_temp(double& kinetic,
                    const int& natom,
                    const int& frozen_freedom,
                    const double* allmass,
                    const ModuleBase::Vector3<double>* vel);

void temp_vector(const int& natom,
                 const ModuleBase::Vector3<double>* vel,
                 const double* allmass,
                 ModuleBase::matrix& t_vector);

double current_step(const int& my_rank, const std::string& file_dir);

} // namespace MD_func

#endif
