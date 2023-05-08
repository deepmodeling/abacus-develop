#ifndef VERLET_H
#define VERLET_H

#include "md_base.h"

class Verlet : public MD_base
{
  public:
    Verlet(MD_parameters& MD_para_in, UnitCell& unit_in);
    ~Verlet();

  private:
    void setup(ModuleESolver::ESolver* p_esolver, const int& my_rank, const std::string& global_readin_dir);
    void first_half(const int& my_rank, std::ofstream& ofs);
    void second_half(const int& my_rank);
    void restart(const int& my_rank, const std::string& global_readin_dir);
    void outputMD(std::ofstream& ofs, const bool& cal_stress, const int& my_rank);
    void apply_thermostat(const int& my_rank);
    void thermalize(const int& nraise, const double& current_temp, const double& target_temp);
    void write_restart(const int& my_rank, const std::string& global_out_dir);
};

#endif