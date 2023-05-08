#ifndef LANGEVIN_H
#define LANGEVIN_H

#include "md_base.h"

class Langevin : public MD_base
{
  public:
    Langevin(MD_parameters& MD_para_in, UnitCell& unit_in);
    ~Langevin();

  private:
    void setup(ModuleESolver::ESolver* p_esolver, const int& my_rank, const std::string& global_readin_dir);
    void first_half(const int& my_rank, std::ofstream& ofs);
    void second_half(const int& my_rank);
    void outputMD(std::ofstream& ofs, const bool& cal_stress, const int& my_rank);
    void write_restart(const int& my_rank, const std::string& global_out_dir);
    void restart(const int& my_rank, const std::string& global_readin_dir);
    void post_force();

    ModuleBase::Vector3<double>* total_force; // total force = true force + Langevin fictitious_force
};

#endif