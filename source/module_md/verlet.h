#ifndef VERLET_H
#define VERLET_H

#include "md_base.h"

class Verlet : public MDrun
{
  public:
    Verlet(MD_parameters &MD_para_in, UnitCell &unit_in);
    ~Verlet();

    void setup(ModuleESolver::ESolver *p_esolver, const int &my_rank, const std::string &global_readin_dir);
    void first_half(const int &my_rank, std::ofstream &ofs);
    void second_half(const int &my_rank);
    void apply_thermostat(const int &my_rank);
    void thermalize(const int &nraise, const double &current_temp, const double &target_temp);
    void outputMD(std::ofstream &ofs, const bool &cal_stress, const int &my_rank);
    void write_restart(const int &my_rank, const std::string &global_out_dir);
    void restart(const int &my_rank, const std::string &global_readin_dir);
};

#endif