#ifndef FIRE_H
#define FIRE_H

#include "md_base.h"

class FIRE : public MD_base
{
  public:
    FIRE(MD_parameters& MD_para_in, UnitCell& unit_in);
    ~FIRE();

  private:
    void setup(ModuleESolver::ESolver* p_esolver, const int& my_rank, const std::string& global_readin_dir);
    void first_half(const int& my_rank, std::ofstream& ofs);
    void second_half(const int& my_rank);
    void outputMD(std::ofstream& ofs, const bool& cal_stress, const int& my_rank);
    void restart(const int& my_rank, const std::string& global_readin_dir);
    void write_restart(const int& my_rank, const std::string& global_out_dir);
    void check_force();
    void check_FIRE();

    double max;         // max force
    double alpha_start; // alpha_start begin
    double alpha;       // alpha begin
    double finc;        // finc begin
    double fdec;        // fdec begin
    double f_alpha;
    int N_min;
    double dt_max;      // dt_max
    int negative_count; // Negative count
};

#endif