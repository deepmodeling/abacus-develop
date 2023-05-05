#ifndef FIRE_H
#define FIRE_H

#include "md_base.h"

class FIRE : public MDrun
{
  public:
    FIRE(MD_parameters &MD_para_in, UnitCell &unit_in);
    ~FIRE();

    void setup(ModuleESolver::ESolver *p_esolver,
               const int &my_rank,
               const std::string &global_readin_dir,
               const double &force_thr);
    void first_half(const int &my_rank);
    void second_half(const int &my_rank, const double &force_thr);
    void outputMD(std::ofstream &ofs, const bool &cal_stress, const int &my_rank);
    void write_restart(const int &my_rank, const std::string &global_out_dir);
    void restart(const int &my_rank, const std::string &global_readin_dir);
    void check_force(const double &force_thr);
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