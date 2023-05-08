#ifndef MSST_H
#define MSST_H

#include "md_base.h"

class MSST : public MD_base
{
  public:
    MSST(MD_parameters& MD_para_in, UnitCell& unit_in);
    ~MSST();

  private:
    void setup(ModuleESolver::ESolver* p_esolver, const int& my_rank, const std::string& global_readin_dir);
    void first_half(const int& my_rank, std::ofstream& ofs);
    void second_half(const int& my_rank);
    void outputMD(std::ofstream& ofs, const bool& cal_stress, const int& my_rank);
    void write_restart(const int& my_rank, const std::string& global_out_dir);
    void restart(const int& my_rank, const std::string& global_readin_dir);
    double extra_term();
    double vel_sum();
    void rescale(std::ofstream& ofs, const double& volume);
    void propagate_vel(const int& my_rank);
    void propagate_voldot();

    ModuleBase::Vector3<double>* old_v;
    ModuleBase::Vector3<double> dilation; // dilation scale
    ModuleBase::Vector3<double> omega;    // time derivative of volume
    double p0;                            // initial pressure
    double v0;                            // initial volume
    double e0;                            // initial energy
    double totmass;                       // total mass of the cell
    double lag_pos;                       // Lagrangian location of cell
    double vsum;                          // sum over v^2
};

#endif