#ifndef NOSE_HOOVER_H
#define NOSE_HOOVER_H

#include "md_base.h"

class Nose_Hoover : public MD_base
{
  public:
    Nose_Hoover(MD_parameters& MD_para_in, UnitCell& unit_in);
    ~Nose_Hoover();

  private:
    void setup(ModuleESolver::ESolver* p_esolver, const int& my_rank, const std::string& global_readin_dir);
    void first_half(const int& my_rank, std::ofstream& ofs);
    void second_half(const int& my_rank);
    void outputMD(std::ofstream& ofs, const bool& cal_stress, const int& my_rank);
    void write_restart(const int& my_rank, const std::string& global_out_dir);
    void restart(const int& my_rank, const std::string& global_readin_dir);

    // perform half-step update of thermostats coupled with particles
    void particle_thermo();

    // perform half-step update of thermostats coupled with barostat
    void baro_thermo();

    // perform half-step update of barostat velocity
    void update_baro();

    // perform half-step update of velocity due to barostat
    void vel_baro();

    // determine target stress
    void target_stress();

    // couple stress component due to md_pcouple
    void couple_stress();

    // perform half-step update of volume
    void update_volume(std::ofstream& ofs);

    const int nc_tchain = 1;
    const int nc_pchain = 1;
    const static int nys = 7;
    double w[nys]; // scale evolution operator

    // thermostats
    int tdof;         // particle degree of freedom
    double t_target;  // target temperature
    double* mass_eta; // mass of thermostats coupled with particles
    double* eta;      // position of thermostats coupled with particles
    double* v_eta;    // velocity of thermostats coupled with particles
    double* g_eta;    // acceleration of thermostats coupled with particles

    // barostat, Voigt notation: x, y, z, yz, xz, xy
    int npt_flag;         // whether NPT ensemble
    double mass_omega[6]; // mass of lattice component
    double v_omega[6];    // velocity of lattice component
    double pstart[6];     // initial stress components
    double pstop[6];      // final stress components
    double pfreq[6];      // Oscillation frequency, used to determine qmass of thermostats coupled with barostat
    int pflag[6];         // control stress components
    int pdim;             // pdim = pflag[0] + pflag[1] + pflag[2], number of barostatted dims
    double p_target[6];   // target stress components
    double p_hydro;       // target hydrostatic target pressure
    double p_current[6];  // current stress after coupled
    double* mass_peta;    // mass of thermostats coupled with barostat
    double* peta;         // position of thermostats coupled with barostat
    double* v_peta;       // velocity of thermostats coupled with barostat
    double* g_peta;       // acceleration of thermostats coupled with barostat
    double mtk_term;      // mtk correction
};

#endif