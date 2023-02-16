#ifndef H_TDDFT_PW_H
#define H_TDDFT_PW_H

#include "pot_base.h"
#include "module_io/input.h"

namespace elecstate
{

class H_TDDFT_pw : public PotBase
{
  public:
    H_TDDFT_pw(const ModulePW::PW_Basis* rho_basis_in, const UnitCell* ucell_in) : ucell_(ucell_in)
    {
        this->dynamic_mode = false;
        this->fixed_mode = true;

        this->rho_basis_ = rho_basis_in;
    }
    ~H_TDDFT_pw(){};

    void cal_fixed_v(double* vl_pseudo) override;

  private:
    // internal time-step,
    //-------hypothesis-------
    // Vext will evolve by time, every time cal_fixed_v() is called, istep++
    //------------------------
    static int istep;

    // parameters
    int stype ; // 0 : length gauge  1: velocity gauge

    std::vector<int> ttype ;
    //  0  Gauss type function.
    //  1  trapezoid type function.
    //  2  Trigonometric functions, sin^2.
    //  3  heaviside function.
    //  4  HHG function.

    int tstart;
    int tend;
    double dt;

    // space domain parameters

    //length gauge
    double lcut1;
    double lcut2;

    // time domain parameters

    // Gauss
    int gauss_count;
    std::vector<double> gauss_omega; // time(a.u.)^-1 
    std::vector<double> gauss_phase ;
    std::vector<double> gauss_sigma ; // time(a.u.)
    std::vector<double> gauss_t0 ;
    std::vector<double> gauss_amp ;  // Ry/bohr

    // trapezoid
    int trape_count;
    std::vector<double> trape_omega ; // time(a.u.)^-1 
    std::vector<double> trape_phase ;
    std::vector<double> trape_t1 ;
    std::vector<double> trape_t2 ;
    std::vector<double> trape_t3 ;
    std::vector<double> trape_amp ; // Ry/bohr

    // Trigonometric
    int trigo_count;
    std::vector<double> trigo_omega1 ; // time(a.u.)^-1 
    std::vector<double> trigo_omega2 ; // time(a.u.)^-1 
    std::vector<double> trigo_phase1 ;
    std::vector<double> trigo_phase2 ;
    std::vector<double> trigo_amp ; // Ry/bohr

    // Heaviside
    int heavi_count;
    std::vector<double> heavi_t0;
    std::vector<double> heavi_amp; // Ry/bohr

    // HHG
    // int hhg_count = 0;
    // std::vector<double> hhg_amp1; // Ry/bohr
    // std::vector<double> hhg_amp2; // Ry/bohr
    // std::vector<double> hhg_omega1; // time(a.u.)^-1 
    // std::vector<double> hhg_omega2; // time(a.u.)^-1 
    // std::vector<double> hhg_phase1;
    // std::vector<double> hhg_phase2;
    // std::vector<double> hhg_t0;
    // std::vector<double> hhg_sigma; // time(a.u.)

    const UnitCell* ucell_ = nullptr;

    std::vector<double> set_parameters(std::string params, double c);
    void read_parameters(Input* in);

    // potential of electric field in space domain : length gauge and velocity gauge
    void cal_v_space(std::vector<double> &vext_space, int direc);
    void cal_v_space_length(std::vector<double> &vext_space, int direc);
    double cal_v_space_length_potential(double i);
    void cal_v_space_velocity(std::vector<double> &vext_space, int direc);

    // potential of electric field in time domain : Gauss , trapezoid, trigonometric, heaviside,  HHG
    double cal_v_time(int t_type);
    double cal_v_time_Gauss();
    double cal_v_time_trapezoid();
    double cal_v_time_trigonometric();
    double cal_v_time_heaviside();
    // double cal_v_time_HHG();
};

} // namespace elecstate

#endif