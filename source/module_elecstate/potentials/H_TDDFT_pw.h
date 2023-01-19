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

    int ttype ;
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
    double gauss_omega; // time(a.u.)^-1 
    double gauss_phase ;
    double gauss_sigma2 ; // time(a.u.)^2
    double gauss_t0 ;
    double gauss_amp ;  // Ry/bohr

    // trapezoid
    double trape_omega ; // time(a.u.)^-1 
    double trape_phase ;
    double trape_t1 ;
    double trape_t2 ;
    double trape_t3 ;
    double trape_amp ; // Ry/bohr

    // Trigonometric
    double trigo_omega1 ; // time(a.u.)^-1 
    double trigo_omega2 ; // time(a.u.)^-1 
    double trigo_phase1 ;
    double trigo_phase2 ;
    double trigo_amp ; // Ry/bohr

    // Heaviside
    int heavi_t0;
    double heavi_amp; // Ry/bohr

    // HHG
    double hhg_amp1; // Ry/bohr
    double hhg_amp2; // Ry/bohr
    double hhg_omega1; // time(a.u.)^-1 
    double hhg_omega2; // time(a.u.)^-1 
    double hhg_phase1;
    double hhg_phase2;
    double hhg_t0;
    double hhg_sigma2; // time(a.u.)^2

    const UnitCell* ucell_ = nullptr;

    void read_parameters(Input* in);

    // potential of electric field in space domain : length gauge and velocity gauge
    void cal_v_space(double* vext_space);
    void cal_v_space_length(double* vext_space);
    double cal_v_space_length_potential(double i);
    void cal_v_space_velocity(double* vext_space);

    // potential of electric field in time domain : Gauss , trapezoid, trigonometric, heaviside,  HHG
    double cal_v_time();
    double cal_v_time_Gauss();
    double cal_v_time_trapezoid();
    double cal_v_time_trigonometric();
    double cal_v_time_heaviside();
    double cal_v_time_HHG();
};

} // namespace elecstate

#endif