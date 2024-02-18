#ifndef H_TDDFT_PW_H
#define H_TDDFT_PW_H

#include "module_io/input.h"
#include "module_io/input_conv.h"
#include "pot_base.h"

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

    /**
     * @brief compute force of electric field
     *
     * @param[in] cell information of cell
     * @param[out] fe force of electric field  F=qE
     */
    static void compute_force(const UnitCell& cell, ModuleBase::matrix& fe);

    // parameters
    static int stype; // 0 : length gauge  1: velocity gauge

    static std::vector<int> ttype;
    //  0  Gauss type function.
    //  1  trapezoid type function.
    //  2  Trigonometric functions, sin^2.
    //  3  heaviside function.
    //  4  HHG function.

    static int tstart;
    static int tend;
    static double dt;
    //cut dt for integral
    static double dt_int;
    static int istep_int;

    // space domain parameters

    //length gauge
    static double lcut1;
    static double lcut2;

    //velocity gauge, vector magnetic potential
    static double At[3];

    // time domain parameters

    // Gauss
    static int gauss_count;
    static std::vector<double> gauss_omega; // time(a.u.)^-1
    static std::vector<double> gauss_phase;
    static std::vector<double> gauss_sigma; // time(a.u.)
    static std::vector<double> gauss_t0;
    static std::vector<double> gauss_amp; // Ry/bohr
    //add for velocity gauge, recut dt into n pieces to make sure the integral is accurate Enough
    //must be even, thus would get odd number of points for simpson integral
    static std::vector<int> gauss_ncut;

    // trapezoid
    static int trape_count;
    static std::vector<double> trape_omega; // time(a.u.)^-1
    static std::vector<double> trape_phase;
    static std::vector<double> trape_t1;
    static std::vector<double> trape_t2;
    static std::vector<double> trape_t3;
    static std::vector<double> trape_amp; // Ry/bohr
    //add for velocity gauge, recut dt into n pieces to make sure the integral is accurate Enough
    static std::vector<int> trape_ncut;

    // Trigonometric
    static int trigo_count;
    static std::vector<double> trigo_omega1; // time(a.u.)^-1
    static std::vector<double> trigo_omega2; // time(a.u.)^-1
    static std::vector<double> trigo_phase1;
    static std::vector<double> trigo_phase2;
    static std::vector<double> trigo_amp; // Ry/bohr
    //add for velocity gauge, recut dt into n pieces to make sure the integral is accurate Enough
    static std::vector<int> trigo_ncut;

    // Heaviside
    static int heavi_count;
    static std::vector<double> heavi_t0;
    static std::vector<double> heavi_amp; // Ry/bohr

    //update At for velocity gauge by intergral of E(t)dt
    static void update_At(void);

  private:
    // internal time-step,
    //-------hypothesis-------
    // Vext will evolve by time, every time cal_fixed_v() is called, istep++
    //------------------------
    static int istep;

    static double amp;

    static double bmod;
    static double bvec[3];

    const UnitCell* ucell_ = nullptr;

    // potential of electric field in space domain : length gauge and velocity gauge
    void cal_v_space(std::vector<double> &vext_space, int direc);
    void cal_v_space_length(std::vector<double> &vext_space, int direc);
    double cal_v_space_length_potential(double i);

    // potential of electric field in time domain : Gauss , trapezoid, trigonometric, heaviside,  HHG
    static double cal_v_time(int t_type, const bool last);
    static double cal_v_time_Gauss(const bool last);
    static double cal_v_time_trapezoid(const bool last);
    static double cal_v_time_trigonometric(const bool last);
    static double cal_v_time_heaviside();
    // double cal_v_time_HHG();

    //get ncut number for At integral
    static int check_ncut(int t_type);

    void prepare(const UnitCell& cell, int& dir);
};

} // namespace elecstate

#endif