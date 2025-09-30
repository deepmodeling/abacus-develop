#ifndef H_TDDFT_PW_H
#define H_TDDFT_PW_H

#include "pot_base.h"
#include "source_io/input_conv.h"
#include "source_io/module_parameter/parameter.h" // PARAM.globalv.global_readin_dir, PARAM.inp.mdp.md_restart

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

        // If it is the first time to create an H_TDDFT_pw instance and is restart calculation,
        // initialize istep using current_step_info
        if (!is_initialized && PARAM.inp.mdp.md_restart)
        {
            int restart_istep = -1;
            std::string file_dir = PARAM.globalv.global_readin_dir;
            current_step_info(file_dir, restart_istep);

            if (restart_istep >= 0)
            {
                H_TDDFT_pw::istep = restart_istep - 1; // Update istep
            }

            is_initialized = true; // Mark as initialized, so that istep will not be initialized again
        }
    }

    ~H_TDDFT_pw() {};

    void cal_fixed_v(double* vl_pseudo) override;

    /**
     * @brief Compute ionic force of electric field
     *
     * @param[in] cell Information of cell
     * @param[out] fe Force of electric field  F = qE
     */
    static void compute_force(const UnitCell& cell, ModuleBase::matrix& fe);

    // parameters
    static int stype; // 0: length gauge; 1: velocity gauge; 2: hybrid gauge

    static std::vector<int> ttype;
    //  0: Gaussian type function.
    //  1: Trapezoid type function.
    //  2: Trigonometric functions, sin^2.
    //  3: Heaviside step function.

    static int tstart;
    static int tend;
    static double dt;
    // cut dt for integral
    static double dt_int;
    static int istep_int;

    // Space domain parameters

    // length gauge
    static double lcut1;
    static double lcut2;

    // velocity gauge, vector potential
    static ModuleBase::Vector3<double> At;
    static ModuleBase::Vector3<double> At_laststep;
    static ModuleBase::Vector3<double> Et;

    // Time domain parameters

    // Gauss
    static int gauss_count;
    static std::vector<double> gauss_omega; // time(a.u.)^-1
    static std::vector<double> gauss_phase;
    static std::vector<double> gauss_sigma; // time(a.u.)
    static std::vector<double> gauss_t0;
    static std::vector<double> gauss_amp; // Ry/bohr
    // add for velocity gauge, recut dt into n pieces to make sure the integral is accurate enough
    // must be even, thus would get odd number of points for Simpson integral
    static std::vector<int> gauss_ncut;

    // Trapezoid
    static int trape_count;
    static std::vector<double> trape_omega; // time(a.u.)^-1
    static std::vector<double> trape_phase;
    static std::vector<double> trape_t1;
    static std::vector<double> trape_t2;
    static std::vector<double> trape_t3;
    static std::vector<double> trape_amp; // Ry/bohr
    // add for velocity gauge, recut dt into n pieces to make sure the integral is accurate enough
    static std::vector<int> trape_ncut;

    // Trigonometric
    static int trigo_count;
    static std::vector<double> trigo_omega1; // time(a.u.)^-1
    static std::vector<double> trigo_omega2; // time(a.u.)^-1
    static std::vector<double> trigo_phase1;
    static std::vector<double> trigo_phase2;
    static std::vector<double> trigo_amp; // Ry/bohr
    // add for velocity gauge, recut dt into n pieces to make sure the integral is accurate enough
    static std::vector<int> trigo_ncut;

    // Heaviside
    static int heavi_count;
    static std::vector<double> heavi_t0;
    static std::vector<double> heavi_amp; // Ry/bohr

    // update At for velocity gauge by intergral of E(t)dt
    static void update_At();

  private:
    static int istep;
    static bool is_initialized; // static flag variable, used to ensure initialization only once

    static double amp;
    static vector<double> global_vext_time;

    const UnitCell* ucell_ = nullptr;

    // Obtain the current MD step information, used for restart calculation
    void current_step_info(const std::string& file_dir, int& istep);

    // Potential of electric field in space domain: for length gauge only
    void cal_v_space(std::vector<double>& vext_space, int direc);
    void cal_v_space_length(std::vector<double>& vext_space, int direc);
    double cal_v_space_length_potential(double i);

    // Potential of electric field in time domain: Gaussian, trapezoid, trigonometric, Heaviside
    static double cal_v_time(int t_type, const bool last);
    static double cal_v_time_Gauss(const bool last);
    static double cal_v_time_trapezoid(const bool last);
    static double cal_v_time_trigonometric(const bool last);
    static double cal_v_time_heaviside(const bool last);

    // Get ncut number for At integral
    static int check_ncut(int t_type);
};

} // namespace elecstate

#endif
