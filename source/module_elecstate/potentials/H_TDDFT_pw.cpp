#include "module_base/timer.h"
#include "potential_new.h"
#include "src_lcao/ELEC_evolve.h"
#include "src_pw/global.h"

namespace elecstate
{

//==========================================================
// this function aims to add external time-dependent potential
// (eg: linear potential) used in tddft
// fuxiang add in 2022-10
//==========================================================

int vext_spatial_domain = 1;
// Here vext_spatial_domain =1 means length gauge and vext_spatial_domain = 2 means velocity gauge in spatial domain.
int vext_time_domain = 1;
// Here vext_time_domain means different types of time functions.
// vext_time_domain = 1  Gauss type function.
// vext_time_domain = 2  trapezoid type function.
// vext_time_domain = 3  Trigonometric functions, sin^2.

int timescale_start = 1; // get the time start;
int timescale_end = 1000; // get the time end;

double td_dt = 0.05; // dt fs

// Gauss
double vext_gauss_omega = 22.13; // eV
double vext_gauss_sigmasquare = 700;
double vext_gauss_timecenter;
double vext_gauss_amplitude = 0.25;

// trapezoid
double vext_trapezoid_omega = 1.60; // eV
double vext_trapezoid_stepcut1 = 1875;
double vext_trapezoid_stepcut2 = 5625;
double vext_trapezoid_stepcut3 = 7500;
double vext_trapezoid_amplitude = 2.74;

// Trigonometric
double vext_trigonometric_omega1 = 1.164656; // eV
double vext_trigonometric_omega2 = 0.029116; // eV
double vext_trigonometric_amplitude = 2.74;

// Heaviside
int vext_heaviside_switch;
double vext_heaviside_amplitude = 1.0;

// vext_spatial_domain = 1;
// vext_time_domain = 1;
// vext_gauss_omega = 22.13;
// td_dt = 0.05;
// vext_gauss_amplitude = 0.25;

// vext_trapezoid_omega = 1.60;
// vext_trapezoid_stepcut1 = 1875;
// vext_trapezoid_stepcut2 = 5625;
// vext_trapezoid_stepcut3 = 7500;
// vext_trapezoid_amplitude = 2.74;

// vext_trigonometric_omega1 = 1.164656;
// vext_trigonometric_omega2 = 0.029116;
// vext_trigonometric_amplitude = 2.74;
// vext_heaviside_amplitude = 1.0;

void Potential::update_for_tddft(int istep)
{
    ModuleBase::TITLE("Potential", "update_for_tddft");
    ModuleBase::timer::tick("Potential", "update_for_tddft");

    if (istep > timescale_end || istep < timescale_start)
    {
        std::cout << "vext = 0! " << std::endl;
        return;
    }

    for (int is = 0; is < GlobalV::NSPIN; is++)
    {
        //====================================================
        // add external linear potential, fuxiang add in 2017/05
        //====================================================
        double* vextold = new double[GlobalC::rhopw->nrxx];
        double* vext = new double[GlobalC::rhopw->nrxx];
        const int yz = GlobalC::rhopw->ny * GlobalC::rhopw->nplane;
        int index, i, j, k;

        if (vext_spatial_domain == 1)
        {
            for (int ir = 0; ir < GlobalC::rhopw->nrxx; ++ir)
            {
                index = ir;
                i = index / yz; // get the z, z is the fastest
                index = index - yz * i; // get (x,y)
                j = index / GlobalC::rhopw->nplane; // get y
                k = index - GlobalC::rhopw->nplane * j + GlobalC::rhopw->startz_current; // get x

                if (ELEC_evolve::td_vext_dire == 1)
                {
                    if (k < GlobalC::rhopw->nx * 0.05)
                    {
                        vextold[ir] = (0.019447 * k / GlobalC::rhopw->nx - 0.001069585) * GlobalC::ucell.lat0;
                        // here 0.019447, 0.001069585 means
                    }
                    else if (k >= GlobalC::rhopw->nx * 0.05 && k < GlobalC::rhopw->nx * 0.95)
                    {
                        vextold[ir] = -0.0019447 * k / GlobalC::rhopw->nx * GlobalC::ucell.lat0;
                    }
                    else if (k >= GlobalC::rhopw->nx * 0.95)
                    {
                        vextold[ir]
                            = (0.019447 * (1.0 * k / GlobalC::rhopw->nx - 1) - 0.001069585) * GlobalC::ucell.lat0;
                    }
                }
                else if (ELEC_evolve::td_vext_dire == 2)
                {
                    if (j < GlobalC::rhopw->nx * 0.05)
                    {
                        vextold[ir] = (0.019447 * j / GlobalC::rhopw->nx - 0.001069585) * GlobalC::ucell.lat0;
                    }
                    else if (j >= GlobalC::rhopw->nx * 0.05 && j < GlobalC::rhopw->nx * 0.95)
                    {
                        vextold[ir] = -0.0019447 * j / GlobalC::rhopw->nx * GlobalC::ucell.lat0;
                    }
                    else if (j >= GlobalC::rhopw->nx * 0.95)
                    {
                        vextold[ir]
                            = (0.019447 * (1.0 * j / GlobalC::rhopw->nx - 1) - 0.001069585) * GlobalC::ucell.lat0;
                    }
                }
                else if (ELEC_evolve::td_vext_dire == 3)
                {
                    if (i < GlobalC::rhopw->nx * 0.05)
                    {
                        vextold[ir] = (0.019447 * i / GlobalC::rhopw->nx - 0.001069585) * GlobalC::ucell.lat0;
                    }
                    else if (i >= GlobalC::rhopw->nx * 0.05 && i < GlobalC::rhopw->nx * 0.95)
                    {
                        vextold[ir] = -0.0019447 * i / GlobalC::rhopw->nx * GlobalC::ucell.lat0;
                    }
                    else if (i >= GlobalC::rhopw->nx * 0.95)
                    {
                        vextold[ir]
                            = (0.019447 * (1.0 * i / GlobalC::rhopw->nx - 1) - 0.001069585) * GlobalC::ucell.lat0;
                    }
                }

                // Gauss
                if (vext_time_domain == 1)
                {
                    // const double w = 22.13;    // eV  //
                    // const double sigmasquare = 700;
                    // const double timecenter = 700;
                    // Notice: these three parameters should be written in INPUT. I will correct soon.

                    double vext_gauss_timenow
                        = (istep - vext_gauss_timecenter) * td_dt * 41.34; // 41.34 is conversion factor of fs-a.u.
                    vext[ir] = vextold[ir] * cos(vext_gauss_omega / 27.2116 * vext_gauss_timenow)
                                     * exp(-vext_gauss_timenow * vext_gauss_timenow * 0.5 / (vext_gauss_sigmasquare))
                                     * vext_gauss_amplitude;
                }

                // HHG of H atom
                if (vext_time_domain == 2)
                {
                    // const double w_h = 0.0588; //a.u.
                    double w_h = vext_trapezoid_omega / 27.2116;
                    // const int stepcut1 = 1875;
                    // const int stepcut2 = 5625;
                    // const int stepcut3 = 7500;
                    //  The parameters should be written in INPUT!
                    if (istep < vext_trapezoid_stepcut1)
                    {
                        vext[ir] = vextold[ir] * vext_trapezoid_amplitude * istep / vext_trapezoid_stepcut1
                                         * cos(w_h * istep * td_dt * 41.34);
                    }
                    else if (istep < vext_trapezoid_stepcut2)
                    {
                        vext[ir]
                            = vextold[ir] * vext_trapezoid_amplitude * cos(w_h * istep * td_dt * 41.34);
                    }
                    else if (istep < vext_trapezoid_stepcut3)
                    {
                        vext[ir] = vextold[ir] * vext_trapezoid_amplitude
                                         * (vext_trapezoid_stepcut3 - istep) / vext_trapezoid_stepcut1
                                         * cos(w_h * istep * td_dt * 41.34);
                    }
                }

                // HHG of H2
                //  Type 3 will be modified into more normolized type soon.
                if (vext_time_domain == 3)
                {
                    const double w_h2 = vext_trigonometric_omega1 / 27.2116; // a.u.
                    const double w_h3 = vext_trigonometric_omega1 / 27.2116; // a.u.
                    const double timenow = (istep)*td_dt * 41.34;
                    // The parameters should be written in INPUT!

                    // vext[ir] =
                    // vextold[ir]*2.74*cos(0.856*timenow)*sin(0.0214*timenow)*sin(0.0214*timenow);
                    // vext[ir] =
                    // vextold[ir]*2.74*cos(0.856*timenow)*sin(0.0214*timenow)*sin(0.0214*timenow)*0.01944;
                    vext[ir] = vextold[ir] * vext_trigonometric_amplitude * cos(w_h2 * timenow)
                                     * sin(w_h3 * timenow) * sin(w_h3 * timenow);
                }

                // Heaviside function
                if (vext_time_domain == 4)
                {
                    if (istep < vext_heaviside_switch)
                        vext[ir] = vextold[ir] * vext_heaviside_amplitude;
                    else if (istep >= vext_heaviside_switch)
                        vext[ir] = vextold[ir] * 0;
                }

                this->v_effective(is, ir) += vext[ir];

                // std::cout << "x: " << k <<"	" << "y: " << j <<"	"<< "z: "<< i <<"	"<< "ir: " << ir <<
                // std::endl; std::cout << "vext: " << vext[ir] << std::endl; std::cout << "vrs: " <<
                // vrs(is,ir) <<std::endl;
            }
        }

        if (vext_spatial_domain == 2)
        {
            for (int ir = 0; ir < GlobalC::rhopw->nrxx; ++ir)
            {
                index = ir;
                i = index / yz; // get the z, z is the fastest
                index = index - yz * i; // get (x,y)
                j = index / GlobalC::rhopw->nplane; // get y
                k = index - GlobalC::rhopw->nplane * j + GlobalC::rhopw->startz_current; // get x

                // vextold[ir] = vext_S[ir] + vext_A2S[ir] + vext_Ar[ir];

                // vext_S = new double[GlobalC::rhopw->nrxx];
                // vext_A2S = new double[GlobalC::rhopw->nrxx];

                // Calculate A(t)
                double A_t;
                A_t = 0.0;

                for (int it = 0; it < istep; ++it)
                {
                    if (vext_time_domain == 1)
                    {
                        double vext_gauss_timenow = (it - vext_gauss_timecenter) * td_dt * 41.34;
                        A_t += -137.124 * cos(vext_gauss_omega / 27.2116 * vext_gauss_timenow)
                               * exp(-vext_gauss_timenow * vext_gauss_timenow * 0.5 / (vext_gauss_sigmasquare))
                               * vext_gauss_amplitude * td_dt * 41.34;
                        // -137.124 is Atomic unit of the speed of light.
                    }
                    if (vext_time_domain == 2)
                    {
                        double w_h = vext_trapezoid_omega / 27.2116;
                        if (it < vext_trapezoid_stepcut1)
                        {
                            A_t += -137.124 * vext_trapezoid_amplitude * it / vext_trapezoid_stepcut1
                                   * cos(w_h * it * td_dt * 41.34) * td_dt * 41.34;
                        }
                        else if (it < vext_trapezoid_stepcut2)
                        {
                            A_t += -137.124 * vext_trapezoid_amplitude * cos(w_h * it * td_dt * 41.34) * td_dt * 41.34;
                        }
                        else if (it < vext_trapezoid_stepcut3)
                        {
                            A_t += vext_trapezoid_amplitude * (vext_trapezoid_stepcut3 - it) / vext_trapezoid_stepcut1
                                   * cos(w_h * it * td_dt * 41.34) * td_dt * 41.34;
                        }
                    }
                    if (vext_time_domain == 3)
                    {
                        const double w_h2 = vext_trigonometric_omega1 / 27.2116; // a.u.
                        const double w_h3 = vext_trigonometric_omega1 / 27.2116; // a.u.
                        const double timenow = (it)*td_dt * 41.34;
                        A_t += vext_trigonometric_amplitude * cos(w_h2 * timenow) * sin(w_h3 * timenow)
                               * sin(w_h3 * timenow) * td_dt * 41.34;
                    }

                    if (vext_time_domain == 4)
                    {
                        if (it < vext_heaviside_switch)
                            A_t += vext_heaviside_amplitude * td_dt * 41.34;
                        else if (istep >= vext_heaviside_switch)
                            A_t += 0;
                    }
                }

                // vextold[ir] = vext_S[ir] + vext_A2S[ir];
                for (int it = 0; it < GlobalC::rhopw->nrxx; ++it)
                {
                    vextold[it] = A_t * A_t;
                }
                // Another term <psi|A▽|psi> will be supplemented at the kinetic energy term.

                // time domian

                // Gauss
                if (vext_time_domain == 1)
                {
                    double vext_gauss_timenow
                        = (istep - vext_gauss_timecenter) * td_dt * 41.34; // 41.34 is conversion factor of fs-a.u.
                    vext[ir] = vextold[ir] * cos(vext_gauss_omega / 27.2116 * vext_gauss_timenow)
                                     * exp(-vext_gauss_timenow * vext_gauss_timenow * 0.5 / (vext_gauss_sigmasquare))
                                     * vext_gauss_amplitude;
                }

                // HHG of H atom
                if (vext_time_domain == 2)
                {
                    // const double w_h = 0.0588; //a.u.
                    double w_h = vext_trapezoid_omega / 27.2116;

                    if (istep < vext_trapezoid_stepcut1)
                    {
                        vext[ir] = vextold[ir] * vext_trapezoid_amplitude * istep / vext_trapezoid_stepcut1
                                         * cos(w_h * istep * td_dt * 41.34);
                    }
                    else if (istep < vext_trapezoid_stepcut2)
                    {
                        vext[ir]
                            = vextold[ir] * vext_trapezoid_amplitude * cos(w_h * istep * td_dt * 41.34);
                    }
                    else if (istep < vext_trapezoid_stepcut3)
                    {
                        vext[ir] = vextold[ir] * vext_trapezoid_amplitude
                                         * (vext_trapezoid_stepcut3 - istep) / vext_trapezoid_stepcut1
                                         * cos(w_h * istep * td_dt * 41.34);
                    }
                }

                // HHG of H2
                if (vext_time_domain == 3)
                {
                    const double w_h2 = vext_trigonometric_omega1 / 27.2116; // a.u.
                    const double w_h3 = vext_trigonometric_omega1 / 27.2116; // a.u.
                    const double timenow = (istep)*td_dt * 41.34;
                    // The parameters should be written in INPUT!

                    // vext[ir] =
                    // vextold[ir]*2.74*cos(0.856*timenow)*sin(0.0214*timenow)*sin(0.0214*timenow);
                    // vext[ir] =
                    // vextold[ir]*2.74*cos(0.856*timenow)*sin(0.0214*timenow)*sin(0.0214*timenow)*0.01944;
                    vext[ir] = vextold[ir] * vext_trigonometric_amplitude * cos(w_h2 * timenow)
                                     * sin(w_h3 * timenow) * sin(w_h3 * timenow);
                }

                // Heaviside function
                if (vext_time_domain == 4)
                {
                    if (istep < vext_heaviside_switch)
                        vext[ir] = vextold[ir] * vext_heaviside_amplitude;
                    else if (istep >= vext_heaviside_switch)
                        vext[ir] = vextold[ir] * 0;
                }

                this->v_effective(is, ir) += vext[ir];

                // std::cout << "x: " << k <<"	" << "y: " << j <<"	"<< "z: "<< i <<"	"<< "ir: " << ir <<
                // std::endl; std::cout << "vext: " << vext[ir] << std::endl; std::cout << "vrs: " <<
                // vrs(is,ir) <<std::endl;
            }
        }

        std::cout << "vext exists" << std::endl;

        delete[] vextold;
        delete[] vext;
    }

    ModuleBase::timer::tick("potential", "update_for_tddft");
    return;
} // end subroutine set_vrs_tddft

} // namespace elecstate
