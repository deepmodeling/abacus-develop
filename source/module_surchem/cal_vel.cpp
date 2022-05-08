#include "../module_base/timer.h"
#include "../module_xc/xc_functional.h"
#include "surchem.h"

void shape_gradn(const double *PS_TOTN_real, PW_Basis &pwb, double *eprime)
{

    double epr_c = 1.0 / sqrt(ModuleBase::TWO_PI) / GlobalV::sigma_k;
    double epr_z = 0;
    double min = 1e-10;
    for (int ir = 0; ir < pwb.nrxx; ir++)
    {
        epr_z = log(max(PS_TOTN_real[ir], min) / GlobalV::nc_k) / sqrt(2) / GlobalV::sigma_k;
        eprime[ir] = epr_c * exp(-pow(epr_z, 2)) / max(PS_TOTN_real[ir], min);
    }
}

void eps_pot(const double *PS_TOTN_real, const complex<double> *phi, PW_Basis &pwb, double *d_eps, double *vwork)
{
    double *eprime = new double[pwb.nrxx];
    ModuleBase::GlobalFunc::ZEROS(eprime, pwb.nrxx);

    shape_gradn(PS_TOTN_real, pwb, eprime);

    for (int ir = 0; ir < pwb.nrxx; ir++)
    {
        eprime[ir] = eprime[ir] * (GlobalV::eb_k - 1);
    }

    ModuleBase::Vector3<double> *nabla_phi = new ModuleBase::Vector3<double>[pwb.nrxx];
    double *phisq = new double[pwb.nrxx];

    // nabla phi
    XC_Functional::grad_rho(phi, nabla_phi);

    for (int ir = 0; ir < pwb.nrxx; ir++)
    {
        phisq[ir] = pow(nabla_phi[ir].x, 2) + pow(nabla_phi[ir].y, 2) + pow(nabla_phi[ir].z, 2);
    }

    for (int ir = 0; ir < pwb.nrxx; ir++)
    {
        vwork[ir] = eprime[ir] * phisq[ir] / (8 * ModuleBase::PI);
    }

    delete[] eprime;
    delete[] nabla_phi;
    delete[] phisq;
}

ModuleBase::matrix surchem::cal_vel(const UnitCell &cell,
                                    PW_Basis &pwb,
                                    const complex<double> *TOTN,
                                    const complex<double> *PS_TOTN,
                                    int nspin)
{
    ModuleBase::TITLE("surchem", "cal_vel");
    ModuleBase::timer::tick("surchem", "cal_vel");

    double *TOTN_real = new double[pwb.nrxx];
    GlobalC::UFFT.ToRealSpace(TOTN, TOTN_real);

    // -4pi * TOTN(G)
    complex<double> *B = new complex<double>[pwb.ngmc];
    for (int ig = 0; ig < pwb.ngmc; ig++)
    {
        B[ig] = -4.0 * ModuleBase::PI * TOTN[ig];
    }

    // Build a nrxx vector to DO FFT .
    double *PS_TOTN_real = new double[pwb.nrxx];
    GlobalC::UFFT.ToRealSpace(PS_TOTN, PS_TOTN_real);

    // build epsilon in real space (nrxx)
    double *epsilon = new double[pwb.nrxx];
    double *epsilon0 = new double[pwb.nrxx];
    cal_epsilon(pwb, PS_TOTN_real, epsilon, epsilon0);

    complex<double> *Sol_phi = new complex<double>[pwb.ngmc];
    complex<double> *Sol_phi0 = new complex<double>[pwb.ngmc];
    int ncgsol = 0;

    double *Vel = new double[pwb.nrxx];
    ModuleBase::GlobalFunc::ZEROS(Vel, pwb.nrxx);

    // Calculate Sol_phi with epsilon.
    ncgsol = 0;
    minimize_cg(cell, pwb, epsilon, B, Sol_phi, ncgsol);

    ncgsol = 0;
    // Calculate Sol_phi0 with epsilon0.
    minimize_cg(cell, pwb, epsilon0, B, Sol_phi0, ncgsol);

    double *phi_tilda_R = new double[pwb.nrxx];
    double *phi_tilda_R0 = new double[pwb.nrxx];
    double *delta_phi_R = new double[pwb.nrxx];

    GlobalC::UFFT.ToRealSpace(Sol_phi, phi_tilda_R);
    GlobalC::UFFT.ToRealSpace(Sol_phi0, phi_tilda_R0);

    // the 1st item of Vel
    for (int i = 0; i < pwb.nrxx; i++)
    {
        delta_phi_R[i] = phi_tilda_R[i] - phi_tilda_R0[i];
        Vel[i] += delta_phi_R[i];
    }

    // calculate Ael
    double Ael = cal_Ael(cell, pwb, TOTN_real, delta_phi_R);

    // the 2nd item of Vel
    double *Vel2 = new double[pwb.nrxx];
    ModuleBase::GlobalFunc::ZEROS(Vel2, pwb.nrxx);

    eps_pot(PS_TOTN_real, Sol_phi, pwb, epsilon, Vel2);

    for (int i = 0; i < pwb.nrxx; i++)
    {
        Vel[i] += Vel2[i];
    }

    ModuleBase::matrix v(nspin, pwb.nrxx);

    if (nspin == 4)
    {
        for (int ir = 0; ir < pwb.nrxx; ir++)
        {
            v(0, ir) += Vel[ir];
        }
    }
    else
    {
        for (int is = 0; is < nspin; is++)
        {
            for (int ir = 0; ir < pwb.nrxx; ir++)
            {
                v(is, ir) += Vel[ir];
            }
        }
    }

    delete[] PS_TOTN_real;
    delete[] Sol_phi;
    delete[] Sol_phi0;
    delete[] B;
    delete[] epsilon;
    delete[] epsilon0;
    delete[] Vel;
    delete[] Vel2;
    delete[] TOTN_real;
    delete[] phi_tilda_R;
    delete[] phi_tilda_R0;
    delete[] delta_phi_R;

    ModuleBase::timer::tick("surchem", "cal_vel");
    return v;
}
