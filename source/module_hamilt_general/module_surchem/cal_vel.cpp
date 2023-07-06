#include "module_base/timer.h"
#include "module_hamilt_general/module_xc/xc_functional.h"
#include "surchem.h"

void shape_gradn(const double* PS_TOTN_real, const ModulePW::PW_Basis* rho_basis, double* eprime)
{

    double epr_c = 1.0 / sqrt(ModuleBase::TWO_PI) / GlobalV::sigma_k;
    double epr_z = 0;
    double min = 1e-10;
    for (int ir = 0; ir < rho_basis->nrxx; ir++)
    {
        epr_z = log(std::max(PS_TOTN_real[ir], min) / GlobalV::nc_k) / sqrt(2) / GlobalV::sigma_k;
        eprime[ir] = epr_c * exp(-pow(epr_z, 2)) / std::max(PS_TOTN_real[ir], min);
    }
}

void eps_pot(const double* PS_TOTN_real,
             const complex<double>* phi,
             const ModulePW::PW_Basis* rho_basis,
             double* d_eps,
             double* vwork)
{
    double *eprime = new double[rho_basis->nrxx];
    ModuleBase::GlobalFunc::ZEROS(eprime, rho_basis->nrxx);

    shape_gradn(PS_TOTN_real, rho_basis, eprime);

    for (int ir = 0; ir < rho_basis->nrxx; ir++)
    {
        eprime[ir] = eprime[ir] * (GlobalV::eb_k - 1);
    }

    ModuleBase::Vector3<double> *nabla_phi = new ModuleBase::Vector3<double>[rho_basis->nrxx];
    double *phisq = new double[rho_basis->nrxx];

    // nabla phi
    XC_Functional::grad_rho(phi, nabla_phi, rho_basis, GlobalC::ucell.tpiba);

    for (int ir = 0; ir < rho_basis->nrxx; ir++)
    {
        phisq[ir] = pow(nabla_phi[ir].x, 2) + pow(nabla_phi[ir].y, 2) + pow(nabla_phi[ir].z, 2);
    }

    for (int ir = 0; ir < rho_basis->nrxx; ir++)
    {
        vwork[ir] = eprime[ir] * phisq[ir] / (8 * ModuleBase::PI);
    }

    delete[] eprime;
    delete[] nabla_phi;
    delete[] phisq;
}

ModuleBase::matrix surchem::cal_vel(const UnitCell& cell,
                                    const ModulePW::PW_Basis* rho_basis,
                                    complex<double>* TOTN,
                                    complex<double>* PS_TOTN,
                                    int nspin)
{
    ModuleBase::TITLE("surchem", "cal_vel");
    ModuleBase::timer::tick("surchem", "cal_vel");

    rho_basis->recip2real(TOTN, TOTN_real);

    // -4pi * TOTN(G)
    complex<double> *B = new complex<double>[rho_basis->npw];
    for (int ig = 0; ig < rho_basis->npw; ig++)
    {
        B[ig] = -4.0 * ModuleBase::PI * TOTN[ig];
    }

    // Build a nrxx vector to DO FFT .
    double *PS_TOTN_real = new double[rho_basis->nrxx];
    rho_basis->recip2real(PS_TOTN, PS_TOTN_real);

    // build epsilon in real space (nrxx)
    double *epsilon = new double[rho_basis->nrxx];
    double *epsilon0 = new double[rho_basis->nrxx];
    cal_epsilon(rho_basis, PS_TOTN_real, epsilon, epsilon0);

    complex<double> *Sol_phi = new complex<double>[rho_basis->npw];
    complex<double> *Sol_phi0 = new complex<double>[rho_basis->npw];
    int ncgsol = 0;

    double *tmp_Vel = new double[rho_basis->nrxx];
    ModuleBase::GlobalFunc::ZEROS(tmp_Vel, rho_basis->nrxx);

    // Calculate Sol_phi with epsilon.
    ncgsol = 0;
    minimize_cg(cell, rho_basis, epsilon, B, Sol_phi, ncgsol);

    ncgsol = 0;
    // Calculate Sol_phi0 with epsilon0.
    minimize_cg(cell, rho_basis, epsilon0, B, Sol_phi0, ncgsol);

    double *phi_tilda_R = new double[rho_basis->nrxx];
    double *phi_tilda_R0 = new double[rho_basis->nrxx];

    rho_basis->recip2real(Sol_phi, phi_tilda_R);
    rho_basis->recip2real(Sol_phi0, phi_tilda_R0);

    // the 1st item of tmp_Vel
    for (int i = 0; i < rho_basis->nrxx; i++)
    {
        delta_phi[i] = phi_tilda_R[i] - phi_tilda_R0[i];
        tmp_Vel[i] += delta_phi[i];
    }

    // calculate Ael
    // double Ael = cal_Ael(cell, pwb);

    // the 2nd item of tmp_Vel
    eps_pot(PS_TOTN_real, Sol_phi, rho_basis, epsilon, epspot);

    for (int i = 0; i < rho_basis->nrxx; i++)
    {
        tmp_Vel[i] += epspot[i];
    }

    // ModuleBase::matrix v(nspin, rho_basis->nrxx);
    ModuleBase::GlobalFunc::ZEROS(Vel.c, nspin * rho_basis->nrxx);

    if (nspin == 4)
    {
        for (int ir = 0; ir < rho_basis->nrxx; ir++)
        {
            Vel(0, ir) += tmp_Vel[ir];
        }
    }
    else
    {
        for (int is = 0; is < nspin; is++)
        {
            for (int ir = 0; ir < rho_basis->nrxx; ir++)
            {
                Vel(is, ir) += tmp_Vel[ir];
            }
        }
    }

    delete[] PS_TOTN_real;
    delete[] Sol_phi;
    delete[] Sol_phi0;
    delete[] B;
    delete[] epsilon;
    delete[] epsilon0;
    delete[] tmp_Vel;
    delete[] phi_tilda_R;
    delete[] phi_tilda_R0;

    ModuleBase::timer::tick("surchem", "cal_vel");
    return Vel;
}
