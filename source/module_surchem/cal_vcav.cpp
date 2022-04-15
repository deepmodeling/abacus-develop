#include "../module_base/timer.h"
#include "../module_xc/xc_functional.h"
#include "surchem.h"

void lapl_rho(PW_Basis &pwb, const std::complex<double> *rhog, double *lapn)
{
    std::complex<double> *gdrtmpg = new std::complex<double>[pwb.ngmc];
    ModuleBase::GlobalFunc::ZEROS(gdrtmpg, pwb.ngmc);
    ModuleBase::GlobalFunc::ZEROS(lapn, pwb.nrxx);

    std::complex<double> *Porter = GlobalC::UFFT.porter;

    // the formula is : rho(r)^prime = \int iG * rho(G)e^{iGr} dG
    for (int ig = 0; ig < pwb.ngmc; ig++)
        gdrtmpg[ig] = rhog[ig];

    // calculate the charge density gradient in reciprocal space.
    ModuleBase::GlobalFunc::ZEROS(Porter, pwb.nrxx);
    for (int ig = 0; ig < pwb.ngmc; ig++)
        Porter[pwb.ig2fftc[ig]]
            = gdrtmpg[ig] * std::complex<double>(pow(pwb.get_G_cartesian_projection(ig, 0), 2), 0.0);
    // bring the gdr from G --> R
    pwb.FFT_chg.FFT3D(Porter, 1);
    // remember to multily 2pi/a0, which belongs to G vectors.
    for (int ir = 0; ir < pwb.nrxx; ir++)
        lapn[ir] = lapn[ir] - Porter[ir].real() * GlobalC::ucell.tpiba2;

    // calculate the charge density gradient in reciprocal space.
    ModuleBase::GlobalFunc::ZEROS(Porter, pwb.nrxx);
    for (int ig = 0; ig < pwb.ngmc; ig++)
        Porter[pwb.ig2fftc[ig]]
            = gdrtmpg[ig] * std::complex<double>(pow(pwb.get_G_cartesian_projection(ig, 1), 2), 0.0);
    // bring the gdr from G --> R
    pwb.FFT_chg.FFT3D(Porter, 1);
    // remember to multily 2pi/a0, which belongs to G vectors.
    for (int ir = 0; ir < pwb.nrxx; ir++)
        lapn[ir] = lapn[ir] - Porter[ir].real() * GlobalC::ucell.tpiba2;

    // calculate the charge density gradient in reciprocal space.
    ModuleBase::GlobalFunc::ZEROS(Porter, pwb.nrxx);
    for (int ig = 0; ig < pwb.ngmc; ig++)
        Porter[pwb.ig2fftc[ig]]
            = gdrtmpg[ig] * std::complex<double>(pow(pwb.get_G_cartesian_projection(ig, 2), 2), 0.0);
    // bring the gdr from G --> R
    pwb.FFT_chg.FFT3D(Porter, 1);
    // remember to multily 2pi/a0, which belongs to G vectors.
    for (int ir = 0; ir < pwb.nrxx; ir++)
        lapn[ir] = lapn[ir] - Porter[ir].real() * GlobalC::ucell.tpiba2;

    delete[] gdrtmpg;
    return;
}

// calculates first derivative of the shape function in realspace
// exp(-(log(n/n_c))^2 /(2 sigma^2)) /(sigma * sqrt(2*pi) )/n
void shape_gradn(const complex<double> *PS_TOTN, PW_Basis &pw, double *eprime)
{

    double *PS_TOTN_real = new double[pw.nrxx];
    ModuleBase::GlobalFunc::ZEROS(PS_TOTN_real, pw.nrxx);
    GlobalC::UFFT.ToRealSpace(PS_TOTN, PS_TOTN_real);

    double epr_c = 1.0 / sqrt(ModuleBase::TWO_PI) / GlobalV::sigma_k;
    double epr_z = 0;
    double min = 1e-10;

    for (int ir = 0; ir < pw.nrxx; ir++)
    {
        epr_z = log(max(PS_TOTN_real[ir], min) / GlobalV::nc_k) / sqrt(2) / GlobalV::sigma_k;
        eprime[ir] = epr_c * exp(-pow(epr_z, 2)) / max(PS_TOTN_real[ir], min);
    }

    delete[] PS_TOTN_real;
}

void createcavity(const UnitCell &ucell, PW_Basis &pwb, const complex<double> *PS_TOTN, double *vwork)
{
    ModuleBase::Vector3<double> *nablan = new ModuleBase::Vector3<double>[pwb.nrxx];
    ModuleBase::GlobalFunc::ZEROS(nablan, pwb.nrxx);
    double *nablan_2 = new double[pwb.nrxx];
    double *sqrt_nablan_2 = new double[pwb.nrxx];
    double *lapn = new double[pwb.nrxx];

    ModuleBase::GlobalFunc::ZEROS(nablan_2, pwb.nrxx);
    ModuleBase::GlobalFunc::ZEROS(sqrt_nablan_2, pwb.nrxx);
    ModuleBase::GlobalFunc::ZEROS(lapn, pwb.nrxx);

    // nabla n
    XC_Functional::grad_rho(PS_TOTN, nablan);

    //  |\nabla n |^2 = nablan_2
    for (int ir = 0; ir < pwb.nrxx; ir++)
    {
        nablan_2[ir] = pow(nablan[ir].x, 2) + pow(nablan[ir].y, 2) + pow(nablan[ir].z, 2);
    }

    // Laplacian of n
    lapl_rho(pwb, PS_TOTN, lapn);

    //-------------------------------------------------------------
    // add -Lap(n)/|\nabla n| to vwork and copy \sqrt(|\nabla n|^2)
    // to sqrt_nablan_2
    //-------------------------------------------------------------

    double tmp = 0;
    double min = 1e-10;
    for (int ir = 0; ir < pwb.nrxx; ir++)
    {
        tmp = sqrt(max(nablan_2[ir], min));
        vwork[ir] = vwork[ir] - (lapn[ir]) / tmp;
        sqrt_nablan_2[ir] = tmp;
    }

    //-------------------------------------------------------------
    // term1 = gamma*A / n, where
    // gamma * A = exp(-(log(n/n_c))^2 /(2 sigma^2)) /(sigma * sqrt(2*pi) )
    //-------------------------------------------------------------
    double *term1 = new double[pwb.nrxx];
    shape_gradn(PS_TOTN, pwb, term1);

    //-------------------------------------------------------------
    // quantum surface area, integral of (gamma*A / n) * |\nabla n|
    //=term1 * sqrt_nablan_2
    //-------------------------------------------------------------
    double qs = 0;

    for (int ir = 0; ir < pwb.nrxx; ir++)
    {
        qs = qs + (term1[ir]) * (sqrt_nablan_2[ir]);

        //   1 / |nabla n|
        sqrt_nablan_2[ir] = 1 / max(sqrt_nablan_2[ir], min);
    }

    //-------------------------------------------------------------
    // cavitation energy
    //-------------------------------------------------------------

    double Ael = surchem::cal_Acav(ucell, pwb, qs);

    //  packs the real array into a complex one
    //  to G space
    complex<double> *inv_gn = new complex<double>[pwb.ngmc];
    GlobalC::UFFT.ToReciSpace(sqrt_nablan_2, inv_gn);

    // \nabla(1 / |\nabla n|), ggn in real space
    ModuleBase::Vector3<double> *ggn = new ModuleBase::Vector3<double>[pwb.nrxx];
    XC_Functional::grad_rho(inv_gn, ggn);

    //-------------------------------------------------------------
    // add -(\nabla n . \nabla(1/ |\nabla n|)) to Vcav in real space
    // and multiply by term1 = gamma*A/n in real space
    //-------------------------------------------------------------
    for (int ir = 0; ir < pwb.nrxx; ir++)
    {
        tmp = nablan[ir].x * ggn[ir].x + nablan[ir].y * ggn[ir].y + nablan[ir].z * ggn[ir].z;
        vwork[ir] = vwork[ir] - tmp;
    }

    for (int ir = 0; ir < pwb.nrxx; ir++)
    {
        vwork[ir] = vwork[ir] * term1[ir] * GlobalV::tau;
    }

    delete[] nablan;
    delete[] nablan_2;
    delete[] sqrt_nablan_2;
    delete[] lapn;
    delete[] term1;
    delete[] inv_gn;
    delete[] ggn;
}

ModuleBase::matrix surchem::cal_vcav(const UnitCell &ucell, PW_Basis &pwb, const complex<double> *PS_TOTN, int nspin)
{
    ModuleBase::TITLE("surchem", "cal_vcav");
    ModuleBase::timer::tick("surchem", "cal_vcav");

    double *Vcav = new double[pwb.nrxx];
    ModuleBase::GlobalFunc::ZEROS(Vcav, pwb.nrxx);

    createcavity(ucell, pwb, PS_TOTN, Vcav);

    ModuleBase::matrix v(nspin, pwb.nrxx);
    if (nspin == 4)
    {
        for (int ir = 0; ir < pwb.nrxx; ir++)
        {
            v(0, ir) += Vcav[ir];
        }
    }
    else
    {
        for (int is = 0; is < nspin; is++)
        {
            for (int ir = 0; ir < pwb.nrxx; ir++)
            {
                v(is, ir) += Vcav[ir];
            }
        }
    }
    delete[] Vcav;
    ModuleBase::timer::tick("surchem", "cal_vcav");
    return v;
}