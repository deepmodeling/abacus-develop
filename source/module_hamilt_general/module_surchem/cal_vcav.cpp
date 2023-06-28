#include "module_base/timer.h"
#include "module_hamilt_general/module_xc/xc_functional.h"
#include "surchem.h"

void lapl_rho(const std::complex<double>* rhog, double* lapn, const ModulePW::PW_Basis* rho_basis)
{
    std::complex<double> *gdrtmpg = new std::complex<double>[rho_basis->npw];
    ModuleBase::GlobalFunc::ZEROS(lapn, rho_basis->nrxx);

    std::complex<double> *aux = new std::complex<double>[rho_basis->nmaxgr];

    // the formula is : rho(r)^prime = \int iG * rho(G)e^{iGr} dG
    for (int ig = 0; ig < rho_basis->npw; ig++)
        gdrtmpg[ig] = rhog[ig];
    for(int i = 0 ; i < 3 ; ++i)
    {
        // calculate the charge density gradient in reciprocal space.
        for (int ig = 0; ig < rho_basis->npw; ig++)
            aux[ig] = gdrtmpg[ig] * pow(rho_basis->gcar[ig][i], 2);
        // bring the gdr from G --> R
        rho_basis->recip2real(aux, aux);
        // remember to multily 2pi/a0, which belongs to G vectors.
        for (int ir = 0; ir < rho_basis->nrxx; ir++)
            lapn[ir] -= aux[ir].real() * GlobalC::ucell.tpiba2;
    }

    delete[] gdrtmpg;
    delete[] aux;
    return;
}

// calculates first derivative of the shape function in realspace
// exp(-(log(n/n_c))^2 /(2 sigma^2)) /(sigma * sqrt(2*pi) )/n
void shape_gradn(const complex<double>* PS_TOTN, const ModulePW::PW_Basis* rho_basis, double* eprime)
{

    double *PS_TOTN_real = new double[rho_basis->nrxx];
    ModuleBase::GlobalFunc::ZEROS(PS_TOTN_real, rho_basis->nrxx);
    rho_basis->recip2real(PS_TOTN, PS_TOTN_real);

    double epr_c = 1.0 / sqrt(ModuleBase::TWO_PI) / GlobalV::sigma_k;
    double epr_z = 0;
    double min = 1e-10;

    for (int ir = 0; ir < rho_basis->nrxx; ir++)
    {
        epr_z = log(std::max(PS_TOTN_real[ir], min) / GlobalV::nc_k) / sqrt(2) / GlobalV::sigma_k;
        eprime[ir] = epr_c * exp(-pow(epr_z, 2)) / std::max(PS_TOTN_real[ir], min);
    }

    delete[] PS_TOTN_real;
}

void surchem::createcavity(const UnitCell& ucell,
                           const ModulePW::PW_Basis* rho_basis,
                           const complex<double>* PS_TOTN,
                           double* vwork)
{
    ModuleBase::Vector3<double> *nablan = new ModuleBase::Vector3<double>[rho_basis->nrxx];
    ModuleBase::GlobalFunc::ZEROS(nablan, rho_basis->nrxx);
    double *nablan_2 = new double[rho_basis->nrxx];
    double *sqrt_nablan_2 = new double[rho_basis->nrxx];
    double *lapn = new double[rho_basis->nrxx];

    ModuleBase::GlobalFunc::ZEROS(nablan_2, rho_basis->nrxx);
    ModuleBase::GlobalFunc::ZEROS(sqrt_nablan_2, rho_basis->nrxx);
    ModuleBase::GlobalFunc::ZEROS(lapn, rho_basis->nrxx);

    // nabla n
    XC_Functional::grad_rho(PS_TOTN, nablan, rho_basis, GlobalC::ucell.tpiba);

    //  |\nabla n |^2 = nablan_2
    for (int ir = 0; ir < rho_basis->nrxx; ir++)
    {
        nablan_2[ir] = pow(nablan[ir].x, 2) + pow(nablan[ir].y, 2) + pow(nablan[ir].z, 2);
    }

    // Laplacian of n
    lapl_rho(PS_TOTN, lapn, rho_basis);

    //-------------------------------------------------------------
    // add -Lap(n)/|\nabla n| to vwork and copy \sqrt(|\nabla n|^2)
    // to sqrt_nablan_2
    //-------------------------------------------------------------

    double tmp = 0;
    double min = 1e-10;
    for (int ir = 0; ir < rho_basis->nrxx; ir++)
    {
        tmp = sqrt(std::max(nablan_2[ir], min));
        vwork[ir] = vwork[ir] - (lapn[ir]) / tmp;
        sqrt_nablan_2[ir] = tmp;
    }

    //-------------------------------------------------------------
    // term1 = gamma*A / n, where
    // gamma * A = exp(-(log(n/n_c))^2 /(2 sigma^2)) /(sigma * sqrt(2*pi) )
    //-------------------------------------------------------------
    double *term1 = new double[rho_basis->nrxx];
    shape_gradn(PS_TOTN, rho_basis, term1);

    //-------------------------------------------------------------
    // quantum surface area, integral of (gamma*A / n) * |\nabla n|
    //=term1 * sqrt_nablan_2
    //-------------------------------------------------------------
    qs = 0;

    for (int ir = 0; ir < rho_basis->nrxx; ir++)
    {
        qs = qs + (term1[ir]) * (sqrt_nablan_2[ir]);

        //   1 / |nabla n|
        sqrt_nablan_2[ir] = 1 / std::max(sqrt_nablan_2[ir], min);
    }

    //-------------------------------------------------------------
    // cavitation energy
    //-------------------------------------------------------------

    // double Ael = cal_Acav(ucell, pwb);

    //  packs the real array into a complex one
    //  to G space
    complex<double> *inv_gn = new complex<double>[rho_basis->npw];
    rho_basis->real2recip(sqrt_nablan_2, inv_gn);
    
    // \nabla(1 / |\nabla n|), ggn in real space
    ModuleBase::Vector3<double> *ggn = new ModuleBase::Vector3<double>[rho_basis->nrxx];
    XC_Functional::grad_rho(inv_gn, ggn, rho_basis, GlobalC::ucell.tpiba);

    //-------------------------------------------------------------
    // add -(\nabla n . \nabla(1/ |\nabla n|)) to Vcav in real space
    // and multiply by term1 = gamma*A/n in real space
    //-------------------------------------------------------------
    for (int ir = 0; ir < rho_basis->nrxx; ir++)
    {
        tmp = nablan[ir].x * ggn[ir].x + nablan[ir].y * ggn[ir].y + nablan[ir].z * ggn[ir].z;
        vwork[ir] = vwork[ir] - tmp;
    }

    for (int ir = 0; ir < rho_basis->nrxx; ir++)
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

ModuleBase::matrix surchem::cal_vcav(const UnitCell& ucell,
                                     const ModulePW::PW_Basis* rho_basis,
                                     complex<double>* PS_TOTN,
                                     int nspin)
{
    ModuleBase::TITLE("surchem", "cal_vcav");
    ModuleBase::timer::tick("surchem", "cal_vcav");

    double *tmp_Vcav = new double[rho_basis->nrxx];
    ModuleBase::GlobalFunc::ZEROS(tmp_Vcav, rho_basis->nrxx);

    createcavity(ucell, rho_basis, PS_TOTN, tmp_Vcav);

     ModuleBase::GlobalFunc::ZEROS(Vcav.c, nspin * rho_basis->nrxx);
    if (nspin == 4)
    {
        for (int ir = 0; ir < rho_basis->nrxx; ir++)
        {
            Vcav(0, ir) += tmp_Vcav[ir];
        }
    }
    else
    {
        for (int is = 0; is < nspin; is++)
        {
            for (int ir = 0; ir < rho_basis->nrxx; ir++)
            {
                Vcav(is, ir) += tmp_Vcav[ir];
            }
        }
    }

    delete[] tmp_Vcav;
    ModuleBase::timer::tick("surchem", "cal_vcav");
    return Vcav;
}