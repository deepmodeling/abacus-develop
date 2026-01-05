#include "source_base/timer.h"
#include "source_hamilt/module_xc/xc_functional.h"
#include "source_io/module_parameter/parameter.h"
#include "surchem.h"


static void local_grad_rho(const UnitCell& ucell,
                           const ModulePW::PW_Basis* rho_basis,
                           const std::complex<double>* rho_G, 
                           ModuleBase::Vector3<double>* grad_R, 
                           std::complex<double>* aux_G,        
                           double* aux_R)                     
{

    for (int i = 0; i < 3; ++i)
    {
        // 1.  iG * rho(G)
        for (int ig = 0; ig < rho_basis->npw; ig++) {
            aux_G[ig] = ModuleBase::IMAG_UNIT * rho_G[ig] * rho_basis->gcar[ig][i];
        }

        // 2. FFT: G -> R

        rho_basis->recip2real(aux_G, aux_R);

        // 3.  2pi/a
        for (int ir = 0; ir < rho_basis->nrxx; ir++) {
            grad_R[ir][i] = aux_R[ir] * ucell.tpiba;
        }
    }
}


void lapl_rho(const double& tpiba2,
              const std::complex<double>* rhog, 
              double* lapn, 
              const ModulePW::PW_Basis* rho_basis)
{
    std::complex<double> *gdrtmpg = new std::complex<double>[rho_basis->npw];
    ModuleBase::GlobalFunc::ZEROS(lapn, rho_basis->nrxx);


    std::complex<double> *aux = new std::complex<double>[rho_basis->nmaxgr];

 
    for (int ig = 0; ig < rho_basis->npw; ig++) {
        gdrtmpg[ig] = rhog[ig];
    }
    
    for(int i = 0 ; i < 3 ; ++i)
    {
 
        ModuleBase::GlobalFunc::ZEROS(aux, rho_basis->nmaxgr);
        
        for (int ig = 0; ig < rho_basis->npw; ig++) {
            aux[ig] = gdrtmpg[ig] * pow(rho_basis->gcar[ig][i], 2);
        }
        
   
        rho_basis->recip2real(aux, aux);
        

        for (int ir = 0; ir < rho_basis->nrxx; ir++) {
            lapn[ir] -= aux[ir].real() * tpiba2;
        }
    }

    delete[] gdrtmpg;
    delete[] aux;
}

// ---------------------------------------------------------
// Helper 3: Shape function 
// ---------------------------------------------------------
void shape_gradn(const std::complex<double>* ps_totn, const ModulePW::PW_Basis* rho_basis, double* eprime)
{
    double *ps_totn_real = new double[rho_basis->nrxx];
    ModuleBase::GlobalFunc::ZEROS(ps_totn_real, rho_basis->nrxx);
    rho_basis->recip2real(ps_totn, ps_totn_real);

    double epr_c = 1.0 / sqrt(ModuleBase::TWO_PI) / PARAM.inp.sigma_k;
    double epr_z = 0;
    double min = 1e-10;

    for (int ir = 0; ir < rho_basis->nrxx; ir++)
    {
        epr_z = log(std::max(ps_totn_real[ir], min) / PARAM.inp.nc_k) / sqrt(2) / PARAM.inp.sigma_k;
        eprime[ir] = epr_c * exp(-pow(epr_z, 2)) / std::max(ps_totn_real[ir], min);
    }

    delete[] ps_totn_real;
}

// ---------------------------------------------------------
// Create Cavity
// ---------------------------------------------------------
void surchem::createcavity(const UnitCell& ucell,
                           const ModulePW::PW_Basis* rho_basis,
                           const std::complex<double>* ps_totn,
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


    std::complex<double> *aux_G = new std::complex<double>[rho_basis->npw];
    double *aux_R = new double[rho_basis->nrxx];

    local_grad_rho(ucell, rho_basis, ps_totn, nablan, aux_G, aux_R);

    // | \nabla n |^2
    for (int ir = 0; ir < rho_basis->nrxx; ir++)
    {
        nablan_2[ir] = pow(nablan[ir].x, 2) + pow(nablan[ir].y, 2) + pow(nablan[ir].z, 2);
    }


    lapl_rho(ucell.tpiba2, ps_totn, lapn, rho_basis);


    double tmp = 0;
    double min = 1e-10;
    for (int ir = 0; ir < rho_basis->nrxx; ir++)
    {
        tmp = sqrt(std::max(nablan_2[ir], min));
        vwork[ir] = vwork[ir] - (lapn[ir]) / tmp;
        sqrt_nablan_2[ir] = tmp;
    }

    // --------------------------------------------------------
    // term1 = gamma*A / n
    // --------------------------------------------------------
    double *term1 = new double[rho_basis->nrxx];
    shape_gradn(ps_totn, rho_basis, term1);

    // --------------------------------------------------------
    // quantum surface area calculation
    // --------------------------------------------------------
    qs = 0;
    for (int ir = 0; ir < rho_basis->nrxx; ir++)
    {
        qs = qs + (term1[ir]) * (sqrt_nablan_2[ir]);
        // 1 / |nabla n|
        sqrt_nablan_2[ir] = 1 / std::max(sqrt_nablan_2[ir], min);
    }

    // Cavitation energy
    this->Acav = PARAM.inp.tau * qs * ucell.omega / rho_basis->nxyz;
    Parallel_Reduce::reduce_pool(this->Acav);


    std::complex<double> *inv_gn = new std::complex<double>[rho_basis->npw];
    rho_basis->real2recip(sqrt_nablan_2, inv_gn);
    

    ModuleBase::Vector3<double> *ggn = new ModuleBase::Vector3<double>[rho_basis->nrxx];
    
 
    local_grad_rho(ucell, rho_basis, inv_gn, ggn, aux_G, aux_R);

    // --------------------------------------------------------
    // add -(\nabla n . \nabla(1/ |\nabla n|)) to Vcav
    // --------------------------------------------------------
    for (int ir = 0; ir < rho_basis->nrxx; ir++)
    {
        tmp = nablan[ir].x * ggn[ir].x + nablan[ir].y * ggn[ir].y + nablan[ir].z * ggn[ir].z;
        vwork[ir] = vwork[ir] - tmp;
    }

    for (int ir = 0; ir < rho_basis->nrxx; ir++)
    {
        vwork[ir] = vwork[ir] * term1[ir] * PARAM.inp.tau;
    }


    delete[] aux_G;     
    delete[] aux_R;    
    delete[] nablan;
    delete[] nablan_2;
    delete[] sqrt_nablan_2;
    delete[] lapn;
    delete[] term1;
    delete[] inv_gn;
    delete[] ggn;
}

// ---------------------------------------------------------
// Main Function: cal_vcav
// ---------------------------------------------------------
void surchem::cal_vcav(const UnitCell& ucell,
                       const ModulePW::PW_Basis* rho_basis,
                       std::complex<double>* ps_totn,
                       int nspin,
                       ModuleBase::matrix& v)
{
    ModuleBase::TITLE("surchem", "cal_vcav");
    ModuleBase::timer::tick("surchem", "cal_vcav");

    double *tmp_Vcav = new double[rho_basis->nrxx];
    ModuleBase::GlobalFunc::ZEROS(tmp_Vcav, rho_basis->nrxx);

    createcavity(ucell, rho_basis, ps_totn, tmp_Vcav);

    if (nspin == 4)
    {
        for (int ir = 0; ir < rho_basis->nrxx; ir++)
        {
            Vcav(0, ir) = tmp_Vcav[ir];
            v(0, ir) += Vcav(0, ir);
        }
    }
    else
    {
        for (int is = 0; is < nspin; is++)
        {
            for (int ir = 0; ir < rho_basis->nrxx; ir++)
            {
                Vcav(is, ir) = tmp_Vcav[ir];
                v(is, ir) += Vcav(is, ir);
            }
        }
    }

    delete[] tmp_Vcav;
    ModuleBase::timer::tick("surchem", "cal_vcav");
    return;
}