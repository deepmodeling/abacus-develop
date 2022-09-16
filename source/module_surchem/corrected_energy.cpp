#include "surchem.h"

double surchem::cal_Ael(const UnitCell &cell, ModulePW::PW_Basis *rho_basis)
{
    double Ael = 0.0;
    for (int ir = 0; ir < rho_basis->nrxx; ir++)
    {
        Ael -= TOTN_real[ir] * delta_phi[ir];
    }
    Parallel_Reduce::reduce_double_pool(Ael);
    Ael = Ael * cell.omega / rho_basis->nxyz;
    // cout << "Ael: " << Ael << endl;
    return Ael;
}

double surchem::cal_Acav(const UnitCell &cell, ModulePW::PW_Basis *rho_basis)
{
    double Acav = 0.0;
    Acav = GlobalV::tau * qs;
    Acav = Acav * cell.omega / rho_basis->nxyz; // unit Ry
    Parallel_Reduce::reduce_double_pool(Acav);
    // cout << "Acav: " << Acav << endl;
    return Acav;
}

void surchem::cal_Acomp(const UnitCell &cell,
                        ModulePW::PW_Basis *rho_basis,
                        const double *const *const rho,
                        vector<double> &res)
{
    double Acomp1 = 0.0; // self
    double Acomp2 = 0.0; // electrons
    double Acomp3 = 0.0; // nuclear

    complex<double> *phi_comp_G = new complex<double>[rho_basis->npw];
    complex<double> *comp_reci = new complex<double>[rho_basis->npw];
    double *phi_comp_R = new double[rho_basis->nrxx];

    ModuleBase::GlobalFunc::ZEROS(phi_comp_G, rho_basis->npw);
    ModuleBase::GlobalFunc::ZEROS(comp_reci, rho_basis->npw);
    ModuleBase::GlobalFunc::ZEROS(phi_comp_R, rho_basis->nrxx);

    // part1: comp & comp
    rho_basis->real2recip(comp_real, comp_reci);
    for (int ig = 0; ig < rho_basis->npw; ig++)
    {
        if (rho_basis->gg[ig] >= 1.0e-12) // LiuXh 20180410
        {
            const double fac = ModuleBase::e2 * ModuleBase::FOUR_PI / (cell.tpiba2 * rho_basis->gg[ig]);
            Acomp1 += (conj(comp_reci[ig]) * comp_reci[ig]).real() * fac;
            phi_comp_G[ig] = fac * comp_reci[ig];
        }
    }
    // 0.5 for double counting
    Parallel_Reduce::reduce_double_pool(Acomp1);
    Acomp1 *= 0.5 * cell.omega;

    // electrons
    double *n_elec_R = new double[rho_basis->nrxx];
    for (int i = 0; i < rho_basis->nrxx; i++)
        n_elec_R[i] = 0.0;
    const int nspin0 = (GlobalV::NSPIN == 2) ? 2 : 1;
    for (int is = 0; is < nspin0; is++)
        for (int ir = 0; ir < rho_basis->nrxx; ir++)
            n_elec_R[ir] += rho[is][ir];

    // nuclear = TOTN_R + n_elec_R
    double *n_nucl_R = new double[rho_basis->nrxx];
    for (int ir = 0; ir < rho_basis->nrxx; ir++)
    {
        n_nucl_R[ir] = TOTN_real[ir] + n_elec_R[ir];
    }

    // part2: electrons
    rho_basis->recip2real(phi_comp_G, phi_comp_R);
    for (int ir = 0; ir < rho_basis->nrxx; ir++)
    {
        Acomp2 += n_elec_R[ir] * phi_comp_R[ir];
    }
    Parallel_Reduce::reduce_double_pool(Acomp2);
    Acomp2 = Acomp2 * cell.omega / rho_basis->nxyz;

    // part3: nuclear
    for (int ir = 0; ir < rho_basis->nrxx; ir++)
    {
        Acomp3 += n_nucl_R[ir] * phi_comp_R[ir];
    }
    Parallel_Reduce::reduce_double_pool(Acomp3);
    Acomp3 = Acomp3 * cell.omega / rho_basis->nxyz;

    delete[] phi_comp_G;
    delete[] phi_comp_R;
    delete[] comp_reci;

    delete[] n_elec_R;
    delete[] n_nucl_R;

    // cout << "Acomp1(self, Ry): " << Acomp1 << endl;
    // cout << "Acomp1(electrons, Ry): " << Acomp2 << endl;
    // cout << "Acomp1(nuclear, Ry): " << Acomp3 << endl;

    res[0] = Acomp1;
    res[1] = Acomp2;
    res[2] = -Acomp3;

    // return Acomp1 + Acomp2 - Acomp3;
}