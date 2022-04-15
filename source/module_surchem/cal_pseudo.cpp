#include "surchem.h"

atom_in surchem::GetAtom;

void surchem::gauss_charge(const UnitCell &cell, PW_Basis &pwb, complex<double> *N)
{
    pwb.setup_structure_factor(); // call strucFac(ntype,ngmc)
    for (int it = 0; it < cell.ntype; it++)
    {
        double RCS = GetAtom.atom_RCS[cell.atoms[it].psd];
        // double RCS = GetAtom.get_RCS(cell.atoms[it].psd);
        double sigma_rc_k = RCS / 2.5;
        double sigma_nc_k = RCS / 10.0;
        for (int ig = 0; ig < pwb.ngmc; ig++)
        {
            // G^2
            double gg = pwb.get_NormG_cartesian(ig);
            gg = gg * cell.tpiba2;

            N[ig].real(N[ig].real()
                       + (GetAtom.atom_Z[cell.atoms[it].psd] - cell.atoms[it].zv) * pwb.strucFac(it, ig).real()
                             * exp(-0.5 * gg * (sigma_rc_k * sigma_rc_k)));
            N[ig].imag(N[ig].imag()
                       + (GetAtom.atom_Z[cell.atoms[it].psd] - cell.atoms[it].zv) * pwb.strucFac(it, ig).imag()
                             * exp(-0.5 * gg * (sigma_rc_k * sigma_rc_k)));
        }
    }
    for (int ig = 0; ig < pwb.ngmc; ig++)
    {
        N[ig] /= cell.omega;
    }
}

void surchem::cal_pseudo(const UnitCell &cell, PW_Basis &pwb, const complex<double> *Porter_g, complex<double> *PS_TOTN)
{
    complex<double> *N = new complex<double>[pwb.ngmc];
    ModuleBase::GlobalFunc::ZEROS(N, pwb.ngmc);
    ModuleBase::GlobalFunc::ZEROS(PS_TOTN, pwb.ngmc);

    gauss_charge(cell, pwb, N);

    for (int ig = 0; ig < pwb.ngmc; ig++)
    {
        PS_TOTN[ig] = N[ig] + Porter_g[ig];
    }

    delete[] N;
}
