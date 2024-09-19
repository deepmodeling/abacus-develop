#include "write_elf.h"
#include "module_io/cube_io.h"

namespace ModuleIO
{
void write_elf(
#ifdef __MPI
    const int& bz,
    const int& nbz,
#endif
    const std::string& fn,
    const int& istep,
    const int& nspin,
    const double* const* rho,
    const double* const* tau,
    ModulePW::PW_Basis* rho_basis,
    const UnitCell* ucell_)
{
    std::vector<double> result(rho_basis->nrxx, 0.);
    // 1) calculate the kinetic energy density of vW KEDF
    std::vector<double> tau_vw(rho_basis->nrxx, 0.);
    for (int is = 0; is < nspin; ++is)
    {
        std::vector<std::vector<double>> nabla_rho(3, std::vector<double>(rho_basis->nrxx, 0.));

        std::complex<double> *recip_rho = new std::complex<double>[rho_basis->npw];
        std::complex<double> *recip_nabla_rho = new std::complex<double>[rho_basis->npw];
        rho_basis->real2recip(rho[is], recip_rho);
        
        std::complex<double> img(0.0, 1.0);
        for (int j = 0; j < 3; ++j)
        {
            for (int ip = 0; ip < rho_basis->npw; ++ip)
            {
                recip_nabla_rho[ip] = img * rho_basis->gcar[ip][j] * recip_rho[ip] * rho_basis->tpiba;
            }

            rho_basis->recip2real(recip_nabla_rho, nabla_rho[j].data());

            for (int ir = 0; ir < rho_basis->nrxx; ++ir)
            {
                tau_vw[ir] += nabla_rho[j][ir] * nabla_rho[j][ir] / (8. * rho[is][ir]) * 2.0; // convert Ha to Ry.
            }
        }
    }

    // 2) calculate the kinetic energy density of TF KEDF
    std::vector<double> tau_TF(rho_basis->nrxx, 0.);
    const double c_tf
        = 3.0 / 10.0 * std::pow(3 * std::pow(M_PI, 2.0), 2.0 / 3.0)
          * 2.0; // 10/3*(3*pi^2)^{2/3}, multiply by 2 to convert unit from Hartree to Ry, finally in Ry*Bohr^(-2)
    if (nspin == 1)
    {
        for (int ir = 0; ir < rho_basis->nrxx; ++ir)
        {
            tau_TF[ir] = c_tf * std::pow(rho[0][ir], 5.0 / 3.0);
        }
    }
    else if (nspin == 2)
    {
        for (int is = 0; is < nspin; ++is)
        {
            for (int ir = 0; ir < rho_basis->nrxx; ++ir)
            {
                tau_TF[ir] += 0.5 * c_tf * std::pow(2.0 * rho[is][ir], 5.0 / 3.0);
            }
        }
    }

    // 3) calculate the enhancement factor F = (tau_KS - tau_vw) / tau_TF
    for (int ir = 0; ir < rho_basis->nrxx; ++ir)
    {
        for (int is = 0; is < nspin; ++is)
        {
            result[ir] += tau[is][ir];
        }
        result[ir] = (result[ir] - tau_vw[ir]) / tau_TF[ir];
    }

    // 4) calculate the ELF = 1 / (1 + F^2)
    for (int ir = 0; ir < rho_basis->nrxx; ++ir)
    {
        result[ir] = 1. / (1. + result[ir] * result[ir]);
    }

    int precision = 9;
    int is = -1;
    double ef_tmp = 0.0;
    int out_fermi = 0;

    ModuleIO::write_cube(
#ifdef __MPI
        bz,
        nbz,
        rho_basis->nplane,
        rho_basis->startz_current,
#endif
        result.data(),
        is,
        nspin,
        istep,
        fn,
        rho_basis->nx,
        rho_basis->ny,
        rho_basis->nz,
        ef_tmp,
        ucell_,
        precision,
        out_fermi);
}
}