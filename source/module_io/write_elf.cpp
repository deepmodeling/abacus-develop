#include "write_elf.h"
#include "module_io/cube_io.h"

namespace ModuleIO
{
void write_elf(
#ifdef __MPI
    const int& bz,
    const int& nbz,
#endif
    const std::string& out_dir,
    const int& istep,
    const int& nspin,
    const double* const* rho,
    const double* const* tau,
    ModulePW::PW_Basis* rho_basis,
    const UnitCell* ucell_,
    const int& precision)
{
    std::vector<std::vector<double>> elf(nspin, std::vector<double>(rho_basis->nrxx, 0.));
    // 1) calculate the kinetic energy density of vW KEDF
    std::vector<std::vector<double>> tau_vw(nspin, std::vector<double>(rho_basis->nrxx, 0.));
    for (int is = 0; is < nspin; ++is)
    {
        std::vector<std::vector<double>> gradient_rho(3, std::vector<double>(rho_basis->nrxx, 0.));

        std::vector<std::complex<double>> recip_rho(rho_basis->npw, 0.0);
        std::vector<std::complex<double>> recip_gradient_rho(rho_basis->npw, 0.0);
        rho_basis->real2recip(rho[is], recip_rho.data());
        
        std::complex<double> img(0.0, 1.0);
        for (int j = 0; j < 3; ++j)
        {
            for (int ip = 0; ip < rho_basis->npw; ++ip)
            {
                recip_gradient_rho[ip] = img * rho_basis->gcar[ip][j] * recip_rho[ip] * rho_basis->tpiba;
            }

            rho_basis->recip2real(recip_gradient_rho.data(), gradient_rho[j].data());

            for (int ir = 0; ir < rho_basis->nrxx; ++ir)
            {
                tau_vw[is][ir] += gradient_rho[j][ir] * gradient_rho[j][ir] / (8. * rho[is][ir]) * 2.0; // convert Ha to Ry.
            }
        }
    }

    // 2) calculate the kinetic energy density of TF KEDF
    std::vector<std::vector<double>> tau_TF(nspin, std::vector<double>(rho_basis->nrxx, 0.));
    const double c_tf
        = 3.0 / 10.0 * std::pow(3 * std::pow(M_PI, 2.0), 2.0 / 3.0)
          * 2.0; // 10/3*(3*pi^2)^{2/3}, multiply by 2 to convert unit from Hartree to Ry, finally in Ry*Bohr^(-2)
    if (nspin == 1)
    {
        for (int ir = 0; ir < rho_basis->nrxx; ++ir)
        {
            tau_TF[0][ir] = c_tf * std::pow(rho[0][ir], 5.0 / 3.0);
        }
    }
    else if (nspin == 2)
    {
        for (int is = 0; is < nspin; ++is)
        {
            for (int ir = 0; ir < rho_basis->nrxx; ++ir)
            {
                tau_TF[is][ir] = 0.5 * c_tf * std::pow(2.0 * rho[is][ir], 5.0 / 3.0);
            }
        }
    }

    // 3) calculate the enhancement factor F = (tau_KS - tau_vw) / tau_TF, and then ELF = 1 / (1 + F^2)
    double eps = 1.0e-5; // suppress the numerical instability in LCAO (Ref: Acta Phys. -Chim. Sin. 2011, 27(12), 2786-2792. doi: 10.3866/PKU.WHXB20112786)
    for (int is = 0; is < nspin; ++is)
    {
        for (int ir = 0; ir < rho_basis->nrxx; ++ir)
        {
            elf[is][ir] = (tau[is][ir] - tau_vw[is][ir] + eps) / tau_TF[is][ir];
            elf[is][ir] = 1. / (1. + elf[is][ir] * elf[is][ir]);
        }
    }

    // 4) output the ELF = 1 / (1 + F^2) to cube file
    double ef_tmp = 0.0;
    int out_fermi = 0;

    if (nspin == 1)
    {
        std::string fn = out_dir + "/ELF.cube";

        int is = -1;
        ModuleIO::write_cube(
    #ifdef __MPI
            bz,
            nbz,
            rho_basis->nplane,
            rho_basis->startz_current,
    #endif
            elf[0].data(),
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
    else if (nspin == 2)
    {
        for (int is = 0; is < nspin; ++is)
        {
            std::string fn_temp = out_dir + "/ELF_SPIN" + std::to_string(is + 1) + ".cube";
            int ispin = is + 1;

            ModuleIO::write_cube(
        #ifdef __MPI
                bz,
                nbz,
                rho_basis->nplane,
                rho_basis->startz_current,
        #endif
                elf[is].data(),
                ispin,
                nspin,
                istep,
                fn_temp,
                rho_basis->nx,
                rho_basis->ny,
                rho_basis->nz,
                ef_tmp,
                ucell_,
                precision,
                out_fermi);   
        }

        std::vector<double> elf_tot(rho_basis->nrxx, 0.0);
        for (int ir = 0; ir < rho_basis->nrxx; ++ir)
        {
            elf_tot[ir] = (tau[0][ir] + tau[1][ir] - tau_vw[0][ir] - tau_vw[1][ir]) / (tau_TF[0][ir] + tau_TF[1][ir]);
            elf_tot[ir] = 1. / (1. + elf_tot[ir] * elf_tot[ir]);
        }
        std::string fn = out_dir + "/ELF.cube";

        int is = -1;
        ModuleIO::write_cube(
    #ifdef __MPI
            bz,
            nbz,
            rho_basis->nplane,
            rho_basis->startz_current,
    #endif
            elf_tot.data(),
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
}