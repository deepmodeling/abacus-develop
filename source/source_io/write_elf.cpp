#include "write_elf.h"
#include "source_io/cube_io.h"

namespace ModuleIO
{
void write_elf(
#ifdef __MPI
    const int& bz,
    const int& nbz,
#endif
    const std::string& out_dir,
    const int& istep_in,
    const int& nspin,
    const double* const* rho,
    const double* const* tau,
    ModulePW::PW_Basis* rho_basis,
    const Parallel_Grid& pgrid,
    const UnitCell* ucell_,
    const int& precision)
{
    // For nspin = 4, we only calculate the total ELF using the rho_total and tau_total,
    // containing in the first channel of rho and tau.
    // What's more, we have not introduced the U(1) and SU(2) gauge invariance corrections
    // proposed by Desmarais J K, Vignale G, Bencheikh K, et al. Physical Review Letters, 2024, 133(13): 136401,
    // where the current density is also included in the ELF calculation.

    int nspin_eff = (nspin == 4) ? 1 : nspin;

    std::vector<std::vector<double>> elf(nspin_eff, std::vector<double>(rho_basis->nrxx, 0.));
    // 1) calculate the kinetic energy density of vW KEDF
    std::vector<std::vector<double>> tau_vw(nspin_eff, std::vector<double>(rho_basis->nrxx, 0.));
    std::vector<double> phi(rho_basis->nrxx, 0.); // phi = sqrt(rho)

    for (int is = 0; is < nspin_eff; ++is)
    {
        for (int ir = 0; ir < rho_basis->nrxx; ++ir)
        {
            phi[ir] = std::sqrt(std::abs(rho[is][ir]));
        }
        
        std::vector<std::vector<double>> gradient_phi(3, std::vector<double>(rho_basis->nrxx, 0.));
        std::vector<std::complex<double>> recip_phi(rho_basis->npw, 0.0);
        std::vector<std::complex<double>> recip_gradient_phi(rho_basis->npw, 0.0);
        
        rho_basis->real2recip(phi.data(), recip_phi.data());
        
        std::complex<double> img(0.0, 1.0);
        for (int j = 0; j < 3; ++j)
        {
            for (int ip = 0; ip < rho_basis->npw; ++ip)
            {
                recip_gradient_phi[ip] = img * rho_basis->gcar[ip][j] * recip_phi[ip] * rho_basis->tpiba;
            }

            rho_basis->recip2real(recip_gradient_phi.data(), gradient_phi[j].data());

            for (int ir = 0; ir < rho_basis->nrxx; ++ir)
            {
                tau_vw[is][ir] += gradient_phi[j][ir] * gradient_phi[j][ir] / 2. * 2.; // convert Ha to Ry.
            }
        }
    }

    // 2) calculate the kinetic energy density of TF KEDF
    std::vector<std::vector<double>> tau_TF(nspin_eff, std::vector<double>(rho_basis->nrxx, 0.));
    const double c_tf
        = 3.0 / 10.0 * std::pow(3 * std::pow(M_PI, 2.0), 2.0 / 3.0)
          * 2.0; // 10/3*(3*pi^2)^{2/3}, multiply by 2 to convert unit from Hartree to Ry, finally in Ry*Bohr^(-2)
    if (nspin == 1 || nspin == 4)
    {
        for (int ir = 0; ir < rho_basis->nrxx; ++ir)
        {
            if (rho[0][ir] > 0.0)
            {
                tau_TF[0][ir] = c_tf * std::pow(rho[0][ir], 5.0 / 3.0);
            }
            else
            {
                tau_TF[0][ir] = 0.0;
            }
        }
    }
    else if (nspin == 2)
    {
        // the spin-scaling law: tau_TF[rho_up, rho_dn] = 1/2 * (tau_TF[2*rho_up] + tau_TF[2*rho_dn])
        for (int is = 0; is < nspin; ++is)
        {
            for (int ir = 0; ir < rho_basis->nrxx; ++ir)
            {
                if (rho[is][ir] > 0.0)
                {
                    tau_TF[is][ir] = 0.5 * c_tf * std::pow(2.0 * rho[is][ir], 5.0 / 3.0);
                }
                else
                {
                    tau_TF[is][ir] = 0.0;
                }
            }
        }
    }

    // 3) calculate the enhancement factor F = (tau_KS - tau_vw) / tau_TF, and then ELF = 1 / (1 + F^2)
    double eps = 1.0e-5; // suppress the numerical instability in LCAO (Ref: Acta Phys. -Chim. Sin. 2011, 27(12), 2786-2792. doi: 10.3866/PKU.WHXB20112786)
    for (int is = 0; is < nspin_eff; ++is)
    {
        for (int ir = 0; ir < rho_basis->nrxx; ++ir)
        {
            if (tau_TF[is][ir] > 1.0e-12)
            {
                elf[is][ir] = (tau[is][ir] - tau_vw[is][ir] + eps) / tau_TF[is][ir];
                elf[is][ir] = 1. / (1. + elf[is][ir] * elf[is][ir]);
            }
            else
            {
                elf[is][ir] = 0.0;
            }
        }
    }

    // 4) output the ELF = 1 / (1 + F^2) to cube file
    double ef_tmp = 0.0;
    int out_fermi = 0;

    if (nspin == 1 || nspin == 4)
    {
        std::string fn = out_dir + "/elf.cube";

        int is = -1;
        ModuleIO::write_vdata_palgrid(pgrid,
            elf[0].data(),
            is,
            nspin,
            istep_in,
            fn,
            ef_tmp,
            ucell_,
            precision,
            out_fermi);   
    }
    else if (nspin == 2)
    {
        for (int is = 0; is < nspin; ++is)
        {
            std::string fn_temp = out_dir + "/elf";

            fn_temp += std::to_string(is + 1) + ".cube";

            int ispin = is + 1;

            ModuleIO::write_vdata_palgrid(pgrid,
                elf[is].data(),
                ispin,
                nspin,
                istep_in,
                fn_temp,
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
        std::string fn = out_dir + "/elf.cube";

        int is = -1;
        ModuleIO::write_vdata_palgrid(pgrid,
            elf_tot.data(),
            is,
            nspin,
            istep_in,
            fn,
            ef_tmp,
            ucell_,
            precision,
            out_fermi);
    }
}
}
