#include "write_elf.h"
#include "source_io/module_output/cube_io.h"
#ifdef _OPENMP
#include <omp.h>
#endif

namespace ModuleIO
{
void write_elf(
    const std::string& out_dir,
    const int& istep_in,
    const int& nspin,
    const double* const* rho,
    const double* const* tau,
    ModulePW::PW_Basis* const rho_basis,
    const Parallel_Grid& pgrid,
    const UnitCell* ucell_,
    const int& precision,
    const std::string& geom_block,
    const bool two_fermi)
{
    ModuleBase::timer::start("ModuleIO", "write_elf");
    // For nspin = 4, we only calculate the total ELF using the
    // rho_total and tau_total, containing in the first channel of
    // rho and tau.
    // What's more, we have not introduced the U(1) and SU(2) gauge
    // invariance corrections proposed by Desmarais J K, Vignale G,
    // Bencheikh K, et al. Physical Review Letters, 2024,
    // 133(13): 136401, where the current density is also included
    // in the ELF calculation.

    const int nspin_eff = (nspin == 4) ? 1 : nspin;

    const int nrxx = rho_basis->nrxx;
    const int npw = rho_basis->npw;

    assert(nrxx>0);
    assert(npw>0);

    std::vector<std::vector<double>> elf(nspin_eff, std::vector<double>(nrxx, 0.));
    // 1) calculate the kinetic energy density of vW KEDF
    std::vector<std::vector<double>> tau_vw(nspin_eff, std::vector<double>(nrxx, 0.));
    std::vector<double> phi(nrxx, 0.);

    for (int is = 0; is < nspin_eff; ++is)
    {
#pragma omp parallel for schedule(static) \
    default(none) firstprivate(nrxx) shared(phi, rho, rho_basis, is)
        for (int ir = 0; ir < nrxx; ++ir)
        {
            phi[ir] = std::sqrt(std::abs(rho[is][ir]));
        }

        std::vector<std::vector<double>> gradient_phi(3, std::vector<double>(nrxx, 0.));
        std::vector<std::complex<double>> recip_phi(npw, 0.0);
        std::vector<std::complex<double>> recip_gradient_phi(npw, 0.0);

        rho_basis->real2recip(phi.data(), recip_phi.data());

        std::complex<double> img(0.0, 1.0);
        for (int j = 0; j < 3; ++j)
        {
#pragma omp parallel for schedule(static) \
    default(none) firstprivate(img, npw) \
    shared(recip_gradient_phi, rho_basis, recip_phi, j)
            for (int ip = 0; ip < npw; ++ip)
            {
                recip_gradient_phi[ip]
                    = img * rho_basis->gcar[ip][j]
                    * recip_phi[ip] * rho_basis->tpiba;
            }

            rho_basis->recip2real(recip_gradient_phi.data(), gradient_phi[j].data());

#pragma omp parallel for schedule(static) \
    default(none) firstprivate(nrxx) \
    shared(tau_vw, gradient_phi, is, j, rho_basis)
            for (int ir = 0; ir < nrxx; ++ir)
            {
                tau_vw[is][ir] += gradient_phi[j][ir]
                    * gradient_phi[j][ir] / 2. * 2.;
            }
        }
    }

    // 2) calculate the kinetic energy density of TF KEDF
    std::vector<std::vector<double>> tau_TF(nspin_eff, std::vector<double>(nrxx, 0.));

    const double c_tf
        = 3.0 / 10.0 * std::pow(3 * std::pow(M_PI, 2.0), 2.0 / 3.0)
          * 2.0; // convert unit from Hartree to Ry
    if (nspin == 1 || nspin == 4)
    {
#pragma omp parallel for schedule(static) \
    default(none) firstprivate(c_tf, nrxx) \
    shared(rho, tau_TF, rho_basis)
        for (int ir = 0; ir < nrxx; ++ir)
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
        // spin-scaling: tau_TF[rho_up,rho_dn]
        //     = 1/2 * (tau_TF[2*rho_up] + tau_TF[2*rho_dn])
        for (int is = 0; is < nspin; ++is)
        {
#pragma omp parallel for schedule(static) \
    default(none) firstprivate(c_tf, nrxx) \
    shared(rho, tau_TF, is, rho_basis)
            for (int ir = 0; ir < nrxx; ++ir)
            {
                if (rho[is][ir] > 0.0)
                {
                    tau_TF[is][ir] = 0.5 * c_tf
                        * std::pow(2.0 * rho[is][ir], 5.0 / 3.0);
                }
                else
                {
                    tau_TF[is][ir] = 0.0;
                }
            }
        }
    }

    // 3) calculate the enhancement factor F = (tau_KS - tau_vw) /
    //    tau_TF, and then ELF = 1 / (1 + F^2)
    const double eps = 1.0e-5;
    for (int is = 0; is < nspin_eff; ++is)
    {
#pragma omp parallel for schedule(static) \
    default(none) firstprivate(eps, nrxx) \
    shared(elf, tau, tau_vw, tau_TF, is, rho_basis)
        for (int ir = 0; ir < nrxx; ++ir)
        {
            if (tau_TF[is][ir] > 1.0e-12)
            {
                elf[is][ir] = (tau[is][ir] - tau_vw[is][ir] + eps)
                    / tau_TF[is][ir];
                elf[is][ir] = 1. / (1. + elf[is][ir] * elf[is][ir]);
            }
            else
            {
                elf[is][ir] = 0.0;
            }
        }
    }

    // 4) output the ELF = 1 / (1 + F^2) to cube file
    const double ef_tmp = 0.0;
    const int out_fermi = 0;

    if (nspin == 1 || nspin == 4)
    {
        std::string fn = out_dir + "elftot" + geom_block + ".cube";

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
            out_fermi,
            two_fermi,
            false);
    }
    else if (nspin == 2)
    {
        for (int is = 0; is < nspin; ++is)
        {
            std::string fn_temp = out_dir + "elf" + "s"
                + std::to_string(is + 1) + geom_block + ".cube";

            const int ispin = is + 1;

            ModuleIO::write_vdata_palgrid(pgrid,
                elf[is].data(),
                ispin,
                nspin,
                istep_in,
                fn_temp,
                ef_tmp,
                ucell_,
                precision,
                out_fermi,
                two_fermi,
                false);
        }

        std::vector<double> elf_tot(nrxx, 0.0);
#pragma omp parallel for schedule(static) \
    default(none) firstprivate(eps, nrxx) \
    shared(elf_tot, tau, tau_vw, tau_TF, rho_basis)
        for (int ir = 0; ir < nrxx; ++ir)
        {
            if (tau_TF[0][ir] + tau_TF[1][ir] > 1.0e-12)
            {
                elf_tot[ir] = (tau[0][ir] + tau[1][ir]
                    - tau_vw[0][ir] - tau_vw[1][ir] + eps)
                    / (tau_TF[0][ir] + tau_TF[1][ir]);
                elf_tot[ir] = 1.
                    / (1. + elf_tot[ir] * elf_tot[ir]);
            }
            else
            {
                elf_tot[ir] = 0.0;
            }
        }
        std::string fn = out_dir + "elftot" + geom_block + ".cube";

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
            out_fermi,
            two_fermi,
            false);
    }
    ModuleBase::timer::end("ModuleIO", "write_elf");
} // end write_elf
} // end namespace ModuleIO
