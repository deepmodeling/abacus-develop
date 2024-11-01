#ifndef GET_PCHG_PW_H
#define GET_PCHG_PW_H

#include "cube_io.h"
#include "module_base/module_device/device.h"
#include "module_base/tool_quit.h"
#include "module_basis/module_pw/pw_basis.h"
#include "module_basis/module_pw/pw_basis_k.h"
#include "module_cell/unitcell.h"
#include "module_elecstate/elecstate.h"
#include "module_elecstate/module_charge/symmetry_rho.h"
#include "module_hamilt_pw/hamilt_pwdft/parallel_grid.h"
#include "module_psi/psi.h"

#include <string>
#include <vector>

namespace ModuleIO
{
template <typename Device>
void get_pchg_pw(const std::vector<int>& bands_to_print,
                 const int nbands,
                 const int nspin,
                 const int nx,
                 const int ny,
                 const int nz,
                 const int nxyz,
                 const int nks,
                 const std::vector<int>& isk,
                 const std::vector<double>& wk,
                 const int pw_big_bz,
                 const int pw_big_nbz,
                 const int ngmc,
                 UnitCell* ucell,
                 const psi::Psi<std::complex<double>>* psi,
                 const ModulePW::PW_Basis* pw_rhod,
                 const ModulePW::PW_Basis_K* pw_wfc,
                 const Device* ctx,
                 Parallel_Grid& Pgrid,
                 const std::string& global_out_dir,
                 const bool if_separate_k)
{
    // bands_picked is a vector of 0s and 1s, where 1 means the band is picked to output
    std::vector<int> bands_picked(nbands, 0);

    // Check if length of bands_to_print is valid
    if (static_cast<int>(bands_to_print.size()) > nbands)
    {
        ModuleBase::WARNING_QUIT("ESolver_KS_PW::get_pchg_pw",
                                 "The number of bands specified by `bands_to_print` in the "
                                 "INPUT file exceeds `nbands`!");
    }

    // Check if all elements in bands_picked are 0 or 1
    for (int value: bands_to_print)
    {
        if (value != 0 && value != 1)
        {
            ModuleBase::WARNING_QUIT("ESolver_KS_PW::get_pchg_pw",
                                     "The elements of `bands_to_print` must be either 0 or 1. "
                                     "Invalid values found!");
        }
    }

    // Fill bands_picked with values from bands_to_print
    // Remaining bands are already set to 0
    int length = std::min(static_cast<int>(bands_to_print.size()), nbands);
    for (int i = 0; i < length; ++i)
    {
        // bands_to_print rely on function parse_expression
        bands_picked[i] = static_cast<int>(bands_to_print[i]);
    }

    std::vector<std::complex<double>> wfcr(nxyz);
    std::vector<std::vector<double>> rho_band(nspin, std::vector<double>(nxyz));

    for (int ib = 0; ib < nbands; ++ib)
    {
        // Skip the loop iteration if bands_picked[ib] is 0
        if (!bands_picked[ib])
        {
            continue;
        }

        for (int is = 0; is < nspin; ++is)
        {
            std::fill(rho_band[is].begin(), rho_band[is].end(), 0.0);
        }

        if (if_separate_k)
        {
            for (int ik = 0; ik < nks; ++ik)
            {
                const int spin_index = isk[ik];
                std::cout << " Calculating band-decomposed charge density for band " << ib + 1 << ", k-point "
                          << ik % (nks / nspin) + 1 << ", spin " << spin_index + 1 << std::endl;

                psi->fix_k(ik);
                pw_wfc->recip_to_real(ctx, &psi[0](ib, 0), wfcr.data(), ik);

                // To ensure the normalization of charge density in multi-k calculation (if if_separate_k is true)
                double wg_sum_k = 0;
                for (int ik_tmp = 0; ik_tmp < nks / nspin; ++ik_tmp)
                {
                    wg_sum_k += wk[ik_tmp];
                }

                double w1 = static_cast<double>(wg_sum_k / ucell->omega);

                for (int i = 0; i < nxyz; ++i)
                {
                    rho_band[spin_index][i] = std::norm(wfcr[i]) * w1;
                }

                std::cout << " Writing cube files...";

                std::stringstream ssc;
                ssc << global_out_dir << "BAND" << ib + 1 << "_K" << ik % (nks / nspin) + 1 << "_SPIN" << spin_index + 1
                    << "_CHG.cube";

                ModuleIO::write_vdata_palgrid(GlobalC::Pgrid,
                    rho_band[spin_index].data(),
                    spin_index,
                    nspin,
                    0,
                    ssc.str(),
                    0.0,
                    ucell);

                std::cout << " Complete!" << std::endl;
            }
        }
        else
        {
            for (int ik = 0; ik < nks; ++ik)
            {
                const int spin_index = isk[ik];
                std::cout << " Calculating band-decomposed charge density for band " << ib + 1 << ", k-point "
                          << ik % (nks / nspin) + 1 << ", spin " << spin_index + 1 << std::endl;

                psi->fix_k(ik);
                pw_wfc->recip_to_real(ctx, &psi[0](ib, 0), wfcr.data(), ik);

                double w1 = static_cast<double>(wk[ik] / ucell->omega);

                for (int i = 0; i < nxyz; ++i)
                {
                    rho_band[spin_index][i] += std::norm(wfcr[i]) * w1;
                }
            }

            // Symmetrize the charge density, otherwise the results are incorrect if the symmetry is on
            std::cout << " Symmetrizing band-decomposed charge density..." << std::endl;
            Symmetry_rho srho;
            for (int is = 0; is < nspin; ++is)
            {
                // Use vector instead of raw pointers
                std::vector<double*> rho_save_pointers(nspin);
                for (int s = 0; s < nspin; ++s)
                {
                    rho_save_pointers[s] = rho_band[s].data();
                }

                std::vector<std::vector<std::complex<double>>> rhog(nspin, std::vector<std::complex<double>>(ngmc));

                // Convert vector of vectors to vector of pointers
                std::vector<std::complex<double>*> rhog_pointers(nspin);
                for (int s = 0; s < nspin; ++s)
                {
                    rhog_pointers[s] = rhog[s].data();
                }

                srho.begin(is,
                           rho_save_pointers.data(),
                           rhog_pointers.data(),
                           ngmc,
                           nullptr,
                           pw_rhod,
                           ucell->symm);
            }

            std::cout << " Writing cube files...";

            for (int is = 0; is < nspin; ++is)
            {
                std::stringstream ssc;
                ssc << global_out_dir << "BAND" << ib + 1 << "_SPIN" << is + 1 << "_CHG.cube";

                ModuleIO::write_vdata_palgrid(GlobalC::Pgrid,
                    rho_band[is].data(),
                    is,
                    nspin,
                    0,
                    ssc.str(),
                    0.0,
                    ucell);
            }

            std::cout << " Complete!" << std::endl;
        }
    }
}
} // namespace ModuleIO

#endif // GET_PCHG_PW_H
