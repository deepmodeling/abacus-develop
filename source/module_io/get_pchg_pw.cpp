#include "get_pchg_pw.h"

#include "module_base/global_variable.h"
#include "module_base/tool_quit.h"
#include "rho_io.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>

namespace ModuleIO
{
void calculate_band_decomposed_charge_density(const std::vector<int>& bands_to_print,
                                              int nbands,
                                              int nspin,
                                              int pw_rho_nxyz,
                                              int kv_nks,
                                              int kv_isk[],
                                              double kv_wk[],
                                              double ucell_omega,
                                              int pw_rho_nx,
                                              int pw_rho_ny,
                                              int pw_rho_nz,
                                              const std::string& global_out_dir,
                                              bool if_separate_k,
                                              std::vector<std::complex<double>>& wfcr,
                                              std::vector<std::vector<double>>& rho_band,
                                              ModulePW::PW_Basis* pw_rhod
                                              // Add other necessary parameters here
)
{
    // bands_picked is a vector of 0s and 1s, where 1 means the band is picked to output
    std::vector<int> bands_picked(nbands, 0);

    // Check if length of bands_to_print is valid
    if (static_cast<int>(bands_to_print.size()) > nbands)
    {
        ModuleBase::WARNING_QUIT("calculate_band_decomposed_charge_density",
                                 "The number of bands specified by `bands_to_print` exceeds `nbands`!");
    }

    // Check if all elements in bands_picked are 0 or 1
    for (int value: bands_to_print)
    {
        if (value != 0 && value != 1)
        {
            ModuleBase::WARNING_QUIT("calculate_band_decomposed_charge_density",
                                     "The elements of `bands_to_print` must be either 0 or 1. Invalid values found!");
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
            for (int ik = 0; ik < kv_nks; ++ik)
            {
                const int spin_index = kv_isk[ik];
                std::cout << " Calculating band-decomposed charge density for band " << ib + 1 << ", k-point "
                          << ik % (kv_nks / nspin) + 1 << ", spin " << spin_index + 1 << std::endl;

                this->psi->fix_k(ik);
                this->pw_wfc->recip_to_real(this->ctx, &psi[0](ib, 0), wfcr.data(), ik);

                // To ensure the normalization of charge density in multi-k calculation (if if_separate_k is true)
                double wg_sum_k = 0;
                for (int ik_tmp = 0; ik_tmp < kv_nks / nspin; ++ik_tmp)
                {
                    wg_sum_k += kv_wk[ik_tmp];
                }

                double w1 = static_cast<double>(wg_sum_k / ucell_omega);

                for (int i = 0; i < pw_rho_nxyz; ++i)
                {
                    rho_band[spin_index][i] = std::norm(wfcr[i]) * w1;
                }

                std::cout << " Writing cube files...";

                std::stringstream ssc;
                ssc << global_out_dir << "BAND" << ib + 1 << "_K" << ik % (kv_nks / nspin) + 1 << "_SPIN"
                    << spin_index + 1 << "_CHG.cube";

                ModuleIO::write_rho(
#ifdef __MPI
                    this->pw_big->bz,
                    this->pw_big->nbz,
                    this->pw_big->nplane,
                    this->pw_big->startz_current,
#endif
                    rho_band[spin_index].data(),
                    spin_index,
                    nspin,
                    0,
                    ssc.str(),
                    pw_rho_nx,
                    pw_rho_ny,
                    pw_rho_nz,
                    0.0,
                    &(GlobalC::ucell),
                    11);

                std::cout << " Complete!" << std::endl;
            }
        }
        else
        {
            for (int ik = 0; ik < kv_nks; ++ik)
            {
                const int spin_index = kv_isk[ik];
                std::cout << " Calculating band-decomposed charge density for band " << ib + 1 << ", k-point "
                          << ik % (kv_nks / nspin) + 1 << ", spin " << spin_index + 1 << std::endl;

                this->psi->fix_k(ik);
                this->pw_wfc->recip_to_real(this->ctx, &psi[0](ib, 0), wfcr.data(), ik);

                double w1 = static_cast<double>(kv_wk[ik] / ucell_omega);

                for (int i = 0; i < pw_rho_nxyz; ++i)
                {
                    rho_band[spin_index][i] += std::norm(wfcr[i]) * w1;
                }
            }

            // Symmetrize the charge density, otherwise the results are incorrect if the symmetry is on
            std::cout << " Symmetrizing band-decomposed charge density..." << std::endl;
            Symmetry_rho srho;
            for (int is = 0; is < nspin; ++is)
            {
                std::vector<double*> rho_save_pointers(nspin);
                for (int s = 0; s < nspin; ++s)
                {
                    rho_save_pointers[s] = rho_band[s].data();
                }

                std::vector<std::vector<std::complex<double>>> rhog(
                    nspin,
                    std::vector<std::complex<double>>(this->pelec->charge->ngmc));

                std::vector<std::complex<double>*> rhog_pointers(nspin);
                for (int s = 0; s < nspin; ++s)
                {
                    rhog_pointers[s] = rhog[s].data();
                }

                srho.begin(is,
                           rho_save_pointers.data(),
                           rhog_pointers.data(),
                           this->pelec->charge->ngmc,
                           nullptr,
                           pw_rhod,
                           GlobalC::Pgrid,
                           GlobalC::ucell.symm);
            }

            std::cout << " Writing cube files...";

            for (int is = 0; is < nspin; ++is)
            {
                std::stringstream ssc;
                ssc << global_out_dir << "BAND" << ib + 1 << "_SPIN" << is + 1 << "_CHG.cube";

                ModuleIO::write_rho(
#ifdef __MPI
                    this->pw_big->bz,
                    this->pw_big->nbz,
                    this->pw_big->nplane,
                    this->pw_big->startz_current,
#endif
                    rho_band[is].data(),
                    is,
                    nspin,
                    0,
                    ssc.str(),
                    pw_rho_nx,
                    pw_rho_ny,
                    pw_rho_nz,
                    0.0,
                    &(GlobalC::ucell),
                    11);
            }

            std::cout << " Complete!" << std::endl;
        }
    }
}
} // namespace ModuleIO
