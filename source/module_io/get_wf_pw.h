#ifndef GET_WF_PW_H
#define GET_WF_PW_H

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
void get_wf_pw(const std::vector<int>& out_wfc_norm,
               const std::vector<int>& out_wfc_re_im,
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
               const Parallel_Grid& pgrid,
               const std::string& global_out_dir)
{
    // bands_picked is a vector of 0s and 1s, where 1 means the band is picked to output
    std::vector<int> bands_picked_norm(nbands, 0);
    std::vector<int> bands_picked_re_im(nbands, 0);

    // Check if length of out_wfc_norm and out_wfc_re_im is valid
    if (static_cast<int>(out_wfc_norm.size()) > nbands || static_cast<int>(out_wfc_re_im.size()) > nbands)
    {
        ModuleBase::WARNING_QUIT("ModuleIO::get_wf_pw",
                                 "The number of bands specified by `out_wfc_norm` or `out_wfc_re_im` in the "
                                 "INPUT file exceeds `nbands`!");
    }

    // Check if all elements in bands_picked are 0 or 1
    for (int value: out_wfc_norm)
    {
        if (value != 0 && value != 1)
        {
            ModuleBase::WARNING_QUIT("ModuleIO::get_wf_pw",
                                     "The elements of `out_wfc_norm` must be either 0 or 1. "
                                     "Invalid values found!");
        }
    }
    for (int value: out_wfc_re_im)
    {
        if (value != 0 && value != 1)
        {
            ModuleBase::WARNING_QUIT("ModuleIO::get_wf_pw",
                                     "The elements of `out_wfc_re_im` must be either 0 or 1. "
                                     "Invalid values found!");
        }
    }

    // Fill bands_picked with values from out_wfc_norm
    // Remaining bands are already set to 0
    int length = std::min(static_cast<int>(out_wfc_norm.size()), nbands);
    for (int i = 0; i < length; ++i)
    {
        // out_wfc_norm rely on function parse_expression
        bands_picked_norm[i] = static_cast<int>(out_wfc_norm[i]);
    }
    length = std::min(static_cast<int>(out_wfc_re_im.size()), nbands);
    for (int i = 0; i < length; ++i)
    {
        bands_picked_re_im[i] = static_cast<int>(out_wfc_re_im[i]);
    }

    std::vector<std::complex<double>> wfcr_norm(nxyz);
    std::vector<std::vector<double>> rho_band_norm(nspin, std::vector<double>(nxyz));

    for (int ib = 0; ib < nbands; ++ib)
    {
        // Skip the loop iteration if bands_picked[ib] is 0
        if (!bands_picked_norm[ib])
        {
            continue;
        }

        for (int is = 0; is < nspin; ++is)
        {
            std::fill(rho_band_norm[is].begin(), rho_band_norm[is].end(), 0.0);
        }
        for (int ik = 0; ik < nks; ++ik)
        {
            const int spin_index = isk[ik];
            std::cout << " Calculating wave function norm for band " << ib + 1 << ", k-point " << ik % (nks / nspin) + 1
                      << ", spin " << spin_index + 1 << std::endl;

            psi->fix_k(ik);
            pw_wfc->recip_to_real(ctx, &psi[0](ib, 0), wfcr_norm.data(), ik);

            // To ensure the normalization of charge density in multi-k calculation
            double wg_sum_k = 0;
            for (int ik_tmp = 0; ik_tmp < nks / nspin; ++ik_tmp)
            {
                wg_sum_k += wk[ik_tmp];
            }

            double w1 = static_cast<double>(wg_sum_k / ucell->omega);

            for (int i = 0; i < nxyz; ++i)
            {
                rho_band_norm[spin_index][i] = std::abs(wfcr_norm[i]) * std::sqrt(w1);
            }

            std::cout << " Writing cube files...";

            std::stringstream ss_file;
            ss_file << global_out_dir << "wf" << ib + 1 << "s" << spin_index + 1 << "k" << ik % (nks / nspin) + 1
                    << ".cube";

            ModuleIO::write_vdata_palgrid(pgrid,
                                          rho_band_norm[spin_index].data(),
                                          spin_index,
                                          nspin,
                                          0,
                                          ss_file.str(),
                                          0.0,
                                          ucell);

            std::cout << " Complete!" << std::endl;
        }
    }

    std::vector<std::complex<double>> wfc_re_im(nxyz);
    std::vector<std::vector<double>> rho_band_re(nspin, std::vector<double>(nxyz));
    std::vector<std::vector<double>> rho_band_im(nspin, std::vector<double>(nxyz));

    for (int ib = 0; ib < nbands; ++ib)
    {
        // Skip the loop iteration if bands_picked[ib] is 0
        if (!bands_picked_re_im[ib])
        {
            continue;
        }

        for (int is = 0; is < nspin; ++is)
        {
            std::fill(rho_band_re[is].begin(), rho_band_re[is].end(), 0.0);
            std::fill(rho_band_im[is].begin(), rho_band_im[is].end(), 0.0);
        }
        for (int ik = 0; ik < nks; ++ik)
        {
            const int spin_index = isk[ik];
            std::cout << " Calculating wave function real and imaginary part for band " << ib + 1 << ", k-point "
                      << ik % (nks / nspin) + 1 << ", spin " << spin_index + 1 << std::endl;

            psi->fix_k(ik);
            pw_wfc->recip_to_real(ctx, &psi[0](ib, 0), wfc_re_im.data(), ik);

            // To ensure the normalization of charge density in multi-k calculation
            double wg_sum_k = 0;
            for (int ik_tmp = 0; ik_tmp < nks / nspin; ++ik_tmp)
            {
                wg_sum_k += wk[ik_tmp];
            }

            double w1 = static_cast<double>(wg_sum_k / ucell->omega);

            for (int i = 0; i < nxyz; ++i)
            {
                rho_band_re[spin_index][i] = std::real(wfc_re_im[i]) * std::sqrt(w1);
                rho_band_im[spin_index][i] = std::imag(wfc_re_im[i]) * std::sqrt(w1);
            }

            std::cout << " Writing cube files...";

            std::stringstream ss_real;
            ss_real << global_out_dir << "wf" << ib + 1 << "s" << spin_index + 1 << "k" << ik % (nks / nspin) + 1
                    << "real.cube";

            ModuleIO::write_vdata_palgrid(pgrid,
                                          rho_band_re[spin_index].data(),
                                          spin_index,
                                          nspin,
                                          0,
                                          ss_real.str(),
                                          0.0,
                                          ucell);

            std::stringstream ss_imag;
            ss_imag << global_out_dir << "wf" << ib + 1 << "s" << spin_index + 1 << "k" << ik % (nks / nspin) + 1
                    << "imag.cube";

            ModuleIO::write_vdata_palgrid(pgrid,
                                          rho_band_im[spin_index].data(),
                                          spin_index,
                                          nspin,
                                          0,
                                          ss_imag.str(),
                                          0.0,
                                          ucell);

            std::cout << " Complete!" << std::endl;
        }
    }
}
} // namespace ModuleIO

#endif // GET_WF_PW_H