#include "get_wf_lcao.h"

#include "source_hamilt/module_gint/gint_env_gamma.h"
#include "source_hamilt/module_gint/gint_env_k.h"
#include "source_cell/cube_io.h"

#include <algorithm>
#include <cassert>
#include <cmath>

Get_wf_lcao::Get_wf_lcao(const psi::Psi<double>& psi, const Parallel_Orbitals& para_orb, const int nspin)
    : psi_gamma_(&psi), psi_k_(nullptr), para_orb_(para_orb), nspin_(nspin), nbands_(para_orb.get_wfc_global_nbands())
{
}

Get_wf_lcao::Get_wf_lcao(const psi::Psi<std::complex<double>>& psi, const Parallel_Orbitals& para_orb, const int nspin)
    : psi_gamma_(nullptr), psi_k_(&psi), para_orb_(para_orb), nspin_(nspin), nbands_(para_orb.get_wfc_global_nbands())
{
}

// For gamma_only
void Get_wf_lcao::begin_gamma(const UnitCell& ucell,
                              const Parallel_Grid& pgrid,
                              const std::vector<int>& out_wfc_norm,
                              const std::vector<int>& out_wfc_re_im,
                              const std::string& global_out_dir,
                              std::ofstream& ofs_running)
{
    ModuleBase::TITLE("Get_wf_lcao", "begin");

    prepare_get_wf(ofs_running);

    assert(psi_gamma_ != nullptr);
    const int nrxx = pgrid.get_nrxx();
    const int nlocal = para_orb_.get_wfc_global_nbasis();
    const int precision = 11;
    const std::vector<int> norm_bands_picked = select_bands(out_wfc_norm);
    const std::vector<int> re_im_bands_picked = select_bands(out_wfc_re_im);

    // LCAO output uses one global process pool after diagonalization.
    // Gamma-only LCAO coefficients and basis functions are real.
    std::vector<double> wfc_gamma(nrxx);
    std::vector<double> wfc_norm(nrxx);
    const std::vector<double> wfc_imag(nrxx, 0.0);

    for (int is = 0; is < nspin_; ++is)
    {
        psi_gamma_->fix_k(is);
        ModuleGint::Gint_env_gamma gint_env(psi_gamma_->get_pointer(), &para_orb_, nbands_, nlocal, wfc_gamma.data());
        for (int ib = 0; ib < nbands_; ++ib)
        {
            // Reconstruct a band only once when both output forms request it.
            if (!norm_bands_picked[ib] && !re_im_bands_picked[ib])
            {
                continue;
            }

            gint_env.cal_env_band(ib);

            if (norm_bands_picked[ib])
            {
                for (int ir = 0; ir < nrxx; ++ir)
                {
                    wfc_norm[ir] = std::abs(wfc_gamma[ir]);
                }

                std::stringstream ss_file;
                ss_file << "wfi" << ib + 1 << "s" << is + 1 << "k1.cube";

                std::stringstream ss_out;
                ss_out << global_out_dir << ss_file.str();

                std::stringstream ss_info;
                ss_info << "Wave func. " << ib + 1 << " spin " << is + 1 << " saved in";

                ModuleBase::GlobalFunc::OUT(ofs_running, ss_info.str(), ss_file.str());

                ModuleIO::write_vdata_palgrid(pgrid, wfc_norm.data(), is, nspin_, 0, ss_out.str(), 0.0, &ucell, precision, 0, false, false);
            }

            if (re_im_bands_picked[ib])
            {
                std::stringstream ss_real;
                ss_real << global_out_dir << "wfi" << ib + 1 << "s" << is + 1 << "k1re.cube";
                ModuleIO::write_vdata_palgrid(pgrid,
                                              wfc_gamma.data(),
                                              is,
                                              nspin_,
                                              0,
                                              ss_real.str(),
                                              0.0,
                                              &ucell,
                                              precision,
                                              0,
                                              false,
                                              false);

                std::stringstream ss_imag;
                ss_imag << global_out_dir << "wfi" << ib + 1 << "s" << is + 1 << "k1im.cube";
                ModuleIO::write_vdata_palgrid(pgrid,
                                              wfc_imag.data(),
                                              is,
                                              nspin_,
                                              0,
                                              ss_imag.str(),
                                              0.0,
                                              &ucell,
                                              precision,
                                              0,
                                              false,
                                              false);
            }
        }
    }
}

// For multi-k
void Get_wf_lcao::begin_k(const UnitCell& ucell,
                          const Parallel_Grid& pgrid,
                          const K_Vectors& kv,
                          const std::vector<int>& out_wfc_norm,
                          const std::vector<int>& out_wfc_re_im,
                          const std::string& global_out_dir,
                          std::ofstream& ofs_running)
{
    ModuleBase::TITLE("Get_wf_lcao", "begin");

    prepare_get_wf(ofs_running);

    assert(psi_k_ != nullptr);

    const int nks = kv.get_nks();
    const int nrxx = pgrid.get_nrxx();
    const int nlocal = para_orb_.get_wfc_global_nbasis();
    // Noncollinear wavefunctions contain two spinor components per spatial grid point.
    const int npol = (nspin_ == 4) ? 2 : 1;
    // Collinear spin channels reuse the same physical k-point numbering.
    const int nks_without_spin = (nspin_ == 2) ? kv.get_nkstot() / 2 : kv.get_nkstot();
    const int precision = 11;
    const std::vector<int> norm_bands_picked = select_bands(out_wfc_norm);
    const std::vector<int> re_im_bands_picked = select_bands(out_wfc_re_im);

    // LCAO k-point parallelism is temporary; output uses the restored global layout.
    std::vector<std::complex<double>> wfc_k(npol * nrxx);
    std::vector<double> wfc_norm(nrxx);
    std::vector<double> wfc_real(nrxx);
    std::vector<double> wfc_imag(nrxx);

    for (int ik = 0; ik < nks; ++ik)
    {
        const int spin_index = kv.isk[ik];
        // Map the local spin-k index to the one-based physical k-point index used in filenames.
        const int k_number = kv.ik2iktot[ik] % nks_without_spin + 1;
        psi_k_->fix_k(ik);

        ModuleGint::Gint_env_k gint_env(psi_k_->get_pointer(), &para_orb_, kv.kvec_d, nbands_, nlocal, ik, npol, wfc_k.data());

        for (int ib = 0; ib < nbands_; ++ib)
        {
            // Reconstruct a band only once when both output forms request it.
            if (!norm_bands_picked[ib] && !re_im_bands_picked[ib])
            {
                continue;
            }

            gint_env.cal_env_band(ib);

            if (norm_bands_picked[ib])
            {
                // The spinor modulus combines both components into one scalar field.
                for (int ir = 0; ir < nrxx; ++ir)
                {
                    wfc_norm[ir] = (npol == 2) ? std::sqrt(std::norm(wfc_k[ir]) + std::norm(wfc_k[nrxx + ir])) : std::abs(wfc_k[ir]);
                }

                std::stringstream ss_file;
                ss_file << "wfi" << ib + 1 << "s" << spin_index + 1 << "k" << k_number << ".cube";

                std::stringstream ss_out;
                ss_out << global_out_dir << ss_file.str();

                std::stringstream ss_info;
                ss_info << "Wave func. " << ib + 1 << " spin " << spin_index + 1 << " k-point " << k_number << " saved in";

                ModuleBase::GlobalFunc::OUT(ofs_running, ss_info.str(), ss_file.str());

                ModuleIO::write_vdata_palgrid(pgrid,
                                              wfc_norm.data(),
                                              spin_index,
                                              nspin_,
                                              0,
                                              ss_out.str(),
                                              0.0,
                                              &ucell,
                                              precision,
                                              0,
                                              false,
                                              false);
            }

            if (re_im_bands_picked[ib])
            {
                for (int ipol = 0; ipol < npol; ++ipol)
                {
                    // Spinors write their two components separately; collinear states use the spin channel.
                    const int component_index = (nspin_ == 4) ? ipol : spin_index;
                    const std::complex<double>* component_wfc = wfc_k.data() + ipol * nrxx;
                    for (int ir = 0; ir < nrxx; ++ir)
                    {
                        wfc_real[ir] = component_wfc[ir].real();
                        wfc_imag[ir] = component_wfc[ir].imag();
                    }

                    std::stringstream ss_real;
                    ss_real << global_out_dir << "wfi" << ib + 1 << "s" << component_index + 1 << "k" << k_number << "re.cube";

                    ModuleIO::write_vdata_palgrid(pgrid,
                                                  wfc_real.data(),
                                                  component_index,
                                                  nspin_,
                                                  0,
                                                  ss_real.str(),
                                                  0.0,
                                                  &ucell,
                                                  precision,
                                                  0,
                                                  false,
                                                  false);

                    std::stringstream ss_imag;
                    ss_imag << global_out_dir << "wfi" << ib + 1 << "s" << component_index + 1 << "k" << k_number << "im.cube";
                    ModuleIO::write_vdata_palgrid(pgrid,
                                                  wfc_imag.data(),
                                                  component_index,
                                                  nspin_,
                                                  0,
                                                  ss_imag.str(),
                                                  0.0,
                                                  &ucell,
                                                  precision,
                                                  0,
                                                  false,
                                                  false);
                }
            }
        }
    }
}

std::vector<int> Get_wf_lcao::select_bands(const std::vector<int>& out_wfc_kb) const
{
    ModuleBase::TITLE("Get_wf_lcao", "select_bands");

    std::vector<int> bands_picked(nbands_, 0);

    // Select bands directly using parameter `out_wfc_norm` or `out_wfc_re_im`
    // Check if length of out_wfc_kb is valid
    if (static_cast<int>(out_wfc_kb.size()) > nbands_)
    {
        ModuleBase::WARNING_QUIT("Get_wf_lcao::select_bands",
                                 "The number of bands specified by `out_wfc_norm` or `out_wfc_re_im` in the INPUT "
                                 "file exceeds `nbands`!");
    }
    // Check if all elements in out_wfc_kb are 0 or 1
    for (int value: out_wfc_kb)
    {
        if (value != 0 && value != 1)
        {
            ModuleBase::WARNING_QUIT("Get_wf_lcao::select_bands",
                                     "The elements of `out_wfc_norm` or `out_wfc_re_im` must be either 0 or 1. Invalid values found!");
        }
    }
    // Fill the selected-band mask with values from out_wfc_kb
    // Remaining bands are already set to 0
    const int length = std::min(static_cast<int>(out_wfc_kb.size()), nbands_);
    std::copy(out_wfc_kb.begin(), out_wfc_kb.begin() + length, bands_picked.begin());

    return bands_picked;
}

void Get_wf_lcao::prepare_get_wf(std::ofstream& ofs_running)
{
    ofs_running << "\n\n";
    ofs_running << " GET_WF CALCULATION BEGINS" << std::endl;

    ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    ofs_running << " |                                                                    |" << std::endl;
    ofs_running << " | Here we use real-space (r) grid integral technique to calculate    |" << std::endl;
    ofs_running << " | the electronic wave function psi(i,r) for each electronic state i. |" << std::endl;
    ofs_running << " | The |psi(i,r)|, Re[psi(i,r)], Im[psi(i,r)] are printed out using   |" << std::endl;
    ofs_running << " | numerical atomic orbitals as basis set.                            |" << std::endl;
    ofs_running << " |                                                                    |" << std::endl;
    ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;

    ofs_running << "\n\n";

    ofs_running << std::setprecision(6);
}
