#include "get_wf_lcao.h"

#include "source_hamilt/module_gint/gint_env_gamma.h"
#include "source_hamilt/module_gint/gint_env_k.h"
#include "source_io/module_output/cube_io.h"

#include <algorithm>
#include <cassert>
#include <cmath>

Get_wf_lcao::Get_wf_lcao(const psi::Psi<double>& psi, const Parallel_Orbitals& para_orb, const int nspin, const double nelec)
    : psi_gamma_(&psi), psi_k_(nullptr), para_orb_(para_orb), nspin_(nspin), nbands_(para_orb.get_wfc_global_nbands()),
      fermi_band_(static_cast<int>((nelec + 1) / 2 + 1.0e-8))
{
}

Get_wf_lcao::Get_wf_lcao(const psi::Psi<std::complex<double>>& psi, const Parallel_Orbitals& para_orb, const int nspin, const double nelec)
    : psi_gamma_(nullptr), psi_k_(&psi), para_orb_(para_orb), nspin_(nspin), nbands_(para_orb.get_wfc_global_nbands()),
      fermi_band_(static_cast<int>((nelec + 1) / 2 + 1.0e-8))
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
    std::vector<double> wfc_gamma(nrxx);

    const std::vector<int> bands_picked = this->select_bands(out_wfc_norm);

    // Calculate out_wfc_norm
    std::vector<double> wfc_norm(nrxx);
    for (int is = 0; is < nspin_; ++is)
    {
        psi_gamma_->fix_k(is);
        ModuleGint::Gint_env_gamma gint_env(psi_gamma_->get_pointer(), &para_orb_, nbands_, nlocal, wfc_gamma.data());
        for (int ib = 0; ib < nbands_; ++ib)
        {
            if (bands_picked[ib])
            {
                gint_env.cal_env_band(ib);
                for (int ir = 0; ir < nrxx; ++ir)
                {
                    wfc_norm[ir] = std::abs(wfc_gamma[ir]);
                }

                // pint out information
                std::stringstream ss_file;
                ss_file << "wfi" << ib + 1 << "s" << is + 1 << "k1.cube";

                std::stringstream ss_out;
                ss_out << global_out_dir << ss_file.str();

                std::stringstream ss_info;
                ss_info << "Wave func. " << ib + 1 << " spin " << is + 1 << " saved in";

                ModuleBase::GlobalFunc::OUT(ofs_running, ss_info.str(), ss_file.str());

                ModuleIO::write_vdata_palgrid(pgrid,
                                              wfc_norm.data(),
                                              is,
                                              nspin_,
                                              0,
                                              ss_out.str(),
                                              0.0,
                                              &(ucell),
                                              11, // default precision
                                              0,
                                              false,
                                              false);
            }
        }
    }

    const std::vector<int> re_im_bands_picked = this->select_bands(out_wfc_re_im);

    // Calculate out_wfc_re_im
    const std::vector<double> wfc_imag(nrxx, 0.0);
    for (int is = 0; is < nspin_; ++is)
    {
        psi_gamma_->fix_k(is);
        ModuleGint::Gint_env_gamma gint_env(psi_gamma_->get_pointer(), &para_orb_, nbands_, nlocal, wfc_gamma.data());
        for (int ib = 0; ib < nbands_; ++ib)
        {
            if (re_im_bands_picked[ib])
            {
                gint_env.cal_env_band(ib);

                // Output real part
                std::stringstream ss_real;
                ss_real << global_out_dir << "wfi" << ib + 1 << "s" << is + 1 << "k1re.cube";
                ModuleIO::write_vdata_palgrid(pgrid, wfc_gamma.data(), is, nspin_, 0, ss_real.str(), 0.0, &(ucell), 11, 0, false, false);

                // Output imaginary part
                std::stringstream ss_imag;
                ss_imag << global_out_dir << "wfi" << ib + 1 << "s" << is + 1 << "k1im.cube";
                ModuleIO::write_vdata_palgrid(pgrid, wfc_imag.data(), is, nspin_, 0, ss_imag.str(), 0.0, &(ucell), 11, 0, false, false);
            }
        }
    }

    return;
}

// For multi-k
void Get_wf_lcao::begin_k(double* const* rho,
                          const ModulePW::PW_Basis_K& pw_wfc,
                          const UnitCell& ucell,
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
    assert(rho != nullptr);
    assert(pgrid.get_nrxx() == pw_wfc.nrxx);

    const int nks = kv.get_nks();
    const int nlocal = para_orb_.get_wfc_global_nbasis();
    const int npol = (nspin_ == 4) ? 2 : 1;

    // for pw_wfc in G space
    psi::Psi<std::complex<double>> psi_g;

    // Reciprocal-space storage used by the current real/imaginary Cube path.
    psi_g.resize(nks, nbands_, pw_wfc.npwk_max);

    const std::vector<int> bands_picked = this->select_bands(out_wfc_norm);

    // Calculate out_wfc_norm
    for (int ik = 0; ik < nks; ++ik)
    {
        const int ispin = kv.isk[ik];
        psi_k_->fix_k(ik);

        ModuleGint::Gint_env_k
            gint_env(psi_k_->get_pointer(), &para_orb_, kv.kvec_c, kv.kvec_d, nbands_, nlocal, ik, nspin_, npol, rho[ispin]);

        for (int ib = 0; ib < nbands_; ++ib)
        {
            if (bands_picked[ib])
            {
                gint_env.cal_env_band(ib);

                // ik0 is the real k-point index, starting from 0
                int ik0 = kv.ik2iktot[ik];
                if (nspin_ == 2)
                {
                    const int half_k = kv.get_nkstot() / 2;
                    if (ik0 >= half_k)
                    {
                        ik0 -= half_k;
                    }
                }

                // pint out information
                std::stringstream ss_file;
                ss_file << "wfi" << ib + 1 << "s" << ispin + 1 << "k" << ik0 + 1 << ".cube";

                std::stringstream ss_out;
                ss_out << global_out_dir << ss_file.str();

                std::stringstream ss_info;
                ss_info << "Wave func. " << ib + 1 << " spin " << ispin + 1 << " k-point " << ik0 + 1 << " saved in";

                ModuleBase::GlobalFunc::OUT(ofs_running, ss_info.str(), ss_file.str());

                ModuleIO::write_vdata_palgrid(pgrid, rho[ispin], ispin, nspin_, 0, ss_out.str(), 0.0, &(ucell), 3, 0, false, false);

                // Prepare the reciprocal-space buffer for real/imaginary Cube output.
                psi_g.fix_k(ik);
                this->set_pw_wfc(pw_wfc, ik, ib, rho, psi_g);
            }
        }
    }

    const std::vector<int> re_im_bands_picked = this->select_bands(out_wfc_re_im);

    // Calculate out_wfc_re_im
    for (int ib = 0; ib < nbands_; ++ib)
    {
        if (re_im_bands_picked[ib])
        {
            for (int ik = 0; ik < nks; ++ik)
            {
                const int ispin = kv.isk[ik];

                psi_g.fix_k(ik);

                // Calculate real-space wave functions
                std::vector<std::complex<double>> wfc_r(pw_wfc.nrxx);
                pw_wfc.recip2real(&psi_g(ib, 0), wfc_r.data(), ik);

                // Extract real and imaginary parts
                std::vector<double> wfc_real(pw_wfc.nrxx);
                std::vector<double> wfc_imag(pw_wfc.nrxx);
                for (int ir = 0; ir < pw_wfc.nrxx; ++ir)
                {
                    wfc_real[ir] = wfc_r[ir].real();
                    wfc_imag[ir] = wfc_r[ir].imag();
                }

                // ik0 is the real k-point index, starting from 0
                int ik0 = kv.ik2iktot[ik];
                if (nspin_ == 2)
                {
                    const int half_k = kv.get_nkstot() / 2;
                    if (ik0 >= half_k)
                    {
                        ik0 -= half_k;
                    }
                }

                // Output real part
                std::stringstream ss_real;
                ss_real << global_out_dir << "wfi" << ib + 1 << "s" << ispin + 1 << "k" << ik0 + 1 << "re.cube";

                ModuleIO::write_vdata_palgrid(pgrid, wfc_real.data(), ispin, nspin_, 0, ss_real.str(), 0.0, &(ucell), 11, 0, false, false);

                // Output imaginary part
                std::stringstream ss_imag;
                ss_imag << global_out_dir << "wfi" << ib + 1 << "s" << ispin + 1 << "k" << ik0 + 1 << "im.cube";
                ModuleIO::write_vdata_palgrid(pgrid, wfc_imag.data(), ispin, nspin_, 0, ss_imag.str(), 0.0, &(ucell), 11, 0, false, false);
            }
        }
    }
    return;
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

    // Check if there are selected bands below the Fermi surface
    bool has_below = false;
    for (int i = 0; i + 1 <= fermi_band_; ++i)
    {
        if (bands_picked[i] == 1)
        {
            has_below = true;
            break;
        }
    }
    if (has_below)
    {
        std::cout << " Plot wave functions below the Fermi surface: band ";
        for (int i = 0; i + 1 <= fermi_band_; ++i)
        {
            if (bands_picked[i] == 1)
            {
                std::cout << i + 1 << " ";
            }
        }
        std::cout << std::endl;
    }

    // Check if there are selected bands above the Fermi surface
    bool has_above = false;
    for (int i = fermi_band_; i < nbands_; ++i)
    {
        if (bands_picked[i] == 1)
        {
            has_above = true;
            break;
        }
    }
    if (has_above)
    {
        std::cout << " Plot wave functions above the Fermi surface: band ";
        for (int i = fermi_band_; i < nbands_; ++i)
        {
            if (bands_picked[i] == 1)
            {
                std::cout << i + 1 << " ";
            }
        }
        std::cout << std::endl;
    }

    return bands_picked;
}

// for each band
void Get_wf_lcao::set_pw_wfc(const ModulePW::PW_Basis_K& pw_wfc,
                             const int ik,
                             const int ib,
                             const double* const* const rho,
                             psi::Psi<std::complex<double>>& wfc_g)
{
    if (ib == 0)
    {
        // once is enough
        ModuleBase::TITLE("Get_wf_lcao", "set_pw_wfc");
    }

    std::vector<std::complex<double>> Porter(pw_wfc.nrxx);
    // here I refer to v_hartree, but I don't know how to deal with NSPIN=4
    const int nspin0 = (nspin_ == 2) ? 2 : 1;
    for (int is = 0; is < nspin0; ++is)
    {
        for (int ir = 0; ir < pw_wfc.nrxx; ++ir)
        {
            Porter[ir] += std::complex<double>(rho[is][ir], 0.0);
        }
    }

    // call FFT
    pw_wfc.real2recip(Porter.data(), &wfc_g(ib, 0), ik);
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
