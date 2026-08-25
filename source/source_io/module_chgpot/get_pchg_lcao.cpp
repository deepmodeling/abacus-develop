#include "get_pchg_lcao.h"

#include "source_estate/module_charge/symm_rho.h"
#include "source_estate/module_dm/cal_dm_psi.h"
#include "source_hamilt/module_gint/gint_interface.h"
#include "source_io/module_output/cube_io.h"

#include <algorithm>
#include <cassert>
#include <map>

Get_pchg_lcao::Get_pchg_lcao(const psi::Psi<double>& psi, const Parallel_Orbitals& para_orb, const int nspin, const double nelec)
    : psi_gamma_(&psi), psi_k_(nullptr), para_orb_(para_orb), nspin_(nspin), nbands_(para_orb.get_wfc_global_nbands()),
      fermi_band_(static_cast<int>((nelec + 1) / 2 + 1.0e-8))
{
}

Get_pchg_lcao::Get_pchg_lcao(const psi::Psi<std::complex<double>>& psi,
                             const Parallel_Orbitals& para_orb,
                             const int nspin,
                             const double nelec)
    : psi_gamma_(nullptr), psi_k_(&psi), para_orb_(para_orb), nspin_(nspin), nbands_(para_orb.get_wfc_global_nbands()),
      fermi_band_(static_cast<int>((nelec + 1) / 2 + 1.0e-8))
{
}

// For gamma_only
void Get_pchg_lcao::begin_gamma(double* const* rho,
                                const ModuleBase::matrix& wg,
                                const UnitCell& ucell,
                                const Parallel_Grid& pgrid,
                                const Grid_Driver& grid_driver,
                                const std::vector<int>& out_pchg,
                                const std::string& global_out_dir,
                                std::ofstream& ofs_running)
{
    ModuleBase::TITLE("Get_pchg_lcao", "begin");

    std::cout << " Calculate |psi(i)|^2 for selected electronic states (gamma only)." << std::endl;

    prepare_get_pchg(ofs_running);

    assert(psi_gamma_ != nullptr);
    assert(rho != nullptr);
    const int nrxx = pgrid.get_nrxx();
    std::vector<double*> rho_pointers(rho, rho + nspin_);
    const std::vector<int> bands_picked = select_bands(out_pchg);

    for (int ib = 0; ib < nbands_; ++ib)
    {
        if (bands_picked[ib])
        {
            // Using new density matrix inplementation (gamma only)
            elecstate::DensityMatrix<double, double> DM(&para_orb_, nspin_);

#ifdef __MPI
            this->idmatrix(ib, wg, DM);
#else
            ModuleBase::WARNING_QUIT("Get_pchg_lcao::begin", "The `pchg` calculation is only available for MPI now!");
#endif

            for (int is = 0; is < nspin_; ++is)
            {
                ModuleBase::GlobalFunc::ZEROS(rho[is], nrxx);
            }

            DM.init_DMR(&grid_driver, &ucell);
            DM.cal_DMR();
            ModuleGint::cal_gint_rho(DM.get_DMR_vector(), nspin_, rho_pointers.data());

            // A solution to replace the original implementation of the following code:
            // pelec->charge->save_rho_before_sum_band();
            // Using std::vector to replace the original double** rho_save
            std::vector<std::vector<double>> rho_save(nspin_, std::vector<double>(nrxx));

            for (int is = 0; is < nspin_; ++is)
            {
                ModuleBase::GlobalFunc::DCOPY(rho[is], rho_save[is].data(), nrxx); // Copy data
            }

            for (int is = 0; is < nspin_; ++is)
            {
                // ssc should be inside the inner loop to reset the string stream each time
                std::stringstream ssc;
                ssc << global_out_dir << "pchgi" << ib + 1 << "s" << is + 1 << ".cube";

                ofs_running << " Writing cube file " << ssc.str() << std::endl;

                const int precision = 6;
                ModuleIO::write_vdata_palgrid(pgrid,
                                              rho_save[is].data(),
                                              is,
                                              nspin_,
                                              0,
                                              ssc.str(),
                                              0.0,
                                              &ucell,
                                              precision,
                                              0,
                                              false,
                                              false);
            }
        }
    }

    return;
}

// For multi-k
void Get_pchg_lcao::begin_k(double* const* rho,
                            std::complex<double>* const* rhog,
                            const ModuleBase::matrix& wg,
                            const ModulePW::PW_Basis& rho_pw,
                            UnitCell& ucell,
                            const Parallel_Grid& pgrid,
                            const Grid_Driver& grid_driver,
                            const K_Vectors& kv,
                            const std::vector<int>& out_pchg,
                            const bool if_separate_k,
                            const std::string& global_out_dir,
                            std::ofstream& ofs_running)
{
    ModuleBase::TITLE("Get_pchg_lcao", "begin");

    std::cout << " Calculate |psi(i)|^2 for selected electronic states (multi-k)." << std::endl;

    prepare_get_pchg(ofs_running);

    assert(psi_k_ != nullptr);
    assert(rho != nullptr);
    assert(rhog != nullptr);
    assert(pgrid.get_nrxx() == rho_pw.nrxx);
    const int nrxx = pgrid.get_nrxx();
    std::vector<double*> rho_pointers(rho, rho + nspin_);
    std::vector<std::complex<double>*> rhog_pointers(rhog, rhog + nspin_);
    const std::vector<int> bands_picked = select_bands(out_pchg);

    for (int ib = 0; ib < nbands_; ++ib)
    {
        if (bands_picked[ib])
        {
            // Using new density matrix inplementation (multi-k)
            const int nspin_dm = std::map<int, int>({{1, 1}, {2, 2}, {4, 1}})[nspin_];
            elecstate::DensityMatrix<std::complex<double>, double> DM(&para_orb_, nspin_dm, kv.kvec_d, kv.get_nks() / nspin_dm);

#ifdef __MPI
            this->idmatrix(ib, wg, DM, kv, if_separate_k);
#else
            ModuleBase::WARNING_QUIT("Get_pchg_lcao::begin", "The `pchg` calculation is only available for MPI now!");
#endif
            // If contribution from different k-points need to be output separately
            if (if_separate_k)
            {
                // For multi-k, loop over all real k-points
                for (int ik = 0; ik < kv.get_nks() / nspin_; ++ik)
                {
                    for (int is = 0; is < nspin_; ++is)
                    {
                        ModuleBase::GlobalFunc::ZEROS(rho[is], nrxx);
                    }

                    DM.init_DMR(&grid_driver, &ucell);
                    DM.cal_DMR(ik);
                    ModuleGint::cal_gint_rho(DM.get_DMR_vector(), nspin_, rho_pointers.data());

                    // Using std::vector to replace the original double** rho_save
                    std::vector<std::vector<double>> rho_save(nspin_, std::vector<double>(nrxx));

                    for (int is = 0; is < nspin_; ++is)
                    {
                        ModuleBase::GlobalFunc::DCOPY(rho[is], rho_save[is].data(), nrxx); // Copy data
                    }

                    for (int is = 0; is < nspin_; ++is)
                    {
                        // ssc should be inside the inner loop to reset the string stream each time
                        std::stringstream ssc;
                        ssc << global_out_dir << "pchgi" << ib + 1 << "s" << is + 1 << "k" << ik + 1 << ".cube";

                        ofs_running << " Writing cube file " << ssc.str() << std::endl;

                        const int precision = 6;
                        ModuleIO::write_vdata_palgrid(pgrid,
                                                      rho_save[is].data(),
                                                      is,
                                                      nspin_,
                                                      0,
                                                      ssc.str(),
                                                      0.0,
                                                      &ucell,
                                                      precision,
                                                      0,
                                                      false,
                                                      false);
                    }
                }
            }
            else
            {
                for (int is = 0; is < nspin_; ++is)
                {
                    ModuleBase::GlobalFunc::ZEROS(rho[is], nrxx);
                }

                DM.init_DMR(&grid_driver, &ucell);
                DM.cal_DMR();
                ModuleGint::cal_gint_rho(DM.get_DMR_vector(), nspin_, rho_pointers.data());
                // Using std::vector to replace the original double** rho_save
                std::vector<std::vector<double>> rho_save(nspin_, std::vector<double>(nrxx));

                for (int is = 0; is < nspin_; ++is)
                {
                    ModuleBase::GlobalFunc::DCOPY(rho[is], rho_save[is].data(), nrxx); // Copy data
                }

                // Symmetrize the charge density, otherwise the results are incorrect if the symmetry is on
                Symmetry_rho srho;
                for (int is = 0; is < nspin_; ++is)
                {
                    std::vector<double*> rho_save_pointers(nspin_);
                    for (int i = 0; i < nspin_; ++i)
                    {
                        rho_save_pointers[i] = rho_save[i].data();
                    }
                    srho.begin(is, rho_save_pointers.data(), rhog_pointers.data(), rho_pw.npw, nullptr, &rho_pw, ucell.symm);
                }

                for (int is = 0; is < nspin_; ++is)
                {
                    // ssc should be inside the inner loop to reset the string stream each time
                    std::stringstream ssc;
                    ssc << global_out_dir << "pchgi" << ib + 1 << "s" << is + 1 << ".cube";

                    ofs_running << " Writing cube file " << ssc.str() << std::endl;

                    const int precision = 6;
                    ModuleIO::write_vdata_palgrid(pgrid,
                                                  rho_save[is].data(),
                                                  is,
                                                  nspin_,
                                                  0,
                                                  ssc.str(),
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

    return;
}

std::vector<int> Get_pchg_lcao::select_bands(const std::vector<int>& out_pchg) const
{
    ModuleBase::TITLE("Get_pchg_lcao", "select_bands");

    std::vector<int> bands_picked(nbands_, 0);

    // Select bands directly using parameter `out_pchg`
    // Check if length of out_pchg is valid
    if (static_cast<int>(out_pchg.size()) > nbands_)
    {
        ModuleBase::WARNING_QUIT("Get_pchg_lcao::select_bands",
                                 "The number of bands specified by `out_pchg` in the INPUT file exceeds `nbands`!");
    }
    // Check if all elements in out_pchg are 0 or 1
    for (int value: out_pchg)
    {
        if (value != 0 && value != 1)
        {
            ModuleBase::WARNING_QUIT("Get_pchg_lcao::select_bands",
                                     "The elements of `out_pchg` must be either 0 or 1. Invalid values found!");
        }
    }
    // Fill the selected-band mask with values from out_pchg
    // Remaining bands are already set to 0
    const int length = std::min(static_cast<int>(out_pchg.size()), nbands_);
    std::copy(out_pchg.begin(), out_pchg.begin() + length, bands_picked.begin());

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
        std::cout << " Plot band-decomposed charge densities below the Fermi surface: band ";
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
        std::cout << " Plot band-decomposed charge densities above the Fermi surface: band ";
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

#ifdef __MPI
// For gamma_only
void Get_pchg_lcao::idmatrix(const int& ib, const ModuleBase::matrix& wg, elecstate::DensityMatrix<double, double>& DM)
{
    ModuleBase::TITLE("Get_pchg_lcao", "idmatrix");
    assert(wg.nr == nspin_);

    for (int is = 0; is < nspin_; ++is)
    {
        std::cout << " Calculating density matrix for band " << ib + 1 << ", spin " << is + 1 << std::endl;

        std::vector<double> wg_local(para_orb_.ncol, 0.0);
        const int ib_local = para_orb_.global2local_col(ib);

        if (ib_local >= 0)
        {
            // For unoccupied bands, use occupation of HOMO
            wg_local[ib_local] = (ib < fermi_band_) ? wg(is, ib) : wg(is, fermi_band_ - 1);
        }

        // wg_wfc(ib,iw) = wg[ib] * wfc(ib,iw);
        psi_gamma_->fix_k(is);

        // psi::Psi<double> wg_wfc(*psi_gamma_, 1, psi_gamma_->get_nbands());
        psi::Psi<double> wg_wfc(1, psi_gamma_->get_nbands(), psi_gamma_->get_nbasis(), psi_gamma_->get_nbasis(), true);
        wg_wfc.set_all_psi(psi_gamma_->get_pointer(), wg_wfc.size());

        for (int ir = 0; ir < wg_wfc.get_nbands(); ++ir)
        {
            BlasConnector::scal(wg_wfc.get_nbasis(), wg_local[ir], wg_wfc.get_pointer() + ir * wg_wfc.get_nbasis(), 1);
        }

        elecstate::psiMulPsiMpi(wg_wfc, *psi_gamma_, DM.get_DMK_pointer(is), para_orb_.desc_wfc, para_orb_.desc);
    }
}

// For multi-k
void Get_pchg_lcao::idmatrix(const int& ib,
                             const ModuleBase::matrix& wg,
                             elecstate::DensityMatrix<std::complex<double>, double>& DM,
                             const K_Vectors& kv,
                             const bool if_separate_k)
{
    ModuleBase::TITLE("Get_pchg_lcao", "idmatrix");
    assert(wg.nr == kv.get_nks());

    // To ensure the normalization of charge density in multi-k calculation (if if_separate_k is true)
    double wg_sum_k = 0;
    double wg_sum_k_homo = 0;
    for (int ik = 0; ik < kv.get_nks() / nspin_; ++ik)
    {
        wg_sum_k += wg(ik, ib);
        wg_sum_k_homo += wg(ik, fermi_band_ - 1);
    }

    for (int ik = 0; ik < kv.get_nks(); ++ik)
    {
        std::cout << " Calculating density matrix for band " << ib + 1 << ", k-point " << ik % (kv.get_nks() / nspin_) + 1 << ", spin "
                  << kv.isk[ik] + 1 << std::endl;

        std::vector<double> wg_local(para_orb_.ncol, 0.0);
        const int ib_local = para_orb_.global2local_col(ib);

        if (ib_local >= 0)
        {
            double wg_value = 0.0;
            if (if_separate_k)
            {
                wg_value = (ib < fermi_band_) ? wg_sum_k : wg_sum_k_homo;
            }
            else
            {
                wg_value = (ib < fermi_band_) ? wg(ik, ib) : wg(ik, fermi_band_ - 1);
            }
            wg_local[ib_local] = wg_value;
        }

        psi_k_->fix_k(ik);

        psi::Psi<std::complex<double>> wg_wfc(1, psi_k_->get_nbands(), psi_k_->get_nbasis(), psi_k_->get_nbasis(), true);
        wg_wfc.set_all_psi(psi_k_->get_pointer(), wg_wfc.size());
        std::complex<double>* wg_wfc_ptr = wg_wfc.get_pointer();
        for (int i = 0; i < wg_wfc.size(); ++i)
        {
            wg_wfc_ptr[i] = std::conj(wg_wfc_ptr[i]);
        }

        for (int ir = 0; ir < wg_wfc.get_nbands(); ++ir)
        {
            BlasConnector::scal(wg_wfc.get_nbasis(), wg_local[ir], wg_wfc.get_pointer() + ir * wg_wfc.get_nbasis(), 1);
        }

        elecstate::psiMulPsiMpi(wg_wfc, *psi_k_, DM.get_DMK_pointer(ik), para_orb_.desc_wfc, para_orb_.desc);
    }
}
#endif // __MPI

void Get_pchg_lcao::prepare_get_pchg(std::ofstream& ofs_running)
{
    ofs_running << "\n\n";
    ofs_running << " GET_PCHG CALCULATION BEGINS" << std::endl;

    ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    ofs_running << " |                                                                    |" << std::endl;
    ofs_running << " |  Here we use real-space (r) grid integral technique to calculate   |" << std::endl;
    ofs_running << " |  the decomposed charge density |psi(i,r)|^2 for each electronic    |" << std::endl;
    ofs_running << " |  state i. The |psi(i,r)|^2 is printed out using numerical atomic   |" << std::endl;
    ofs_running << " |  orbitals as basis set.                                            |" << std::endl;
    ofs_running << " |                                                                    |" << std::endl;
    ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;

    ofs_running << "\n\n";

    ofs_running << std::setprecision(6);
}
