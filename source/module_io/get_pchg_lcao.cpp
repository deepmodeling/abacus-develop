#include "get_pchg_lcao.h"

#include "module_io/cube_io.h"
#include "source_estate/module_charge/symmetry_rho.h"
#include "source_estate/module_dm/cal_dm_psi.h"
#include "module_hamilt_lcao/module_gint/temp_gint/gint_interface.h"

Get_pchg_lcao::Get_pchg_lcao(psi::Psi<double>* psi_gamma_in, const Parallel_Orbitals* ParaV_in)
    : psi_gamma(psi_gamma_in), ParaV(ParaV_in)
{
}

Get_pchg_lcao::Get_pchg_lcao(psi::Psi<std::complex<double>>* psi_k_in, const Parallel_Orbitals* ParaV_in)
    : psi_k(psi_k_in), ParaV(ParaV_in)
{
}

Get_pchg_lcao::~Get_pchg_lcao()
{
}

// For gamma_only
void Get_pchg_lcao::begin(Gint_Gamma& gg,
                          double** rho,
                          const ModuleBase::matrix& wg,
                          const std::vector<double>& ef_all_spin,
                          const int rhopw_nrxx,
                          const std::vector<int>& out_pchg,
                          const int nbands,
                          const double nelec,
                          const int nspin,
                          const UnitCell* ucell_in,
                          const Parallel_Grid& pgrid,
                          const Grid_Driver* GridD_in,
                          const K_Vectors& kv,
                          const std::string& global_out_dir,
                          std::ofstream& ofs_running)
{
    ModuleBase::TITLE("Get_pchg_lcao", "begin");

    std::cout << " Calculate |psi(i)|^2 for selected electronic states (gamma only)." << std::endl;

    // if ucell is odd, it's correct,
    // if ucell is even, it's also correct.
    // +1.0e-8 in case like (2.999999999+1)/2
    const int fermi_band = static_cast<int>((nelec + 1) / 2 + 1.0e-8);

    prepare_get_pchg(ofs_running);

    // Set this->bands_picked_
    select_bands(out_pchg, nbands, fermi_band);

    for (int ib = 0; ib < nbands; ++ib)
    {
        if (bands_picked_[ib])
        {
            // Using new density matrix inplementation (gamma only)
            elecstate::DensityMatrix<double, double> DM(this->ParaV, nspin);

#ifdef __MPI
            this->idmatrix(ib, nspin, nelec, wg, DM, kv);
#else
            ModuleBase::WARNING_QUIT("Get_pchg_lcao::begin", "The `pchg` calculation is only available for MPI now!");
#endif

            for (int is = 0; is < nspin; ++is)
            {
                ModuleBase::GlobalFunc::ZEROS(rho[is], rhopw_nrxx);
            }

            DM.init_DMR(GridD_in, ucell_in);
            DM.cal_DMR();
#ifdef __OLD_GINT
            gg.initialize_pvpR(*ucell_in, GridD_in, nspin);
            gg.transfer_DM2DtoGrid(DM.get_DMR_vector());
            Gint_inout inout(rho, Gint_Tools::job_type::rho, nspin);
            gg.cal_gint(&inout);
#else
            ModuleGint::cal_gint_rho(DM.get_DMR_vector(), nspin, rho);
#endif

            // A solution to replace the original implementation of the following code:
            // pelec->charge->save_rho_before_sum_band();
            // Using std::vector to replace the original double** rho_save
            std::vector<std::vector<double>> rho_save(nspin, std::vector<double>(rhopw_nrxx));

            for (int is = 0; is < nspin; ++is)
            {
                ModuleBase::GlobalFunc::DCOPY(rho[is], rho_save[is].data(), rhopw_nrxx); // Copy data
            }

            for (int is = 0; is < nspin; ++is)
            {
                // ssc should be inside the inner loop to reset the string stream each time
                std::stringstream ssc;
                ssc << global_out_dir << "pchgi" << ib + 1 << "s" << is + 1 << ".cube";

                ofs_running << " Writing cube file " << ssc.str() << std::endl;

                // Use a const vector to store efermi for all spins, replace the original implementation:
                // const double ef_tmp = pelec->eferm.get_efval(is);
                double ef_spin = ef_all_spin[is];
                ModuleIO::write_vdata_palgrid(pgrid, rho_save[is].data(), is, nspin, 0, ssc.str(), ef_spin, ucell_in);
            }
        }
    }

    return;
}

// For multi-k
void Get_pchg_lcao::begin(Gint_k& gk,
                          double** rho,
                          std::complex<double>** rhog,
                          const ModuleBase::matrix& wg,
                          const std::vector<double>& ef_all_spin,
                          const ModulePW::PW_Basis* rho_pw,
                          const int rhopw_nrxx,
                          const std::vector<int>& out_pchg,
                          const int nbands,
                          const double nelec,
                          const int nspin,
                          UnitCell* ucell_in,
                          const Parallel_Grid& pgrid,
                          const Grid_Driver* GridD_in,
                          const K_Vectors& kv,
                          const std::string& global_out_dir,
                          std::ofstream& ofs_running,
                          const bool if_separate_k,
                          const int chr_ngmc)
{
    ModuleBase::TITLE("Get_pchg_lcao", "begin");

    std::cout << " Calculate |psi(i)|^2 for selected electronic states (multi-k)." << std::endl;

    const int fermi_band = static_cast<int>((nelec + 1) / 2 + 1.0e-8);

    prepare_get_pchg(ofs_running);

    // Set this->bands_picked_
    select_bands(out_pchg, nbands, fermi_band);

    for (int ib = 0; ib < nbands; ++ib)
    {
        if (bands_picked_[ib])
        {
            // Using new density matrix inplementation (multi-k)
            const int nspin_dm = std::map<int, int>({{1, 1}, {2, 2}, {4, 1}})[nspin];
            elecstate::DensityMatrix<std::complex<double>, double> DM(this->ParaV,
                                                                      nspin_dm,
                                                                      kv.kvec_d,
                                                                      kv.get_nks() / nspin_dm);

#ifdef __MPI
            this->idmatrix(ib, nspin, nelec, wg, DM, kv, if_separate_k);
#else
            ModuleBase::WARNING_QUIT("Get_pchg_lcao::begin", "The `pchg` calculation is only available for MPI now!");
#endif
            // If contribution from different k-points need to be output separately
            if (if_separate_k)
            {
                // For multi-k, loop over all real k-points
                for (int ik = 0; ik < kv.get_nks() / nspin; ++ik)
                {
                    for (int is = 0; is < nspin; ++is)
                    {
                        ModuleBase::GlobalFunc::ZEROS(rho[is], rhopw_nrxx);
                    }

                    DM.init_DMR(GridD_in, ucell_in);
                    DM.cal_DMR(ik);
#ifdef __OLD_GINT
                    gk.initialize_pvpR(*ucell_in, GridD_in, nspin);
                    gk.transfer_DM2DtoGrid(DM.get_DMR_vector());
                    Gint_inout inout(rho, Gint_Tools::job_type::rho, nspin);
                    gk.cal_gint(&inout);
#else
                    ModuleGint::cal_gint_rho(DM.get_DMR_vector(), nspin, rho);
#endif
                

                    // Using std::vector to replace the original double** rho_save
                    std::vector<std::vector<double>> rho_save(nspin, std::vector<double>(rhopw_nrxx));

                    for (int is = 0; is < nspin; ++is)
                    {
                        ModuleBase::GlobalFunc::DCOPY(rho[is], rho_save[is].data(), rhopw_nrxx); // Copy data
                    }

                    for (int is = 0; is < nspin; ++is)
                    {
                        // ssc should be inside the inner loop to reset the string stream each time
                        std::stringstream ssc;
                        ssc << global_out_dir << "pchgi" << ib + 1 << "s" << is + 1 << "k" << ik + 1 << ".cube";

                        ofs_running << " Writing cube file " << ssc.str() << std::endl;

                        double ef_spin = ef_all_spin[is];
                        ModuleIO::write_vdata_palgrid(pgrid,
                                                      rho_save[is].data(),
                                                      is,
                                                      nspin,
                                                      0,
                                                      ssc.str(),
                                                      ef_spin,
                                                      ucell_in);
                    }
                }
            }
            else
            {
                for (int is = 0; is < nspin; ++is)
                {
                    ModuleBase::GlobalFunc::ZEROS(rho[is], rhopw_nrxx);
                }

                DM.init_DMR(GridD_in, ucell_in);
                DM.cal_DMR();
#ifdef __OLD_GINT
                gk.initialize_pvpR(*ucell_in, GridD_in, nspin);
                gk.transfer_DM2DtoGrid(DM.get_DMR_vector());
                Gint_inout inout(rho, Gint_Tools::job_type::rho, nspin);
                gk.cal_gint(&inout);
#else
                ModuleGint::cal_gint_rho(DM.get_DMR_vector(), nspin, rho);
#endif
                // Using std::vector to replace the original double** rho_save
                std::vector<std::vector<double>> rho_save(nspin, std::vector<double>(rhopw_nrxx));

                for (int is = 0; is < nspin; ++is)
                {
                    ModuleBase::GlobalFunc::DCOPY(rho[is], rho_save[is].data(), rhopw_nrxx); // Copy data
                }

                // Symmetrize the charge density, otherwise the results are incorrect if the symmetry is on
                Symmetry_rho srho;
                for (int is = 0; is < nspin; ++is)
                {
                    std::vector<double*> rho_save_pointers(nspin);
                    for (int i = 0; i < nspin; ++i)
                    {
                        rho_save_pointers[i] = rho_save[i].data();
                    }
                    srho.begin(is, rho_save_pointers.data(), rhog, chr_ngmc, nullptr, rho_pw, ucell_in->symm);
                }

                for (int is = 0; is < nspin; ++is)
                {
                    // ssc should be inside the inner loop to reset the string stream each time
                    std::stringstream ssc;
                    ssc << global_out_dir << "pchgi" << ib + 1 << "s" << is + 1 << ".cube";

                    ofs_running << " Writing cube file " << ssc.str() << std::endl;

                    double ef_spin = ef_all_spin[is];
                    ModuleIO::write_vdata_palgrid(pgrid,
                                                  rho_save[is].data(),
                                                  is,
                                                  nspin,
                                                  0,
                                                  ssc.str(),
                                                  ef_spin,
                                                  ucell_in);
                }
            }
        }
    }

    return;
}

void Get_pchg_lcao::select_bands(const std::vector<int>& out_pchg, const int nbands, const int fermi_band)
{
    ModuleBase::TITLE("Get_pchg_lcao", "select_bands");

    int bands_below = 0;
    int bands_above = 0;

    this->bands_picked_.resize(nbands);
    ModuleBase::GlobalFunc::ZEROS(bands_picked_.data(), nbands);

    // Select bands directly using parameter `out_pchg`
    // Check if length of out_pchg is valid
    if (static_cast<int>(out_pchg.size()) > nbands)
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
    // Fill bands_picked_ with values from out_pchg
    // Remaining bands are already set to 0
    const int length = std::min(static_cast<int>(out_pchg.size()), nbands);
    std::copy(out_pchg.begin(), out_pchg.begin() + length, bands_picked_.begin());

    // Check if there are selected bands below the Fermi surface
    bool has_below = false;
    for (int i = 0; i + 1 <= fermi_band; ++i)
    {
        if (bands_picked_[i] == 1)
        {
            has_below = true;
            break;
        }
    }
    if (has_below)
    {
        std::cout << " Plot band-decomposed charge densities below the Fermi surface: band ";
        for (int i = 0; i + 1 <= fermi_band; ++i)
        {
            if (bands_picked_[i] == 1)
            {
                std::cout << i + 1 << " ";
            }
        }
        std::cout << std::endl;
    }

    // Check if there are selected bands above the Fermi surface
    bool has_above = false;
    for (int i = fermi_band; i < nbands; ++i)
    {
        if (bands_picked_[i] == 1)
        {
            has_above = true;
            break;
        }
    }
    if (has_above)
    {
        std::cout << " Plot band-decomposed charge densities above the Fermi surface: band ";
        for (int i = fermi_band; i < nbands; ++i)
        {
            if (bands_picked_[i] == 1)
            {
                std::cout << i + 1 << " ";
            }
        }
        std::cout << std::endl;
    }
}

#ifdef __MPI
// For gamma_only
void Get_pchg_lcao::idmatrix(const int& ib,
                             const int nspin,
                             const double& nelec,
                             const ModuleBase::matrix& wg,
                             elecstate::DensityMatrix<double, double>& DM,
                             const K_Vectors& kv)
{
    ModuleBase::TITLE("Get_pchg_lcao", "idmatrix");
    assert(wg.nr == nspin);

    const int fermi_band = static_cast<int>((nelec + 1) / 2 + 1.0e-8);

    for (int is = 0; is < nspin; ++is)
    {
        std::cout << " Calculating density matrix for band " << ib + 1 << ", spin " << is + 1 << std::endl;

        std::vector<double> wg_local(this->ParaV->ncol, 0.0);
        const int ib_local = this->ParaV->global2local_col(ib);

        if (ib_local >= 0)
        {
            // For unoccupied bands, use occupation of HOMO
            wg_local[ib_local] = (ib < fermi_band) ? wg(is, ib) : wg(is, fermi_band - 1);
        }

        // wg_wfc(ib,iw) = wg[ib] * wfc(ib,iw);
        this->psi_gamma->fix_k(is);

        // psi::Psi<double> wg_wfc(*this->psi_gamma, 1, this->psi_gamma->get_nbands());
        psi::Psi<double> wg_wfc(1,
                                this->psi_gamma->get_nbands(),
                                this->psi_gamma->get_nbasis(),
                                this->psi_gamma->get_nbasis(),
                                true);
        wg_wfc.set_all_psi(this->psi_gamma->get_pointer(), wg_wfc.size());

        for (int ir = 0; ir < wg_wfc.get_nbands(); ++ir)
        {
            BlasConnector::scal(wg_wfc.get_nbasis(), wg_local[ir], wg_wfc.get_pointer() + ir * wg_wfc.get_nbasis(), 1);
        }

        elecstate::psiMulPsiMpi(wg_wfc,
                                *(this->psi_gamma),
                                DM.get_DMK_pointer(is),
                                this->ParaV->desc_wfc,
                                this->ParaV->desc);
    }
}

// For multi-k
void Get_pchg_lcao::idmatrix(const int& ib,
                             const int nspin,
                             const double& nelec,
                             const ModuleBase::matrix& wg,
                             elecstate::DensityMatrix<std::complex<double>, double>& DM,
                             const K_Vectors& kv,
                             const bool if_separate_k)
{
    ModuleBase::TITLE("Get_pchg_lcao", "idmatrix");
    assert(wg.nr == kv.get_nks());

    const int fermi_band = static_cast<int>((nelec + 1) / 2 + 1.0e-8);

    // To ensure the normalization of charge density in multi-k calculation (if if_separate_k is true)
    double wg_sum_k = 0;
    double wg_sum_k_homo = 0;
    for (int ik = 0; ik < kv.get_nks() / nspin; ++ik)
    {
        wg_sum_k += wg(ik, ib);
        wg_sum_k_homo += wg(ik, fermi_band - 1);
    }

    for (int ik = 0; ik < kv.get_nks(); ++ik)
    {
        std::cout << " Calculating density matrix for band " << ib + 1 << ", k-point "
                  << ik % (kv.get_nks() / nspin) + 1 << ", spin " << kv.isk[ik] + 1 << std::endl;

        std::vector<double> wg_local(this->ParaV->ncol, 0.0);
        const int ib_local = this->ParaV->global2local_col(ib);

        if (ib_local >= 0)
        {
            double wg_value;
            if (if_separate_k)
            {
                wg_value = (ib < fermi_band) ? wg_sum_k : wg_sum_k_homo;
            }
            else
            {
                wg_value = (ib < fermi_band) ? wg(ik, ib) : wg(ik, fermi_band - 1);
            }
            wg_local[ib_local] = wg_value;
        }

        this->psi_k->fix_k(ik);

        psi::Psi<std::complex<double>> wg_wfc(1,
                                              this->psi_k->get_nbands(),
                                              this->psi_k->get_nbasis(),
                                              this->psi_k->get_nbasis(),
                                              true);
        wg_wfc.set_all_psi(this->psi_k->get_pointer(), wg_wfc.size());

        for (int ir = 0; ir < wg_wfc.get_nbands(); ++ir)
        {
            BlasConnector::scal(wg_wfc.get_nbasis(), wg_local[ir], wg_wfc.get_pointer() + ir * wg_wfc.get_nbasis(), 1);
        }

        elecstate::psiMulPsiMpi(wg_wfc,
                                *(this->psi_k),
                                DM.get_DMK_pointer(ik),
                                this->ParaV->desc_wfc,
                                this->ParaV->desc);
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