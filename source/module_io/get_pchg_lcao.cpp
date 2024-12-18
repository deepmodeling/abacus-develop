#include "get_pchg_lcao.h"

#include "module_base/blas_connector.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/parallel_common.h"
#include "module_base/scalapack_connector.h"
#include "module_elecstate/module_charge/symmetry_rho.h"
#include "module_elecstate/module_dm/cal_dm_psi.h"
#include "module_elecstate/module_dm/density_matrix.h"
#include "module_hamilt_lcao/module_gint/gint.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_io/cube_io.h"

IState_Charge::IState_Charge(psi::Psi<double>* psi_gamma_in, const Parallel_Orbitals* ParaV_in)
    : psi_gamma(psi_gamma_in), ParaV(ParaV_in)
{
}

IState_Charge::IState_Charge(psi::Psi<std::complex<double>>* psi_k_in, const Parallel_Orbitals* ParaV_in)
    : psi_k(psi_k_in), ParaV(ParaV_in)
{
}

IState_Charge::~IState_Charge()
{
}

// For gamma_only
void IState_Charge::begin(Gint_Gamma& gg,
                          double** rho,
                          const ModuleBase::matrix& wg,
                          const std::vector<double>& ef_all_spin,
                          const int rhopw_nrxx,
                          const int rhopw_nplane,
                          const int rhopw_startz_current,
                          const int rhopw_nx,
                          const int rhopw_ny,
                          const int rhopw_nz,
                          const int bigpw_bz,
                          const int bigpw_nbz,
                          const bool gamma_only_local,
                          const int nbands_istate,
                          const std::vector<int>& out_pchg,
                          const int nbands,
                          const double nelec,
                          const int nspin,
                          const int nlocal,
                          const std::string& global_out_dir,
                          std::ofstream& ofs_warning,
                          const UnitCell* ucell_in,
                          const Parallel_Grid& pgrid,
                          const Grid_Driver* GridD_in,
                          const K_Vectors& kv)
{
    ModuleBase::TITLE("IState_Charge", "begin");

    std::cout << " Calculate |psi(i)|^2 for selected bands (band-decomposed charge densities, gamma only)."
              << std::endl;

    // Determine the mode based on the input parameters
    int mode = 0;
    // mode = 1: select bands below and above the Fermi surface using parameter `nbands_istate`
    if (nbands_istate > 0 && static_cast<int>(out_pchg.size()) == 0)
    {
        mode = 1;
    }
    // mode = 2: select bands directly using parameter `out_pchg`
    else if (static_cast<int>(out_pchg.size()) > 0)
    {
        // If out_pchg is not empty, set mode to 2
        mode = 2;
        std::cout << " Notice: INPUT parameter `nbands_istate` overwritten by `out_pchg`!" << std::endl;
    }

    // if ucell is odd, it's correct,
    // if ucell is even, it's also correct.
    // +1.0e-8 in case like (2.999999999+1)/2
    const int fermi_band = static_cast<int>((nelec + 1) / 2 + 1.0e-8);
    std::cout << " number of electrons = " << nelec << std::endl;
    std::cout << " number of occupied bands = " << fermi_band << std::endl;

    // Set this->bands_picked_ according to the mode
    select_bands(nbands_istate, out_pchg, nbands, nelec, mode, fermi_band);

    for (int ib = 0; ib < nbands; ++ib)
    {
        if (bands_picked_[ib])
        {
            // Using new density matrix inplementation (gamma only)
            elecstate::DensityMatrix<double, double> DM(this->ParaV, nspin);

#ifdef __MPI
            this->idmatrix(ib, nspin, nelec, nlocal, wg, DM, kv);
#else
            ModuleBase::WARNING_QUIT("IState_Charge::begin", "The `pchg` calculation is only available for MPI now!");
#endif

            for (int is = 0; is < nspin; ++is)
            {
                ModuleBase::GlobalFunc::ZEROS(rho[is], rhopw_nrxx);
            }

            std::cout << " Performing grid integral over real space grid for band " << ib + 1 << "..." << std::endl;

            DM.init_DMR(GridD_in, ucell_in);
            DM.cal_DMR();
            gg.initialize_pvpR(*ucell_in, GridD_in, PARAM.inp.nspin);
            gg.transfer_DM2DtoGrid(DM.get_DMR_vector());
            Gint_inout inout(rho, Gint_Tools::job_type::rho, PARAM.inp.nspin);
            gg.cal_gint(&inout);

            // A solution to replace the original implementation of the following code:
            // pelec->charge->save_rho_before_sum_band();
            // Using std::vector to replace the original double** rho_save
            std::vector<std::vector<double>> rho_save(nspin, std::vector<double>(rhopw_nrxx));

            for (int is = 0; is < nspin; ++is)
            {
                ModuleBase::GlobalFunc::DCOPY(rho[is], rho_save[is].data(), rhopw_nrxx); // Copy data
            }

            std::cout << " Writing cube files...";

            for (int is = 0; is < nspin; ++is)
            {
                // ssc should be inside the inner loop to reset the string stream each time
                std::stringstream ssc;
                ssc << global_out_dir << "BAND" << ib + 1 << "_GAMMA" << "_SPIN" << is + 1 << "_CHG.cube";

                // Use a const vector to store efermi for all spins, replace the original implementation:
                // const double ef_tmp = pelec->eferm.get_efval(is);
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

            std::cout << " Complete!" << std::endl;
        }
    }

    return;
}

// For multi-k
void IState_Charge::begin(Gint_k& gk,
                          double** rho,
                          std::complex<double>** rhog,
                          const ModuleBase::matrix& wg,
                          const std::vector<double>& ef_all_spin,
                          const ModulePW::PW_Basis* rho_pw,
                          const int rhopw_nrxx,
                          const int rhopw_nplane,
                          const int rhopw_startz_current,
                          const int rhopw_nx,
                          const int rhopw_ny,
                          const int rhopw_nz,
                          const int bigpw_bz,
                          const int bigpw_nbz,
                          const bool gamma_only_local,
                          const int nbands_istate,
                          const std::vector<int>& out_pchg,
                          const int nbands,
                          const double nelec,
                          const int nspin,
                          const int nlocal,
                          const std::string& global_out_dir,
                          std::ofstream& ofs_warning,
                          UnitCell* ucell_in,
                          const Parallel_Grid& pgrid,
                          const Grid_Driver* GridD_in,
                          const K_Vectors& kv,
                          const bool if_separate_k,
                          Parallel_Grid* Pgrid,
                          const int ngmc)
{
    ModuleBase::TITLE("IState_Charge", "begin");

    std::cout << " Calculate |psi(i)|^2 for selected bands (band-decomposed charge densities, multi-k)." << std::endl;

    int mode = 0;
    if (nbands_istate > 0 && static_cast<int>(out_pchg.size()) == 0)
    {
        mode = 1;
    }
    else if (static_cast<int>(out_pchg.size()) > 0)
    {
        // If out_pchg is not empty, set mode to 2
        mode = 2;
        std::cout << " Notice: INPUT parameter `nbands_istate` overwritten by `out_pchg`!" << std::endl;
    }
    else
    {
        mode = 3;
    }

    const int fermi_band = static_cast<int>((nelec + 1) / 2 + 1.0e-8);
    std::cout << " number of electrons = " << nelec << std::endl;
    std::cout << " number of occupied bands = " << fermi_band << std::endl;

    // Set this->bands_picked_ according to the mode
    select_bands(nbands_istate, out_pchg, nbands, nelec, mode, fermi_band);

    for (int ib = 0; ib < nbands; ++ib)
    {
        if (bands_picked_[ib])
        {
            // Using new density matrix inplementation (multi-k)
            const int nspin_dm = std::map<int, int>({ {1,1},{2,2},{4,1} })[nspin];
            elecstate::DensityMatrix<std::complex<double>, double> DM(this->ParaV, nspin_dm, kv.kvec_d, kv.get_nks() / nspin_dm);

#ifdef __MPI
            this->idmatrix(ib, nspin, nelec, nlocal, wg, DM, kv, if_separate_k);
#else
            ModuleBase::WARNING_QUIT("IState_Charge::begin", "The `pchg` calculation is only available for MPI now!");
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

                    std::cout << " Performing grid integral over real space grid for band " << ib + 1 << ", k-point "
                              << ik + 1 << "..." << std::endl;

                    DM.init_DMR(GridD_in, ucell_in);
                    DM.cal_DMR(ik);
                    gk.initialize_pvpR(*ucell_in, GridD_in, PARAM.inp.nspin);
                    gk.transfer_DM2DtoGrid(DM.get_DMR_vector());
                    Gint_inout inout(rho, Gint_Tools::job_type::rho, PARAM.inp.nspin);
                    gk.cal_gint(&inout);

                    // Using std::vector to replace the original double** rho_save
                    std::vector<std::vector<double>> rho_save(nspin, std::vector<double>(rhopw_nrxx));

                    for (int is = 0; is < nspin; ++is)
                    {
                        ModuleBase::GlobalFunc::DCOPY(rho[is], rho_save[is].data(), rhopw_nrxx); // Copy data
                    }

                    std::cout << " Writing cube files...";

                    for (int is = 0; is < nspin; ++is)
                    {
                        // ssc should be inside the inner loop to reset the string stream each time
                        std::stringstream ssc;
                        ssc << global_out_dir << "BAND" << ib + 1 << "_K" << ik + 1 << "_SPIN" << is + 1 << "_CHG.cube";

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

                    std::cout << " Complete!" << std::endl;
                }
            }
            else
            {
                for (int is = 0; is < nspin; ++is)
                {
                    ModuleBase::GlobalFunc::ZEROS(rho[is], rhopw_nrxx);
                }

                std::cout << " Performing grid integral over real space grid for band " << ib + 1 << "..." << std::endl;

                DM.init_DMR(GridD_in, ucell_in);
                DM.cal_DMR();
                gk.initialize_pvpR(*ucell_in, GridD_in, PARAM.inp.nspin);
                gk.transfer_DM2DtoGrid(DM.get_DMR_vector());
                Gint_inout inout(rho, Gint_Tools::job_type::rho, PARAM.inp.nspin);
                gk.cal_gint(&inout);

                // Using std::vector to replace the original double** rho_save
                std::vector<std::vector<double>> rho_save(nspin, std::vector<double>(rhopw_nrxx));

                for (int is = 0; is < nspin; ++is)
                {
                    ModuleBase::GlobalFunc::DCOPY(rho[is], rho_save[is].data(), rhopw_nrxx); // Copy data
                }

                // Symmetrize the charge density, otherwise the results are incorrect if the symmetry is on
                std::cout << " Symmetrizing band-decomposed charge density..." << std::endl;
                Symmetry_rho srho;
                for (int is = 0; is < nspin; ++is)
                {
                    std::vector<double*> rho_save_pointers(nspin);
                    for (int i = 0; i < nspin; ++i)
                    {
                        rho_save_pointers[i] = rho_save[i].data();
                    }
                    srho.begin(is, rho_save_pointers.data(), rhog, ngmc, nullptr, rho_pw, ucell_in->symm);
                }

                std::cout << " Writing cube files...";

                for (int is = 0; is < nspin; ++is)
                {
                    // ssc should be inside the inner loop to reset the string stream each time
                    std::stringstream ssc;
                    ssc << global_out_dir << "BAND" << ib + 1 << "_SPIN" << is + 1 << "_CHG.cube";

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

                std::cout << " Complete!" << std::endl;
            }
        }
    }

    return;
}

void IState_Charge::select_bands(const int nbands_istate,
                                 const std::vector<int>& out_pchg,
                                 const int nbands,
                                 const double nelec,
                                 const int mode,
                                 const int fermi_band)
{
    ModuleBase::TITLE("IState_Charge", "select_bands");

    int bands_below = 0;
    int bands_above = 0;

    this->bands_picked_.resize(nbands);
    ModuleBase::GlobalFunc::ZEROS(bands_picked_.data(), nbands);

    // mode = 1: select bands below and above the Fermi surface using parameter `nbands_istate`
    if (mode == 1)
    {
        bands_below = nbands_istate;
        bands_above = nbands_istate;

        std::cout << " Plot band-decomposed charge densities below the Fermi surface with " << bands_below << " bands."
                  << std::endl;

        std::cout << " Plot band-decomposed charge densities above the Fermi surface with " << bands_above << " bands."
                  << std::endl;

        for (int ib = 0; ib < nbands; ++ib)
        {
            if (ib >= fermi_band - bands_below)
            {
                if (ib < fermi_band + bands_above)
                {
                    bands_picked_[ib] = 1;
                }
            }
        }
    }
    // mode = 2: select bands directly using parameter `out_pchg`
    else if (mode == 2)
    {
        // Check if length of out_pchg is valid
        if (static_cast<int>(out_pchg.size()) > nbands)
        {
            ModuleBase::WARNING_QUIT("IState_Charge::select_bands",
                                     "The number of bands specified by `out_pchg` in the INPUT file exceeds `nbands`!");
        }
        // Check if all elements in out_pchg are 0 or 1
        for (int value: out_pchg)
        {
            if (value != 0 && value != 1)
            {
                ModuleBase::WARNING_QUIT("IState_Charge::select_bands",
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
    else
    {
        ModuleBase::WARNING_QUIT("IState_Charge::select_bands", "Invalid mode! Please check the code.");
    }
}

#ifdef __MPI
// For gamma_only
void IState_Charge::idmatrix(const int& ib,
                             const int nspin,
                             const double& nelec,
                             const int nlocal,
                             const ModuleBase::matrix& wg,
                             elecstate::DensityMatrix<double, double>& DM,
                             const K_Vectors& kv)
{
    ModuleBase::TITLE("IState_Charge", "idmatrix");
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
        psi::Psi<double> wg_wfc(*this->psi_gamma, 1);

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
void IState_Charge::idmatrix(const int& ib,
                             const int nspin,
                             const double& nelec,
                             const int nlocal,
                             const ModuleBase::matrix& wg,
                             elecstate::DensityMatrix<std::complex<double>, double>& DM,
                             const K_Vectors& kv,
                             const bool if_separate_k)
{
    ModuleBase::TITLE("IState_Charge", "idmatrix");
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
        psi::Psi<std::complex<double>> wg_wfc(*this->psi_k, 1);

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
#endif
