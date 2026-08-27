#include "get_pchg_lcao.h"

#include "source_estate/module_charge/symm_rho.h"
#include "source_estate/module_dm/cal_dm_psi.h"
#include "source_hamilt/module_gint/gint_interface.h"
#include "source_io/module_output/cube_io.h"

#include <algorithm>
#include <cassert>

Get_pchg_lcao::Get_pchg_lcao(const psi::Psi<double>& psi, const Parallel_Orbitals& para_orb, const int nspin)
    : psi_gamma_(&psi), psi_k_(nullptr), para_orb_(para_orb), nspin_(nspin), nbands_(para_orb.get_wfc_global_nbands())
{
}

Get_pchg_lcao::Get_pchg_lcao(const psi::Psi<std::complex<double>>& psi, const Parallel_Orbitals& para_orb, const int nspin)
    : psi_gamma_(nullptr), psi_k_(&psi), para_orb_(para_orb), nspin_(nspin), nbands_(para_orb.get_wfc_global_nbands())
{
}

// For gamma_only
void Get_pchg_lcao::begin_gamma(const UnitCell& ucell,
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
    const int nrxx = pgrid.get_nrxx();
    const int precision = 11;
    std::vector<std::vector<double>> rho(nspin_, std::vector<double>(nrxx));
    std::vector<double*> rho_pointers(nspin_);
    for (int is = 0; is < nspin_; ++is)
    {
        rho_pointers[is] = rho[is].data();
    }
    const std::vector<int> bands_picked = select_bands(out_pchg);

    for (int ib = 0; ib < nbands_; ++ib)
    {
        if (bands_picked[ib])
        {
            ModuleBase::matrix state_weights(nspin_, nbands_);
            const double spin_degeneracy = nspin_ == 1 ? 2.0 : 1.0;
            for (int is = 0; is < nspin_; ++is)
            {
                state_weights(is, ib) = spin_degeneracy;
            }

            elecstate::DensityMatrix<double, double> DM(&para_orb_, nspin_);
            elecstate::cal_dm_psi(&para_orb_, state_weights, *psi_gamma_, DM);

            for (int is = 0; is < nspin_; ++is)
            {
                std::fill(rho[is].begin(), rho[is].end(), 0.0);
            }

            DM.init_DMR(&grid_driver, &ucell);
            DM.cal_DMR();
            ModuleGint::cal_gint_rho(DM.get_DMR_vector(), nspin_, rho_pointers.data());

            for (int is = 0; is < nspin_; ++is)
            {
                // ssc should be inside the inner loop to reset the string stream each time
                std::stringstream ssc;
                ssc << global_out_dir << "pchgi" << ib + 1 << "s" << is + 1 << ".cube";

                ofs_running << " Writing cube file " << ssc.str() << std::endl;

                ModuleIO::write_vdata_palgrid(pgrid, rho[is].data(), is, nspin_, 0, ssc.str(), 0.0, &ucell, precision, 0, false, false);
            }
        }
    }

    return;
}

// For multi-k
void Get_pchg_lcao::begin_k(const ModulePW::PW_Basis& rho_pw,
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
    assert(pgrid.get_nrxx() == rho_pw.nrxx);
    const int nrxx = pgrid.get_nrxx();
    const int precision = 11;
    std::vector<std::vector<double>> rho(nspin_, std::vector<double>(nrxx));
    std::vector<double*> rho_pointers(nspin_);
    for (int is = 0; is < nspin_; ++is)
    {
        rho_pointers[is] = rho[is].data();
    }

    const bool needs_symmetry = !if_separate_k && ModuleSymmetry::Symmetry::symm_flag == 1;
    std::vector<std::vector<std::complex<double>>> rhog;
    std::vector<std::complex<double>*> rhog_pointers;
    if (needs_symmetry)
    {
        rhog.resize(nspin_, std::vector<std::complex<double>>(rho_pw.npw));
        rhog_pointers.resize(nspin_);
        for (int is = 0; is < nspin_; ++is)
        {
            rhog_pointers[is] = rhog[is].data();
        }
    }
    const std::vector<int> bands_picked = select_bands(out_pchg);

    for (int ib = 0; ib < nbands_; ++ib)
    {
        if (bands_picked[ib])
        {
            ModuleBase::matrix state_weights(kv.get_nks(), nbands_);
            const double spin_degeneracy = nspin_ == 1 ? 2.0 : 1.0;
            for (int ik = 0; ik < kv.get_nks(); ++ik)
            {
                state_weights(ik, ib) = if_separate_k ? spin_degeneracy : kv.wk[ik];
            }

            const int nspin_dm = nspin_ == 2 ? 2 : 1;
            const int nk_output = kv.get_nks() / nspin_dm;
            elecstate::DensityMatrix<std::complex<double>, double> DM(&para_orb_, nspin_dm, kv.kvec_d, nk_output);
            elecstate::cal_dm_psi(&para_orb_, state_weights, *psi_k_, DM);

            // If contribution from different k-points need to be output separately
            if (if_separate_k)
            {
                // For multi-k, loop over all real k-points
                for (int ik = 0; ik < nk_output; ++ik)
                {
                    for (int is = 0; is < nspin_; ++is)
                    {
                        std::fill(rho[is].begin(), rho[is].end(), 0.0);
                    }

                    DM.init_DMR(&grid_driver, &ucell);
                    DM.cal_DMR(ik);
                    ModuleGint::cal_gint_rho(DM.get_DMR_vector(), nspin_, rho_pointers.data());

                    for (int is = 0; is < nspin_; ++is)
                    {
                        // ssc should be inside the inner loop to reset the string stream each time
                        std::stringstream ssc;
                        ssc << global_out_dir << "pchgi" << ib + 1 << "s" << is + 1 << "k" << ik + 1 << ".cube";

                        ofs_running << " Writing cube file " << ssc.str() << std::endl;

                        ModuleIO::write_vdata_palgrid(pgrid,
                                                      rho[is].data(),
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
                    std::fill(rho[is].begin(), rho[is].end(), 0.0);
                }

                DM.init_DMR(&grid_driver, &ucell);
                DM.cal_DMR();
                ModuleGint::cal_gint_rho(DM.get_DMR_vector(), nspin_, rho_pointers.data());

                // Symmetrize the charge density, otherwise the results are incorrect if the symmetry is on
                if (needs_symmetry)
                {
                    Symmetry_rho srho;
                    if (nspin_ == 4)
                    {
                        srho.begin(0, rho_pointers.data(), rhog_pointers.data(), rho_pw.npw, nullptr, &rho_pw, ucell.symm);
                        srho.begin_soc(rho_pointers.data(), rhog_pointers.data(), &rho_pw, ucell.symm);
                    }
                    else
                    {
                        for (int is = 0; is < nspin_; ++is)
                        {
                            srho.begin(is, rho_pointers.data(), rhog_pointers.data(), rho_pw.npw, nullptr, &rho_pw, ucell.symm);
                        }
                    }
                }

                for (int is = 0; is < nspin_; ++is)
                {
                    // ssc should be inside the inner loop to reset the string stream each time
                    std::stringstream ssc;
                    ssc << global_out_dir << "pchgi" << ib + 1 << "s" << is + 1 << ".cube";

                    ofs_running << " Writing cube file " << ssc.str() << std::endl;

                    ModuleIO::write_vdata_palgrid(pgrid, rho[is].data(), is, nspin_, 0, ssc.str(), 0.0, &ucell, precision, 0, false, false);
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

    return bands_picked;
}

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
