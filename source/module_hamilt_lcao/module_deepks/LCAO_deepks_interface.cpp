#ifdef __DEEPKS
#include "LCAO_deepks_interface.h"

#include "LCAO_deepks_io.h" // mohan add 2024-07-22
#include "module_base/global_variable.h"
#include "module_base/tool_title.h"
#include "module_elecstate/cal_dm.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"
#include "module_parameter/parameter.h"

template <typename TK, typename TR>
LCAO_Deepks_Interface<TK, TR>::LCAO_Deepks_Interface(std::shared_ptr<LCAO_Deepks> ld_in) : ld(ld_in)
{
}

template <typename TK, typename TR>
void LCAO_Deepks_Interface<TK, TR>::out_deepks_labels(const double& etot,
                                                      const int& nks,
                                                      const int& nat,
                                                      const int& nlocal,
                                                      const ModuleBase::matrix& ekb,
                                                      const std::vector<ModuleBase::Vector3<double>>& kvec_d,
                                                      const UnitCell& ucell,
                                                      const LCAO_Orbitals& orb,
                                                      const Grid_Driver& GridD,
                                                      const Parallel_Orbitals* ParaV,
                                                      const psi::Psi<TK>& psi,
                                                      const elecstate::DensityMatrix<TK, double>* dm,
                                                      hamilt::HamiltLCAO<TK, TR>* p_ham)
{
    ModuleBase::TITLE("LCAO_Deepks_Interface", "out_deepks_labels");
    ModuleBase::timer::tick("LCAO_Deepks_Interface", "out_deepks_labels");

    // define TH for different types
    using TH = std::conditional_t<std::is_same<TK, double>::value, ModuleBase::matrix, ModuleBase::ComplexMatrix>;

    const int my_rank = GlobalV::MY_RANK;
    const int nspin = PARAM.inp.nspin;

    // calculating deepks correction to bandgap and save the results
    if (PARAM.inp.deepks_out_labels)
    {
        // mohan updated 2024-07-25
        const std::string file_etot = PARAM.globalv.global_out_dir + "deepks_etot.npy";
        const std::string file_ebase = PARAM.globalv.global_out_dir + "deepks_ebase.npy";

        LCAO_deepks_io::save_npy_e(etot, file_etot, my_rank);

        if (PARAM.inp.deepks_scf)
        {
            /// ebase :no deepks E_delta including
            LCAO_deepks_io::save_npy_e(etot - ld->E_delta, file_ebase, my_rank);
        }
        else // deepks_scf = 0; base calculation
        {
            /// no scf, e_tot=e_base
            LCAO_deepks_io::save_npy_e(etot, file_ebase, my_rank);
        }

        std::vector<std::vector<TK>>* h_delta = nullptr;
        if constexpr (std::is_same<TK, double>::value)
        {
            h_delta = &ld->H_V_delta;
        }
        else
        {
            h_delta = &ld->H_V_delta_k;
        }

        if (PARAM.inp.deepks_bandgap)
        {
            const int nocc = (PARAM.inp.nelec + 1) / 2;
            std::vector<double> o_tot(nks);
            for (int iks = 0; iks < nks; ++iks)
            {
                // record band gap for each k point (including spin)
                o_tot[iks] = ekb(iks, nocc) - ekb(iks, nocc - 1);
            }

            const std::string file_otot = PARAM.globalv.global_out_dir + "deepks_otot.npy";
            LCAO_deepks_io::save_npy_o(o_tot, file_otot, nks, my_rank);

            if (PARAM.inp.deepks_scf)
            {
                ModuleBase::matrix wg_hl;
                std::vector<TH> dm_bandgap;

                // Calculate O_delta
                if constexpr (std::is_same<TK, double>::value) // for gamma only
                {
                    wg_hl.create(nspin, PARAM.inp.nbands);
                    dm_bandgap.resize(nspin);
                    for (int is = 0; is < nspin; ++is)
                    {
                        wg_hl.zero_out();
                        wg_hl(is, nocc - 1) = -1.0;
                        wg_hl(is, nocc) = 1.0;
                        elecstate::cal_dm(ParaV, wg_hl, psi, dm_bandgap);
                    }
                }
                else // for multi-k
                {
                    wg_hl.create(nks, PARAM.inp.nbands);
                    dm_bandgap.resize(nks);
                    wg_hl.zero_out();
                    for (int ik = 0; ik < nks; ik++)
                    {
                        wg_hl(ik, nocc - 1) = -1.0;
                        wg_hl(ik, nocc) = 1.0;
                    }
                    elecstate::cal_dm(ParaV, wg_hl, psi, dm_bandgap);
                }

                std::vector<double> o_delta(nks, 0.0);

                // calculate and save orbital_precalc: [nks,NAt,NDscrpt]
                torch::Tensor orbital_precalc;
                std::vector<torch::Tensor> gevdm;
                ld->cal_gevdm(nat, gevdm);
                DeePKS_domain::cal_orbital_precalc<TK, TH>(dm_bandgap,
                                                           ld->lmaxd,
                                                           ld->inlmax,
                                                           nat,
                                                           nks,
                                                           ld->inl_l,
                                                           kvec_d,
                                                           ld->phialpha,
                                                           gevdm,
                                                           ld->inl_index,
                                                           ucell,
                                                           orb,
                                                           *ParaV,
                                                           GridD,
                                                           orbital_precalc);
                DeePKS_domain::cal_o_delta<TK, TH>(dm_bandgap, *h_delta, o_delta, *ParaV, nks);

                // save obase and orbital_precalc
                LCAO_deepks_io::save_npy_orbital_precalc(nat,
                                                         nks,
                                                         ld->des_per_atom,
                                                         orbital_precalc,
                                                         PARAM.globalv.global_out_dir,
                                                         my_rank);
                const std::string file_obase = PARAM.globalv.global_out_dir + "deepks_obase.npy";
                std::vector<double> o_base(nks);
                for (int iks = 0; iks < nks; ++iks)
                {
                    o_base[iks] = o_tot[iks] - o_delta[iks];
                }
                LCAO_deepks_io::save_npy_o(o_base, file_obase, nks, my_rank);
            }    // end deepks_scf == 1
            else // deepks_scf == 0
            {
                const std::string file_obase = PARAM.globalv.global_out_dir + "deepks_obase.npy";
                LCAO_deepks_io::save_npy_o(o_tot, file_obase, nks, my_rank); // no scf, o_tot=o_base
            }                                                                // end deepks_scf == 0
        }                                                                    // end bandgap label

        // save H(R) matrix
        if (true) // should be modified later!
        {
            const std::string file_hr = PARAM.globalv.global_out_dir + "deepks_hr.npy";
            const hamilt::HContainer<TR>& hR = *(p_ham->getHR());

            // How to save H(R)?
        }

        if (PARAM.inp.deepks_v_delta)
        {
            std::vector<TH> h_tot(nks);
            std::vector<std::vector<TK>> h_mat(nks, std::vector<TK>(ParaV->nloc));
            for (int ik = 0; ik < nks; ik++)
            {
                h_tot[ik].create(nlocal, nlocal);
                p_ham->updateHk(ik);
                const TK* hk_ptr = p_ham->getHk();
                for (int i = 0; i < ParaV->nloc; i++)
                {
                    h_mat[ik][i] = hk_ptr[i];
                }
            }

            DeePKS_domain::collect_h_mat<TK, TH>(*ParaV, h_mat, h_tot, nlocal, nks);

            const std::string file_htot = PARAM.globalv.global_out_dir + "deepks_htot.npy";
            LCAO_deepks_io::save_npy_h<TK, TH>(h_tot, file_htot, nlocal, nks, my_rank);

            if (PARAM.inp.deepks_scf)
            {
                std::vector<TH> v_delta(nks);
                std::vector<TH> h_base(nks);
                for (int ik = 0; ik < nks; ik++)
                {
                    v_delta[ik].create(nlocal, nlocal);
                    h_base[ik].create(nlocal, nlocal);
                }
                DeePKS_domain::collect_h_mat<TK, TH>(*ParaV, *h_delta, v_delta, nlocal, nks);

                // save v_delta and h_base
                const std::string file_hbase = PARAM.globalv.global_out_dir + "deepks_hbase.npy";
                for (int ik = 0; ik < nks; ik++)
                {
                    h_base[ik] = h_tot[ik] - v_delta[ik];
                }
                LCAO_deepks_io::save_npy_h<TK, TH>(h_base, file_hbase, nlocal, nks, my_rank);

                const std::string file_vdelta = PARAM.globalv.global_out_dir + "deepks_vdelta.npy";
                LCAO_deepks_io::save_npy_h<TK, TH>(v_delta, file_vdelta, nlocal, nks, my_rank);

                if (PARAM.inp.deepks_v_delta == 1) // v_delta_precalc storage method 1
                {
                    std::vector<torch::Tensor> gevdm;
                    ld->cal_gevdm(nat, gevdm);

                    torch::Tensor v_delta_precalc;
                    DeePKS_domain::cal_v_delta_precalc<TK>(nlocal,
                                                           ld->lmaxd,
                                                           ld->inlmax,
                                                           nat,
                                                           nks,
                                                           ld->inl_l,
                                                           kvec_d,
                                                           ld->phialpha,
                                                           gevdm,
                                                           ld->inl_index,
                                                           ucell,
                                                           orb,
                                                           *ParaV,
                                                           GridD,
                                                           v_delta_precalc);

                    LCAO_deepks_io::save_npy_v_delta_precalc<TK>(nat,
                                                                 nks,
                                                                 nlocal,
                                                                 ld->des_per_atom,
                                                                 v_delta_precalc,
                                                                 PARAM.globalv.global_out_dir,
                                                                 my_rank);
                }
                else if (PARAM.inp.deepks_v_delta == 2) // v_delta_precalc storage method 2
                {
                    torch::Tensor phialpha_out;
                    DeePKS_domain::prepare_phialpha<TK>(nlocal,
                                                        ld->lmaxd,
                                                        ld->inlmax,
                                                        nat,
                                                        nks,
                                                        kvec_d,
                                                        ld->phialpha,
                                                        ucell,
                                                        orb,
                                                        *ParaV,
                                                        GridD,
                                                        phialpha_out);

                    LCAO_deepks_io::save_npy_phialpha<TK>(nat,
                                                          nks,
                                                          nlocal,
                                                          ld->inlmax,
                                                          ld->lmaxd,
                                                          phialpha_out,
                                                          PARAM.globalv.global_out_dir,
                                                          my_rank);
                    std::vector<torch::Tensor> gevdm;
                    ld->cal_gevdm(nat, gevdm);

                    torch::Tensor gevdm_out;
                    DeePKS_domain::prepare_gevdm(nat, ld->lmaxd, ld->inlmax, orb, gevdm, gevdm_out);

                    LCAO_deepks_io::save_npy_gevdm(nat,
                                                   ld->inlmax,
                                                   ld->lmaxd,
                                                   gevdm_out,
                                                   PARAM.globalv.global_out_dir,
                                                   my_rank);
                }
            }
            else // deepks_scf == 0
            {
                const std::string file_hbase = PARAM.globalv.global_out_dir + "deepks_hbase.npy";
                LCAO_deepks_io::save_npy_h<TK, TH>(h_tot, file_hbase, nlocal, nks, my_rank);
            }
        } // end v_delta label

    } // end deepks_out_labels

    // DeePKS PDM and descriptor
    if (PARAM.inp.deepks_out_labels || PARAM.inp.deepks_scf)
    {
        // this part is for integrated test of deepks
        // so it is printed no matter even if deepks_out_labels is not used
        // when deepks_scf is on, the init pdm should be same as the out pdm, so we should not recalculate the pdm
        if (!PARAM.inp.deepks_scf)
        {
            ld->cal_projected_DM<TK>(dm, ucell, orb, GridD);
        }

        ld->check_projected_dm(); // print out the projected dm for NSCF calculaiton

        std::vector<torch::Tensor> descriptor;
        DeePKS_domain::cal_descriptor(nat,
                                      ld->inlmax,
                                      ld->inl_l,
                                      ld->pdm,
                                      descriptor,
                                      ld->des_per_atom); // final descriptor
        DeePKS_domain::check_descriptor(ld->inlmax,
                                        ld->des_per_atom,
                                        ld->inl_l,
                                        ucell,
                                        PARAM.globalv.global_out_dir,
                                        descriptor);

        if (PARAM.inp.deepks_out_labels)
        {
            LCAO_deepks_io::save_npy_d(nat,
                                       ld->des_per_atom,
                                       ld->inlmax,
                                       ld->inl_l,
                                       PARAM.inp.deepks_equiv,
                                       descriptor,
                                       PARAM.globalv.global_out_dir,
                                       GlobalV::MY_RANK); // libnpy needed
        }
    }

    /// print out deepks information to the screen
    if (PARAM.inp.deepks_scf)
    {
        ld->cal_e_delta_band(dm->get_DMK_vector(), nks);
        std::cout << "E_delta_band = " << std::setprecision(8) << ld->e_delta_band << " Ry"
                  << " = " << std::setprecision(8) << ld->e_delta_band * ModuleBase::Ry_to_eV << " eV" << std::endl;
        std::cout << "E_delta_NN = " << std::setprecision(8) << ld->E_delta << " Ry"
                  << " = " << std::setprecision(8) << ld->E_delta * ModuleBase::Ry_to_eV << " eV" << std::endl;
    }
}

template class LCAO_Deepks_Interface<double, double>;
template class LCAO_Deepks_Interface<std::complex<double>, double>;
template class LCAO_Deepks_Interface<std::complex<double>, std::complex<double>>;

#endif
