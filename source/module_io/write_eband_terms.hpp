#ifndef _WRITE_EBAND_SEPARATE_TERMS_HPP_
#define _WRITE_EBAND_SEPARATE_TERMS_HPP_
#include "module_io/write_vxc.hpp"
#include "module_hamilt_lcao/hamilt_lcaodft/operator_lcao/ekinetic_new.h"
#include "module_hamilt_lcao/hamilt_lcaodft/operator_lcao/nonlocal_new.h"
#include "module_basis/module_nao/two_center_bundle.h"
namespace ModuleIO
{
    template <typename TK, typename TR>
    void write_eband_terms(const int nspin,
        const int nbasis,
        const int drank,
        const Parallel_Orbitals* pv,
        const psi::Psi<TK>& psi,
        const UnitCell& ucell,
        Structure_Factor& sf,
        const ModulePW::PW_Basis& rho_basis,
        const ModulePW::PW_Basis& rhod_basis,
        const ModuleBase::matrix& vloc,
        const Charge& chg,
        Gint_Gamma& gint_gamma, // mohan add 2024-04-01
        Gint_k& gint_k,         // mohan add 2024-04-01
        const K_Vectors& kv,
        const ModuleBase::matrix& wg,
        Grid_Driver& gd,
        const std::vector<double>& orb_cutoff,
        const TwoCenterBundle& two_center_bundle
#ifdef __EXX
        , std::vector<std::map<int, std::map<TAC, RI::Tensor<double>>>>* Hexxd = nullptr
        , std::vector<std::map<int, std::map<TAC, RI::Tensor<std::complex<double>>>>>* Hexxc = nullptr
#endif
    )
    {
        // 0. prepare
        const int& nbands = wg.nc;
        const int& nspin0 = (nspin == 2) ? 2 : 1;
        double etxc = 0.0;
        double vtxc = 0.0;
        Parallel_2D p2d;
        set_para2d_MO(*pv, nbands, p2d);
        typename TGint<TK>::type* gint = nullptr;
        set_gint_pointer<TK>(gint_gamma, gint_k, gint);
        auto if_gamma_fix = [](hamilt::HContainer<TR>& hR) {if (std::is_same<TK, double>::value) { hR.fix_gamma(); }};
        auto all_band_energy = [&wg](const int ik, const std::vector<double>& e_orb)->double
            {
                double e = 0;
                for (int i = 0; i < e_orb.size(); ++i) { e += e_orb[i] * wg(ik, i); }
                return e;
            };
        auto all_k_all_band_energy = [&wg, &all_band_energy](const std::vector<std::vector<double>>& e_orb)->double
            {
                double e = 0;
                for (int ik = 0; ik < e_orb.size(); ++ik) { e += all_band_energy(ik, e_orb[ik]); }
                return e;
            };

        // 1. kinetic
        if (PARAM.inp.t_in_h)
        {
            hamilt::HS_Matrix_K<TK> kinetic_k_ao(pv, 1);
            hamilt::HContainer<TR> kinetic_R_ao(pv);
            if_gamma_fix(kinetic_R_ao);
            hamilt::EkineticNew<hamilt::OperatorLCAO<TK, TR>> kinetic_op(&kinetic_k_ao, kv.kvec_d,
                &kinetic_R_ao, &ucell, orb_cutoff, &gd, two_center_bundle.kinetic_orb.get());
            kinetic_op.contributeHR();
            std::vector<std::vector<double>> e_orb_kinetic;
            for (int ik = 0;ik < kv.get_nks();++ik)
            {
                kinetic_k_ao.set_zero_hk();
                dynamic_cast<hamilt::OperatorLCAO<TK, TR>*>(&kinetic_op)->contributeHk(ik);
                e_orb_kinetic.emplace_back(orbital_energy(ik, nbands,
                    cVc(kinetic_k_ao.get_hk(), &psi(ik, 0, 0), nbasis, nbands, *pv, p2d), p2d));
            }
            write_orb_energy(kv, nspin0, nbands, e_orb_kinetic, "kinetic", "");
            // ======test=======
            // std::cout << "e_kinetic:" << all_k_all_band_energy(e_orb_kinetic) << std::endl;
            // ======test=======
        }

        // 2. pp: local
        if (PARAM.inp.vl_in_h)
        {
            elecstate::Potential pot_local(&rhod_basis, &rho_basis, &ucell, &vloc, &sf, &etxc, &vtxc);
            pot_local.pot_register({ "local" });
            pot_local.update_from_charge(&chg, &ucell);
            hamilt::HS_Matrix_K<TK> v_pp_local_k_ao(pv, 1);
            hamilt::HContainer<TR> v_pp_local_R_ao(pv);
            if_gamma_fix(v_pp_local_R_ao);
            std::vector<std::vector<double>> e_orb_pp_local;
            hamilt::Veff<hamilt::OperatorLCAO<TK, TR>> v_pp_local_op(gint, &v_pp_local_k_ao, kv.kvec_d, &pot_local, &v_pp_local_R_ao, &ucell, orb_cutoff, &gd, nspin);
            v_pp_local_op.contributeHR();
            for (int ik = 0;ik < kv.get_nks();++ik)
            {
                v_pp_local_k_ao.set_zero_hk();
                dynamic_cast<hamilt::OperatorLCAO<TK, TR>*>(&v_pp_local_op)->contributeHk(ik);
                e_orb_pp_local.emplace_back(orbital_energy(ik, nbands,
                    cVc(v_pp_local_k_ao.get_hk(), &psi(ik, 0, 0), nbasis, nbands, *pv, p2d), p2d));
            }
            write_orb_energy(kv, nspin0, nbands, e_orb_pp_local, "vpp_local", "");
            // ======test=======
            // std::cout << "e_pp_local:" << all_k_all_band_energy(e_orb_pp_local) << std::endl;
            // ======test=======
        }

        // 3. pp: nonlocal
        if (PARAM.inp.vnl_in_h)
        {
            hamilt::HS_Matrix_K<TK> v_pp_nonlocal_k_ao(pv, 1);
            hamilt::HContainer<TR> v_pp_nonlocal_R_ao(pv);
            if_gamma_fix(v_pp_nonlocal_R_ao);
            std::vector<std::vector<double>> e_orb_pp_nonlocal;
            hamilt::NonlocalNew<hamilt::OperatorLCAO<TK, TR>> v_pp_nonlocal_op(&v_pp_nonlocal_k_ao, kv.kvec_d,
                &v_pp_nonlocal_R_ao, &ucell, orb_cutoff, &gd, two_center_bundle.overlap_orb_beta.get());
            v_pp_nonlocal_op.contributeHR();
            for (int ik = 0;ik < kv.get_nks();++ik)
            {
                v_pp_nonlocal_k_ao.set_zero_hk();
                dynamic_cast<hamilt::OperatorLCAO<TK, TR>*>(&v_pp_nonlocal_op)->contributeHk(ik);
                e_orb_pp_nonlocal.emplace_back(orbital_energy(ik, nbands,
                    cVc(v_pp_nonlocal_k_ao.get_hk(), &psi(ik, 0, 0), nbasis, nbands, *pv, p2d), p2d));
            }
            write_orb_energy(kv, nspin0, nbands, e_orb_pp_nonlocal, "vpp_nonlocal", "");
            // ======test=======
            // std::cout << "e_pp_nonlocal:" << all_k_all_band_energy(e_orb_pp_nonlocal) << std::endl;
            // ======test=======
        }

        // 4. hartree
        if (PARAM.inp.vh_in_h)
        {
            elecstate::Potential pot_hartree(&rhod_basis, &rho_basis, &ucell, &vloc, &sf, &etxc, &vtxc);
            pot_hartree.pot_register({ "hartree" });
            pot_hartree.update_from_charge(&chg, &ucell);
            std::vector<hamilt::HContainer<TR>> v_hartree_R_ao(nspin0, hamilt::HContainer<TR>(pv));
            for (int is = 0; is < nspin0; ++is)
            {
                v_hartree_R_ao[is].set_zero();
                if_gamma_fix(v_hartree_R_ao[is]);
            }
            hamilt::HS_Matrix_K<TK> v_hartree_k_ao(pv, 1);
            std::vector<hamilt::Veff<hamilt::OperatorLCAO<TK, TR>>*> v_hartree_op(nspin0);
            for (int is = 0; is < nspin0; ++is)
            {
                v_hartree_op[is] = new hamilt::Veff<hamilt::OperatorLCAO<TK, TR>>(gint,
                    &v_hartree_k_ao, kv.kvec_d, &pot_hartree, &v_hartree_R_ao[is], &ucell, orb_cutoff, &gd, nspin);
                v_hartree_op[is]->contributeHR();
            }
            std::vector<std::vector<double>> e_orb_hartree;
            for (int ik = 0;ik < kv.get_nks();++ik)
            {
                int is = kv.isk[ik];
                v_hartree_k_ao.set_zero_hk();
                dynamic_cast<hamilt::OperatorLCAO<TK, TR>*>(v_hartree_op[is])->contributeHk(ik);
                e_orb_hartree.emplace_back(orbital_energy(ik, nbands,
                    cVc(v_hartree_k_ao.get_hk(), &psi(ik, 0, 0), nbasis, nbands, *pv, p2d), p2d));
            }
            for (auto& op : v_hartree_op) { delete op; }
            write_orb_energy(kv, nspin0, nbands, e_orb_hartree, "vhartree", "");
            // ======test=======
            // std::cout << "e_hartree:" << all_k_all_band_energy(e_orb_hartree) << std::endl;
            // ======test=======
        }

        // 5. xc (including exx)
        if (!PARAM.inp.out_mat_xc)  // avoid duplicate output
        {
            write_Vxc<TK, TR>(nspin, nbasis, drank, pv, psi, ucell, sf, rho_basis, rhod_basis, vloc, chg, gint_gamma, gint_k, kv, orb_cutoff, wg, gd
#ifdef __EXX
                , Hexxd, Hexxc
#endif
            );
        }
    }
}
#endif
