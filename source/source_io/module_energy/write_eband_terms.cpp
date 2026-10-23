#include "write_eband_terms.h"

#include "source_io/module_hs/vxc_op_tools.h"
#include "source_io/module_parameter/parameter.h"
#include "source_lcao/module_operator_lcao/ekinetic.h"
#include "source_lcao/module_operator_lcao/nonlocal.h"
#include "source_lcao/module_operator_lcao/veff_lcao.h"
#include "source_cell/module_neighbor/sltk_grid_driver.h"

#include "source_base/constants.h"
#include "source_psi/psi.h"
#include "source_cell/unitcell.h"
#include "source_cell/klist.h"
#include "source_estate/module_pot/potential_new.h"
#include "source_basis/module_nao/two_center_bundle.h"
#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_estate/module_charge/charge.h"
#include "source_pw/module_pwdft/stru_fac.h"
#include "source_hamilt/module_surchem/surchem.h"
#include "source_basis/module_pw/pw_basis.h"

#include <string>

namespace ModuleIO
{
template <typename TK, typename TR>
void write_eband_terms(const int nspin,
                       const int istep,
                       const double etot,
                       const int my_rank,
                       const int nbasis,
                       const Parallel_Orbitals* pv,
                       const psi::Psi<TK>& psi,
                       const UnitCell& ucell,
                       Structure_Factor& sf,
                       surchem& solvent,
                       const ModulePW::PW_Basis& rho_basis,
                       const ModulePW::PW_Basis& rhod_basis,
                       const ModuleBase::matrix& vloc,
                       const Charge& chg,
                       const K_Vectors& kv,
                       const ModuleBase::matrix& wg,
                       Grid_Driver& gd,
                       const std::vector<double>& orb_cutoff,
                       const TwoCenterBundle& two_center_bundle
)
    {
        // 0. prepare
        const int& nbands = wg.nc;
        const int nspin0 = (nspin == 2) ? 2 : 1;
        const double etot_eV = etot * ModuleBase::Ry_to_eV;
        const std::string time_stamp = local_time_stamp();
        double etxc = 0.0;
        double vtxc = 0.0;

        Parallel_2D p2d;

        set_para2d_MO(*pv, nbands, p2d);

        auto if_gamma_fix = [](hamilt::HContainer<TR>& hR)
        {
            if (std::is_same<TK, double>::value)
            {
                hR.fix_gamma();
            }
        };

        // Only rank 0 writes; the orbital energies above are reduced on all ranks.
        auto write_term = [&](const std::vector<std::vector<double>>& e_orb,
                              const std::string& out_name,
                              const std::string& term_symbol)
        {
            if (my_rank != 0)
            {
                return;
            }
            const std::string header = build_eband_header(istep, etot_eV, term_symbol, time_stamp);
            const std::string& out_dir = PARAM.globalv.global_out_dir;
            write_orb_energy(kv, nspin0, nbands, e_orb, out_name, header, out_dir);
        };

        auto all_band_energy = [&wg](const int ik, const std::vector<double>& e_orb)->double
        {
            double e = 0;
            for (int i = 0; i < e_orb.size(); ++i) { e += e_orb[i] * wg(ik, i); }
            return e;
        };

        auto all_k_all_band_energy = [&wg, &all_band_energy](const std::vector<std::vector<double>>& e_orb)->double
        {
            double e = 0;
            for (int ik = 0; ik < e_orb.size(); ++ik)
            {
                e += all_band_energy(ik, e_orb[ik]);
            }
            return e;
        };

        // 1. kinetic
        if (PARAM.inp.t_in_h)
        {
            hamilt::HS_Matrix_K<TK> kinetic_k_ao(pv, 1);
            hamilt::HContainer<TR> kinetic_R_ao(pv);
            if_gamma_fix(kinetic_R_ao);

            hamilt::EKinetic<hamilt::OperatorLCAO<TK, TR>> kinetic_op(&kinetic_k_ao, kv.kvec_d,
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

            write_term(e_orb_kinetic, "eband_kin.txt", "T");
        }

        // 2. pp: local
        if (PARAM.inp.vl_in_h)
        {
            elecstate::Potential pot_local(&rhod_basis, &rho_basis, &ucell, &vloc, &sf, &solvent, &etxc, &vtxc);
            pot_local.pot_register({ "local" });
            pot_local.update_from_charge(&chg, &ucell);
            hamilt::HS_Matrix_K<TK> v_pp_local_k_ao(pv, 1);
            hamilt::HContainer<TR> v_pp_local_R_ao(pv);
            if_gamma_fix(v_pp_local_R_ao);
            std::vector<std::vector<double>> e_orb_pp_local;

            hamilt::Veff<hamilt::OperatorLCAO<TK, TR>> v_pp_local_op(
                    &v_pp_local_k_ao,
                    kv.kvec_d,
                    &pot_local,
                    &v_pp_local_R_ao,
                    &ucell,
                    orb_cutoff,
                    &gd,
                    nspin);

            v_pp_local_op.contributeHR();
            for (int ik = 0;ik < kv.get_nks();++ik)
            {
                v_pp_local_k_ao.set_zero_hk();
                dynamic_cast<hamilt::OperatorLCAO<TK, TR>*>(&v_pp_local_op)->contributeHk(ik);
                e_orb_pp_local.emplace_back(orbital_energy(ik, nbands,
                    cVc(v_pp_local_k_ao.get_hk(), &psi(ik, 0, 0), nbasis, nbands, *pv, p2d), p2d));
            }
            write_term(e_orb_pp_local, "eband_vl.txt", "V_loc");
        }

        // 3. pp: nonlocal
        if (PARAM.inp.vnl_in_h)
        {
            hamilt::HS_Matrix_K<TK> v_pp_nonlocal_k_ao(pv, 1);
            hamilt::HContainer<TR> v_pp_nonlocal_R_ao(pv);
            if_gamma_fix(v_pp_nonlocal_R_ao);
            std::vector<std::vector<double>> e_orb_pp_nonlocal;
            hamilt::Nonlocal<hamilt::OperatorLCAO<TK, TR>> v_pp_nonlocal_op(&v_pp_nonlocal_k_ao, kv.kvec_d,
                &v_pp_nonlocal_R_ao, &ucell, orb_cutoff, &gd, two_center_bundle.overlap_orb_beta.get());
            v_pp_nonlocal_op.contributeHR();
            for (int ik = 0;ik < kv.get_nks();++ik)
            {
                v_pp_nonlocal_k_ao.set_zero_hk();
                dynamic_cast<hamilt::OperatorLCAO<TK, TR>*>(&v_pp_nonlocal_op)->contributeHk(ik);
                e_orb_pp_nonlocal.emplace_back(orbital_energy(ik, nbands,
                    cVc(v_pp_nonlocal_k_ao.get_hk(), &psi(ik, 0, 0), nbasis, nbands, *pv, p2d), p2d));
            }
            write_term(e_orb_pp_nonlocal, "eband_vnl.txt", "V_nonlocal");
        }

        // 4. hartree
        if (PARAM.inp.vh_in_h)
        {
            elecstate::Potential pot_hartree(&rhod_basis, &rho_basis, &ucell, &vloc, &sf, &solvent, &etxc, &vtxc);
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
                v_hartree_op[is] = new hamilt::Veff<hamilt::OperatorLCAO<TK, TR>>(
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
            write_term(e_orb_hartree, "eband_hart.txt", "V_H");
        }
        // NOTE: the XC term (eband_vxc.txt etc.) is written by write_Vxc from
        // ctrl_runner_lcao.cpp for both out_mat_xc and out_eband,
        // so it is intentionally not handled here.
    }

// Explicit template instantiations
template void write_eband_terms<double, double>(
    const int, const int, const double, const int, const int, const Parallel_Orbitals*, const psi::Psi<double>&,
    const UnitCell&, Structure_Factor&, surchem&, const ModulePW::PW_Basis&, const ModulePW::PW_Basis&,
    const ModuleBase::matrix&, const Charge&, const K_Vectors&, const ModuleBase::matrix&, Grid_Driver&,
    const std::vector<double>&, const TwoCenterBundle&
);

template void write_eband_terms<std::complex<double>, double>(
    const int, const int, const double, const int, const int, const Parallel_Orbitals*, const psi::Psi<std::complex<double>>&,
    const UnitCell&, Structure_Factor&, surchem&, const ModulePW::PW_Basis&, const ModulePW::PW_Basis&,
    const ModuleBase::matrix&, const Charge&, const K_Vectors&, const ModuleBase::matrix&, Grid_Driver&,
    const std::vector<double>&, const TwoCenterBundle&
);

template void write_eband_terms<std::complex<double>, std::complex<double>>(
    const int, const int, const double, const int, const int, const Parallel_Orbitals*, const psi::Psi<std::complex<double>>&,
    const UnitCell&, Structure_Factor&, surchem&, const ModulePW::PW_Basis&, const ModulePW::PW_Basis&,
    const ModuleBase::matrix&, const Charge&, const K_Vectors&, const ModuleBase::matrix&, Grid_Driver&,
    const std::vector<double>&, const TwoCenterBundle&
);

} // namespace ModuleIO
