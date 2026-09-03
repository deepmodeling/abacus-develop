#include "lr_force.h"
#include "cal_hs_grad.h"
#include "pulay_force_hcontainer.h"
#include "source_lcao/pulay_fs.h"   // only for gint terms
#include "source_lcao/module_lr/utils/lr_util_hcontainer.h"
#include "source_lcao/module_lr/utils/lr_util_print.h"
#ifdef __EXX
#include "source_lcao/module_lr/operator_casida/operator_lr_exx.h"
#endif
namespace LR
{
    template<typename TK>
    ModuleBase::matrix LR_Force<TK>::cal_force_overlap_edm(const elecstate::DensityMatrix<TK, double>& edm)
    {
        // const double* dS[3] = { dSloc_x,  dSloc_y,  dSloc_z };
        std::vector<hamilt::HContainer<double>>  dS = cal_hs_grad('S', this->ucell_, this->pv_, this->gd_, this->two_center_bundle_);
        // test: output dS
        // std::cout << "dS in 3 directions:\n";
        // for (int i = 0;i < 3;++i) { LR_Util::print_HR(dS.at(i), this->ucell_.nat, "dS" + std::to_string(i)); }
        ModuleBase::matrix foverlap = PulayForceStress::cal_pulay_fs(edm, this->ucell_, dS, -1.);
        if (PARAM.inp.test_force)
        {
            ModuleIO::print_force(GlobalV::ofs_running, this->ucell_, "OVERLAP     FORCE (eV/Angstrom)", foverlap, false);
        }
        return foverlap;
    }

    template<typename TK>
    ModuleBase::matrix LR_Force<TK>::reproduce_force_gs(const K_Vectors& kv,
        const elecstate::DensityMatrix<TK, double>& dm_gs,
        const elecstate::DensityMatrix<TK, double>& edm_gs)
    {
        const int& nspin = PARAM.inp.nspin;
        // local + Hartree + xc term, including Hellmann-Feynman and Pulay
        ModuleBase::matrix f_gs_hf_pulay = cal_force_hamilt_gs_dm_relaxed_diff(dm_gs, dm_gs, true); // pw(vl_dvl+ewald)+vnl+t_dphi+vl_dphi
        // edm term
        ModuleBase::matrix f_nonortho = cal_force_overlap_edm(edm_gs); // overlap
#ifdef __EXX
        if (gs_is_hybrid())
        {
            const auto& Ds_gs = LR_Util::get_exx_Ds_gs(dm_gs, ucell_, kv, pv_);
            const auto& Ds_gs_2 = LR_Util::get_exx_Ds_gs(dm_gs, ucell_, kv, pv_);
            // at nspin=1 there is only one channel, and `get_exx_Ds_gs` already returns 0.5*D
            // (= D_up), so both halves below reuse channel 0.
            const int is_2nd = (Ds_gs.size() > 1) ? 1 : 0;
            ModuleBase::matrix f_gs_exx(ucell_.nat, 3);
            // test the two function using the two spin channels respectively
            // 0.5 is from dE = 0.5 dTr[D(HD)]. No 0.5 in excited-state calculateion of dTr[(T+Z)(HD)]
            f_gs_exx += cal_force_exx_gs_dm_relaxed_diff(Ds_gs.at(0), Ds_gs_2.at(0), alpha_, std::to_string(0)) * 0.5;    // test passed， = 0.5 groud-state EXX force
            f_gs_exx += cal_force_exx_dm_trans(Ds_gs.at(is_2nd), alpha_, std::to_string(is_2nd)) * 0.5;
            if (PARAM.inp.test_force)
                ModuleIO::print_force(GlobalV::ofs_running, ucell_, "EXX GS FORCE reproduce (eV/Angstrom)", f_gs_exx, false);
            f_gs_hf_pulay += f_gs_exx;
        }
#endif
        return f_gs_hf_pulay + f_nonortho;
    }

    template<typename TK>
    ModuleBase::matrix LR_Force<TK>::reproduce_force_gs_loc(
        const elecstate::DensityMatrix<TK, double>& dm_gs,
        const elecstate::Potential& pot_gs)
    {
        const Charge chr_gs = dm_to_charge(dm_gs);
        //  local pp (Pulay) + Hartree + xc (grid integration)
        ModuleBase::matrix fvl_dphi(this->ucell_.nat, 3);
        ModuleBase::matrix stress_tmp;  // no use now, only for passing into interfaces
        PulayForceStress::cal_pulay_fs(dm_gs.get_DMR_vector().size()/*nspin*/, fvl_dphi, stress_tmp,
            dm_gs, this->ucell_, &pot_gs, true, false);
        return fvl_dphi;
    }

    template<typename TK>
    void LR_Force<TK>::cal_H2_sz_center2_deriv(const std::vector<double>& orb_cutoffs, const K_Vectors& kv)
    {
        GlobalV::ofs_running << "  ==== Test H2_SZ_CENTER2_DERIV dtau(Sij) and dtau(hij) ====" << std::endl;
        const std::vector<ModuleBase::Vector3<double>>& kvd_test = { ModuleBase::Vector3<double>(0.0, 0.0, 0.0) };
        auto init_dm_eff = [&, this](const int i, const int j) -> elecstate::DensityMatrix<TK, double>
            {   // dm_{ij}=1, other elements = 0, i,j = 0,1
                std::vector<TK> dm_2d(4, 0.0);
                std::cout << "i<<1 + j =" << ((i << 1) + j) << std::endl;
                dm_2d[i * 2 + j] = 1.0;
                LR_Util::matsym(dm_2d.data(), 2);   //symmetrization is a must for calling 1-electron Pulay force funcs
                elecstate::DensityMatrix<TK, double> dm(&this->pv_, 1, kvd_test, 1);
                dm.set_DMK_pointer(0, dm_2d.data());
                LR_Util::initialize_DMR(dm, this->pv_, this->ucell_, this->gd_, orb_cutoffs);
                dm.cal_DMR();
                return dm;
            };
        for (auto&& i : { 0, 1 })
            for (auto&& j : { 0, 1 })
            {
                elecstate::DensityMatrix<TK, double> dm_ij = init_dm_eff(i, j);
                elecstate::Potential pot_hij = dm_to_hxc_potential(dm_ij);
                // 1. dtau(S_ij)
                {
                    std::vector<hamilt::HContainer<double>>  dS = cal_hs_grad('S', this->ucell_, this->pv_, this->gd_, this->two_center_bundle_);   // (dr i|j)
                    ModuleBase::matrix foverlap = PulayForceStress::cal_pulay_fs(dm_ij, this->ucell_, dS, -1.); //dtau(i|j), related to (dr i|j) (1, -1 or 2)
                    ModuleIO::print_force(GlobalV::ofs_running, this->ucell_,
                        "H2_SZ_CENTER2_dtau_S(" + std::to_string(i) + std::to_string(j) + ") FORCE (Ry/au)",
                        foverlap, true);   // F_S_ij = dtau(S_ij)
                }
                // 2. dtau(h_ij), h = T + Vl + Vnl
                {
                    ModuleBase::matrix stress_tmp;  // dummy
                    // kinetic (Pulay term only)
                    std::vector<hamilt::HContainer<double>> dT = cal_hs_grad('T', this->ucell_, this->pv_, this->gd_, this->two_center_bundle_);
                    ModuleBase::matrix ft_dphi = PulayForceStress::cal_pulay_fs(dm_ij, this->ucell_, dT);

                    // local pp Hellmann-Feynman term (which does not depend on the charge density if Hxc is not included)
                    const Charge chr_dummy = dm_to_charge(dm_ij);
                    ModuleBase::matrix fvl_dvl = PARAM.inp.vl_in_h ?
                        ForcePWTerms<double>()(this->ucell_, chr_dummy, this->rhopw_, this->locpp_, this->sf_, /*with_ewald=*/ false) :
                        ModuleBase::matrix(this->ucell_.nat, 3);
                    // local pp Pulay term
                    ModuleBase::matrix fvl_dphi(this->ucell_.nat, 3);
                    elecstate::Potential pot_loc = this->local_potential();
                    PulayForceStress::cal_pulay_fs(dm_ij.get_DMR_vector().size()/*nspin*/, fvl_dphi, stress_tmp,
                        dm_ij, this->ucell_, &pot_loc, true, false);

                    // nonlocal pp term (Hellmann-Feynman + Pulay)
                    ModuleBase::matrix fvnl = cal_force_nonlocal(this->ucell_, this->kvec_d_, this->gd_, this->two_center_bundle_, dm_ij);

                    ModuleIO::print_force(GlobalV::ofs_running, this->ucell_,
                        "H2_SZ_CENTER2_dtau_h1e(" + std::to_string(i) + std::to_string(j) + ") FORCE (Ry/au)",
                        (ft_dphi + fvl_dvl + fvl_dphi + fvnl) * (-1), true);   // F_h_ij = -dtau(h_ij)
                }
            }
    }

    template<typename TK>
    void LR_Force<TK>::cal_H2_sz_center4(const std::vector<double>& orb_cutoffs,
        const K_Vectors& kv, const bool is_grad)
    {
        const std::string label = is_grad ? "dtau(ij | kl)" : "(ij | kl)";;
        GlobalV::ofs_running << "  ==== Test H2_SZ_CENTER4_HXC " << label << " ====" << std::endl;
        const std::vector<ModuleBase::Vector3<double>>& kvd_test = { ModuleBase::Vector3<double>(0.0, 0.0, 0.0) };
        auto init_dm_eff = [&, this](const int i, const int j, const bool symmetrize = false) -> elecstate::DensityMatrix<TK, double>
            {   // dm_{ij}=1, other elements = 0, i,j = 0,1
                std::vector<TK> dm_2d(4, 0.0);
                std::cout<<"i<<1 + j =" << ((i<<1) + j) << std::endl;
                dm_2d[i * 2 + j] = 1.0;
                if (symmetrize) { LR_Util::matsym(dm_2d.data(), 2); }   //symmetrization is a must for calling 1-electron Pulay force funcs
                elecstate::DensityMatrix<TK, double> dm(&this->pv_, 1, kvd_test, 1);
                dm.set_DMK_pointer(0, dm_2d.data());
                LR_Util::initialize_DMR(dm, this->pv_, this->ucell_, this->gd_, orb_cutoffs);
                dm.cal_DMR();
                return dm;
            };
#ifdef __EXX
        std::vector<std::tuple<std::set<int>, std::set<int>>> judge = RI_2D_Comm::get_2D_judge(ucell_, pv_);
#endif

        for (auto&& i : { 0, 1 })
            for (auto&& j : { 0, 1 })
            {
                elecstate::DensityMatrix<TK, double> dm_ij = init_dm_eff(i, j, false);
                elecstate::Potential pot_hxc_ij = dm_to_hxc_potential(dm_ij);
                elecstate::DensityMatrix<TK, double> dm_ij_sym = init_dm_eff(i, j, true);
                for (auto&& k : { 0, 1 })
                    for (auto&& l : { 0, 1 })
                    {
                        // 1. build dm(kl)
                        elecstate::DensityMatrix<TK, double> dm_kl = init_dm_eff(k, l, false);
                        if (is_grad)
                        {
                            // 2. pulay term + Hellmann-Feynman term
                            elecstate::Potential pot_hxc_kl = dm_to_hxc_potential(dm_kl);
                            elecstate::DensityMatrix<TK, double> dm_kl_sym = init_dm_eff(k, l, true);
                            ModuleBase::matrix fhartree_pulay(this->ucell_.nat, 3), fhartree_h_f(this->ucell_.nat, 3);
                            ModuleBase::matrix stress_tmp;  // dummy
                            PulayForceStress::cal_pulay_fs(1/*nspin*/, fhartree_pulay, stress_tmp, dm_ij_sym, this->ucell_, &pot_hxc_kl, true, false);  // Pulay term
                            PulayForceStress::cal_pulay_fs(1/*nspin*/, fhartree_h_f, stress_tmp, dm_kl_sym, this->ucell_, &pot_hxc_ij, true, false);  // Hellmann-Feynman term
                            ModuleIO::print_force(GlobalV::ofs_running, this->ucell_,
                                "H2_SZ_CENTER4_HXC_dtau(" + std::to_string(i) + std::to_string(j) + "|" + std::to_string(k) + std::to_string(l) + ") FORCE (Ry/au)",
                                -(fhartree_pulay + fhartree_h_f), true);   // F_Hxc_ijkl = -dtau(ij|kl)
                            ModuleIO::print_force(GlobalV::ofs_running, this->ucell_,
                                "H2_SZ_CENTER4_HXC_Pulay_dtau(" + std::to_string(i) + std::to_string(j) + "|" + std::to_string(k) + std::to_string(l) + ") FORCE (Ry/au)",
                                -fhartree_pulay, true);   // F_Hxc_ijkl = -dtau(ij|kl)
                            ModuleIO::print_force(GlobalV::ofs_running, this->ucell_,
                                "H2_SZ_CENTER4_HXC_H-F_dtau(" + std::to_string(i) + std::to_string(j) + "|" + std::to_string(k) + std::to_string(l) + ") FORCE (Ry/au)",
                                -fhartree_h_f, true);   // F_Hxc_ijkl = -dtau(ij|kl)
#ifdef __EXX
                            if (!this->exx_lri_.expired())
                            {   // match the Gint result with LibRI
                                auto ds_kl = LR_Util::get_exx_Ds_spin1(dm_kl, ucell_, kv, pv_); // returns ds_kl*0.5
                                auto ds_ij = LR_Util::get_exx_Ds_spin1(dm_ij, ucell_, kv, pv_); // returns ds_ij*0.5
                                // calulates F=0.5*d(ik|jl) (only one spin channel). 0.5 is the 2-electron integral prefactor.
                                // 4 cancels the two 0.5s in Ds, induced by `split_m2D_ktoR`. 
                                // No spin factor or two-electron-energy factor (1/2) are hard-coded in this function.
                                ModuleBase::matrix f_exx = this->cal_force_exx_gs_dm_relaxed_diff(ds_kl, ds_ij, alpha_ * 4.0, "");
                                ModuleIO::print_force(GlobalV::ofs_running, ucell_,
                                    "H2_SZ_CENTER4_EXX_dtau(" + std::to_string(i) + std::to_string(j) + "|" + std::to_string(k) + std::to_string(l) + ") FORCE (Ry/au)",
                                    f_exx , true);   //d(ik|jl)=F
                            }
#endif
                        }
                        else
                        {
                            //2. build charge & potential
                            elecstate::Potential pot_hxc_kl = dm_to_hxc_potential(dm_kl);
                            const Charge& charge_ij = this->dm_to_charge(dm_ij);
                            // 3. cal energy
                            double e_hxc = std::inner_product(charge_ij.rho[0],
                                charge_ij.rho[0] + this->rhopw_.nrxx,
                                pot_hxc_kl.get_eff_v(0), 0.0) * 0.5 * this->ucell_.omega / static_cast<double>(this->rhopw_.nrxx);
                            GlobalV::ofs_running << "  H2_SZ_CENTER4_COULOMB ("
                                << std::to_string(i) + std::to_string(j) + "|" + std::to_string(k) + std::to_string(l)
                                << ") by Gint: " << std::setprecision(15) << e_hxc * 2 << std::endl; // 2 for testing (ij|kl) instead of real Coulomb energy 0.5*(ij|kl)
#ifdef __EXX
                            if (!this->exx_lri_.expired())
                            {   // match the Gint result with LibRI
                                auto ds_kl = LR_Util::get_exx_Ds_spin1(dm_kl, ucell_, kv, pv_); // returns ds_kl*0.5
                                auto ds_ij = LR_Util::get_exx_Ds_spin1(dm_ij, ucell_, kv, pv_); // returns ds_ij*0.5
                                auto lri = this->exx_lri_.lock();
                                lri->get().set_Ds(std::move(ds_kl), lri->get_info().dm_threshold);
                                lri->get().cal_Hs();
                                lri->Hexxs[0] = RI::Communicate_Tensors_Map_Judge::comm_map2_first(
                                    lri->get_mpi_comm(), std::move(lri->get().Hs), std::get<0>(judge[0]), std::get<1>(judge[0]));
                                lri->post_process_Hexx(lri->Hexxs[0]);
                                TK e_exx = this->alpha_ * lri->get().post_2D.cal_energy(ds_ij, lri->Hexxs[0]) * 2.0; // 4 is to cancel two 0.5^2 in split_m2D_ktoR(nspin=1)`, and 0.5 for Fock energy
                                GlobalV::ofs_running << "  H2_SZ_CENTER4_COULOMB ("
                                    << std::to_string(i) + std::to_string(k) + "|" + std::to_string(j) + std::to_string(l)
                                    << ") by LibRI: " << std::setprecision(15) << -e_exx * 2.0 << " where alpha = " << this->alpha_ << std::endl;  //-2 for Fock energy -> integral
                            }
#endif  
                        }

                    }
            }
    }
}



template class LR::LR_Force<double>;
template class LR::LR_Force<std::complex<double>>;