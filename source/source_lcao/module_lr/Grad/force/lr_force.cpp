#include "lr_force.h"
#include "cal_hs_grad.h"
#include "pulay_force_hcontainer.h"
#include "source_lcao/pulay_fs.h"   // only for gint terms
#include "source_hamilt/module_gint/gint_interface.h"
#include "source_lcao/module_lr/utils/lr_util.h"
// #include "source_lcao/module_lr/utils/lr_util_hcontainer.h"
namespace LR
{
    template<typename TK>
    Charge LR_Force<TK>::dm_to_charge(const elecstate::DensityMatrix<TK, double>& dm)
    {
        const int& nspin_dm = dm.get_DMR_vector().size();
        const int& nspin_global = PARAM.inp.nspin;
        Charge chr;
        chr.set_rhopw(const_cast<ModulePW::PW_Basis*>(&this->rhopw_));
        chr.allocate(nspin_global, /*kin_den=*/false); //chr still needs global nspin, because Forces (PW) depends on it
        // So huge a Charge class...
        // 1. Using a (private) `allocate_rho` to control whether to delete will definately cause memory leak here. No need for such judgement. 
        // 2. Charge-dependent interfaces need refactor: only rhopw_ and rho dependence are enough.

        ModuleGint::cal_gint_rho(dm.get_DMR_vector(), nspin_dm, chr.rho, false);
        // if (nspin_dm == 1 && nspin_global == 2), chr.rho[1][irxx]=0 has been set in Charge::allocate()
        return chr;
    }

    template<typename TK>
    elecstate::Potential LR_Force<TK>::dm_to_hxc_potential(const elecstate::DensityMatrix<TK, double>& dm)
    {
        double etxc = 0.0, vtxc = 0.0;
        elecstate::Potential pot(&this->rhodpw_, &this->rhopw_, &this->ucell_,
            &this->locpp_.vloc, const_cast<Structure_Factor*>(&this->sf_),
            nullptr/*surchem*/, &etxc, &vtxc);
        PARAM.inp.vh_in_h ? pot.pot_register({ "hartree", "xc" }) : pot.pot_register({ "xc" });
        const Charge& charge = this->dm_to_charge(dm);
        pot.init_pot(&charge); // call update_from_charge inside
        return pot;
    }

    template<typename TK>
    elecstate::Potential LR_Force<TK>::local_potential()
    {
        elecstate::Potential pot(&this->rhodpw_, &this->rhopw_, &this->ucell_,
            &this->locpp_.vloc, const_cast<Structure_Factor*>(&this->sf_),
            nullptr/*surchem*/, nullptr/*etxc*/, nullptr/*vtxc*/);
        pot.pot_register({ "local" });
        pot.init_pot(nullptr);
        return pot;
    }

    template<typename TK>
    ModuleBase::matrix LR_Force<TK>::cal_force_hamilt_gs_dm_relaxed_diff(const elecstate::DensityMatrix<TK, double>& relax_diff_dm,
        const elecstate::DensityMatrix<TK, double>& dm_gs,
        const bool reproduce_gs)
    {
        const bool with_ewald = reproduce_gs;
        const Charge chr_diff_relaxed = dm_to_charge(relax_diff_dm);

        // 1. local pp (Hellmann-Feynman)(fvl_dvl) + ewald + core correction (+ self-consistent charge)
        ModuleBase::matrix f_pw = PARAM.inp.vl_in_h ?
            ForcePWTerms<double>()(this->ucell_, chr_diff_relaxed, this->rhopw_, this->locpp_, this->sf_, with_ewald) :
            ModuleBase::matrix(this->ucell_.nat, 3);

        // 2. nonlocal pp (Hellmann-Feynman + Pulay)
        ModuleBase::matrix fvnl = cal_force_nonlocal(this->ucell_, this->kvec_d_, this->gd_, this->two_center_bundle_, relax_diff_dm);

        // // 3. local pp (Pulay) + Hartree + xc (grid integration)
        // ModuleBase::matrix fvl_dphi(this->ucell_.nat, 3);
        // ModuleBase::matrix stress_tmp;  // no use now, only for passing into interfaces
        // PulayForceStress::cal_pulay_fs(relax_diff_dm.get_DMR_vector().size()/*nspin*/, fvl_dphi, stress_tmp,
        //     relax_diff_dm, this->ucell_, &pot_gs, true, false);

        // 3.1. local pp (Pulay)
        ModuleBase::matrix fvl_dphi(this->ucell_.nat, 3);
        ModuleBase::matrix stress_tmp;  // no use now, only for passing into interfaces
        elecstate::Potential pot_loc = this->local_potential();
        PulayForceStress::cal_pulay_fs(relax_diff_dm.get_DMR_vector().size()/*nspin*/, fvl_dphi, stress_tmp,
            relax_diff_dm, this->ucell_, &pot_loc, true, false);

        // 3.2. Hartree + xc (Pulay) 
        //  method 1
        // ModuleBase::matrix fgs_dphi(this->ucell_.nat, 3);
        // PulayForceStress::cal_pulay_fs(relax_diff_dm.get_DMR_vector().size()/*nspin*/, fgs_dphi, stress_tmp,
        //     relax_diff_dm, this->ucell_, &pot_gs, true, false);
        // ModuleBase::matrix fhxc_dphi = (fgs_dphi - fvl_dphi) * 0.5; // avoid double count of hxc Pulay term
        // method 2 
        ModuleBase::matrix fhxc_dphi(this->ucell_.nat, 3);
        elecstate::Potential pot_hxc = this->dm_to_hxc_potential(dm_gs);
        // `cal_pulay_fs` calculates 1*Pulay-term. 
        // For ground-state DFT, Pulay term = Hellmann-Feynman term, F = 1/2(Pulay + H-F) = Pulay, so directly call it once gives correct result.
        PulayForceStress::cal_pulay_fs(relax_diff_dm.get_DMR_vector().size()/*nspin*/, fhxc_dphi, stress_tmp,
            relax_diff_dm, this->ucell_, &pot_hxc, true, false);
        if (reproduce_gs) {fhxc_dphi *= 0.5;} // avoid double count

        // 3.3 Hartree + xc (Hellmann-Feynman)
        ModuleBase::matrix fhxc_dvhxc(this->ucell_.nat, 3);
        elecstate::Potential pot_hxc_relaxed_diff = this->dm_to_hxc_potential(relax_diff_dm);
        //`cal_pulay_fs` calculates only one spin channel because `relax_diff_dm` has only one.
        PulayForceStress::cal_pulay_fs(1/*nspin*/, fhxc_dvhxc, stress_tmp,
            dm_gs, this->ucell_, &pot_hxc_relaxed_diff, true, false);
        if(!reproduce_gs) {fhxc_dvhxc *= 2;} // for the two channels of the ground-state dm. 

        // 4. kinetic (Pulay)
        std::vector<hamilt::HContainer<double>> dT = cal_hs_grad('T', this->ucell_, this->pv_, this->gd_, this->two_center_bundle_);
        ModuleBase::matrix ft_dphi = PulayForceStress::cal_pulay_fs(relax_diff_dm, this->ucell_, dT);

        if (PARAM.inp.test_force)
        {
            ModuleIO::print_force(GlobalV::ofs_running, this->ucell_, "PW      FORCE (eV/Angstrom)", f_pw, false);
            ModuleIO::print_force(GlobalV::ofs_running, this->ucell_, "NONLOCAL     FORCE (eV/Angstrom)", fvnl, false);
            ModuleIO::print_force(GlobalV::ofs_running, this->ucell_, "KINETIC     FORCE (eV/Angstrom)", ft_dphi, false);
            ModuleIO::print_force(GlobalV::ofs_running, this->ucell_, "LOCAL-PP Pulay FORCE (eV/Angstrom)", fvl_dphi, false);
            ModuleIO::print_force(GlobalV::ofs_running, this->ucell_, "HARTREE+XC Pulay FORCE (eV/Angstrom)", fhxc_dphi, false);
            ModuleIO::print_force(GlobalV::ofs_running, this->ucell_, "HARTREE+XC Hellmann-Feynman FORCE (eV/Angstrom)", fhxc_dvhxc, false);
        }

        // from the formula, we do not need the non-ortho term (overlap*edm) here.
        return f_pw + fvnl + ft_dphi + fvl_dphi + fhxc_dphi + fhxc_dvhxc;
    }
    template<typename TK>
    ModuleBase::matrix LR_Force<TK>::cal_force_hxc_dmtrans(const elecstate::DensityMatrix<TK, double>& dm_trans, const PotHxcLR& pot_hxc)
    {
        // `dm_trans` (D^X) carries the singlet spin normalization (sqrt(2) per channel), 
        // so D^X in pot_hxc and cal_pulay_fs together already contribute a factor 2.
        // So cal_pulay_fs here returns 2*Pulay = Pulay + Hellmann-Feynman force. *2 is not needed here.
        return PulayForceStress::cal_pulay_fs(dm_trans, this->ucell_, &pot_hxc);
    }

#ifdef __EXX
    template<typename TK>
    ModuleBase::matrix LR_Force<TK>::cal_force_exx_dm_trans(
        const std::map<int, std::map<TAC, RI::Tensor<TK>>>& dm_trans,
        const double& alpha,
        const std::string& spin_suffix)
    {
        ModuleBase::matrix f_exx_dmtrans(this->ucell_.nat, 3);
        auto& exx_lri_kernel = this->exx_lri_.lock()->get();
        exx_lri_kernel.set_Ds(dm_trans, this->exx_lri_.lock()->get_info().dm_threshold, spin_suffix);
        exx_lri_kernel.cal_Hs({ "", "", spin_suffix });
        exx_lri_kernel.cal_force({ "", "", spin_suffix, "", "" });// using dm_trans，Pulay term only
        for (std::size_t idim = 0; idim < 3; ++idim)
            for (const auto& force_item : exx_lri_kernel.force[idim])
                f_exx_dmtrans(force_item.first, idim) = std::real(force_item.second);
        const double fac = -2.0 * alpha; //-2 is the same as post_process_Hexx, Hartree to Ry (which didn't act on Hs)
        const double pulay_to_total_sym = 2.0;  // Pulay -> Pulay + Hellmann-Feynman, only when Ds_left and Ds_right are equal
        return f_exx_dmtrans * fac * pulay_to_total_sym; // dm_trans (DX) already contain the spin channel (sqrt(2) times of up/down channel DX)
        // return f_exx_dmtrans * fac * 2; // 2 is the same in post_process_Eexx at nspin=1 ( up->up + down->down, 2 spin-conserving transitions)
    }

    template<typename TK>
    ModuleBase::matrix LR_Force<TK>::cal_force_exx_gs_dm_relaxed_diff(
        const std::map<int, std::map<TAC, RI::Tensor<TK>>>& dm_gs,
        const std::map<int, std::map<TAC, RI::Tensor<TK>>>& relaxed_diff_dm,
        const double& alpha,
        const std::string& spin_suffix)
    {
        ModuleBase::matrix f_exx_gs_diff(this->ucell_.nat, 3);
        auto& exx_lri_kernel = this->exx_lri_.lock()->get();
        ExxForceTwoDM<int, int, 3, TK> lr_exx_kernel(std::move(exx_lri_kernel));

        auto add_force_from_kernel = [&]() {
            for (std::size_t idim = 0; idim < 3; ++idim)
                for (const auto& force_item : lr_exx_kernel.force[idim])
                    f_exx_gs_diff(force_item.first, idim) += std::real(force_item.second);
            };

        auto transpose_dm = [](const std::map<int, std::map<TAC, RI::Tensor<TK>>>& dm)
            -> std::map<int, std::map<TAC, RI::Tensor<TK>>>
            {
                std::map<int, std::map<TAC, RI::Tensor<TK>>> dm_transpose;
                for (const auto& pair0 : dm)
                    for (const auto& pair1 : pair0.second)
                    {
                        const int& iat0 = pair0.first;
                        const int& iat1 = pair1.first.first;
                        const auto& R = pair1.first.second;
                        dm_transpose[iat1][{iat0, { -R[0], -R[1], -R[2] }}] = pair1.second.transpose();
                    }
                return dm_transpose;
            };

        // `cal_force` calculates Pulay term. 
        // If D_IJ = D_KL(H - F = Pulay), it caluclates 0.5 * d(ik | jl).
        // Multiply spin factor (outside) on it gives the final result.
        // The spin factor is not hard-coded in this function.

        // 1. Pulay term
        lr_exx_kernel.set_Ds(dm_gs, this->exx_lri_.lock()->get_info().dm_threshold, spin_suffix);
        lr_exx_kernel.cal_Hs({ "", "", spin_suffix });  // using dm_gs as D_KL
        lr_exx_kernel.cal_force(relaxed_diff_dm, { "", "", spin_suffix , "", "" }); // using relaxed_diff_dm as D_IJ
        add_force_from_kernel();

        // 2. Hellmann-Feynman term
        const auto& dm_gs_transpose = transpose_dm(dm_gs);
        const auto& relaxed_diff_dm_transpose = transpose_dm(relaxed_diff_dm);
        lr_exx_kernel.set_Ds(relaxed_diff_dm_transpose, this->exx_lri_.lock()->get_info().dm_threshold, spin_suffix);
        lr_exx_kernel.cal_Hs({ "", "", spin_suffix });  // using relaxed_diff_dm as D_KL
        lr_exx_kernel.cal_force(dm_gs_transpose, { "", "", spin_suffix , "", "" }); // using dm_gs as D_IJ
        add_force_from_kernel();

        // move back 
        exx_lri_kernel = std::move(lr_exx_kernel);

        // -2 * 0.5 * alpha
        // -2 is the same as post_process_Hexx (a.u. to Ry, which didn't act on Hs)
        // 0.5 is the 2-electron integral prefactor，used in ground-state energy/force where two density matrix are identical
        // But the LR-grad Lagrangian/force 2-e term here Tr[(T+Z)H[D]] does not have 1/2 factor.
        const double fac = -2 * alpha;
        return f_exx_gs_diff * fac;
    }
#endif
}

template class LR::LR_Force<double>;
template class LR::LR_Force<std::complex<double>>;