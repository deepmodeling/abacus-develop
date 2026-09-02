#include "lr_force.h"
#include "cal_hs_grad.h"
#include "pulay_force_hcontainer.h"
#include "source_lcao/pulay_fs.h"   // only for gint terms
#include "source_hamilt/module_gint/gint_interface.h"
#include "source_lcao/module_lr/utils/lr_util.h"
// #include "source_lcao/module_lr/utils/lr_util_hcontainer.h"
namespace LR
{
    /// The LR density matrices ($D^X$, $T+D^Z$, EDM) carry one channel in the closed-shell
    /// singlet/triplet algorithm and two independent channels in the open-shell one.
    template<typename TK>
    inline bool is_openshell_dm(const elecstate::DensityMatrix<TK, double>& dm)
    {
        return dm.get_DMR_vector().size() == 2;
    }

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
        const bool reproduce_gs,
        const PotHxcLR* pot_hxc_gs)
    {
        const bool with_ewald = reproduce_gs;
        // Two independent spin channels in `relax_diff_dm` <=> open-shell (spin-unrestricted) LR.
        // The closed-shell singlet/triplet algorithm always builds a single-channel LR density
        // matrix, even at nspin=2.
        const bool openshell = is_openshell_dm(relax_diff_dm);
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
        // The potential here must be the *linear response* of $V^\text{Hxc}$ to the difference
        // density, i.e. $v_H[\rho^{T+Z}] + f_{xc}[\rho^\text{gs}]\,\rho^{T+Z}$ -- NOT
        // $v_\text{Hxc}[\rho^{T+Z}]$. This term is
        //     $\sum_{\kappa\lambda}(T{+}D^Z)_{\kappa\lambda}\int\phi_\kappa\phi_\lambda\,
        //      f_{xc}\sum_{\alpha\beta}D^\text{gs}_{\alpha\beta}(\phi_\alpha\phi_\beta)^x$,
        // the half of $\partial_x V^\text{Hxc}$ whose basis derivative falls on the *ground-state*
        // pair. Hartree is linear in the density so feeding it $\rho^{T+Z}$ happens to be right;
        // xc is not -- $\rho^{T+Z}$ is not even positive everywhere, while LDA has
        // $v_{xc}\propto-\rho^{1/3}$.
        //
        // `pot_hxc_gs` supplies exactly this object (Hartree weight 1, xc = $(f_{uu}+f_{ud})/2$ at
        // nspin=2), with no extra factor. Verified on H2/SZ TDRPA@LDA, where $K^T\equiv0$ makes the
        // triplet gradient identical to $d(\varepsilon_a-\varepsilon_i)/dx$: analytic 28.5982 vs the
        // KS-gap finite difference 28.598156. It used to be off by -3.5 eV/Ang.
        //
        // `reproduce_gs` is the exception: there `relax_diff_dm` *is* the ground-state density
        // matrix and the term being checked is the true ground-state force, for which
        // $v_\text{Hxc}[\rho^\text{gs}]$ is the correct potential.
        if (reproduce_gs || pot_hxc_gs == nullptr)
        {
            elecstate::Potential pot_hxc_relaxed_diff = this->dm_to_hxc_potential(relax_diff_dm);
            //`cal_pulay_fs` calculates only one spin channel because `relax_diff_dm` has only one.
            PulayForceStress::cal_pulay_fs(1/*nspin*/, fhxc_dvhxc, stress_tmp,
                dm_gs, this->ucell_, &pot_hxc_relaxed_diff, true, false);
            fhxc_dvhxc *= gs_dm_channel_factor();
        }
        else if (openshell)
        {
            // $v_\sigma=\sum_{\sigma'}f^{\sigma\sigma'}\rho^{T+Z}_{\sigma'}$ (+ the full Hartree),
            // contracted with $D^\text{gs}_\sigma$. No factor 2 and no `gs_dm_channel_factor`:
            // the two ground-state channels are summed explicitly by `cal_gint_fvl` below, and
            // each carries occupation 1 rather than 2.
            constexpr int nspin_dm = 2;
            std::vector<ModuleBase::matrix> v_lin(nspin_dm, ModuleBase::matrix(1, this->rhopw_.nrxx));
            for (int sl = 0; sl < nspin_dm; ++sl)
            {
                for (int sr = 0; sr < nspin_dm; ++sr)
                {
                    double* rho_in[1] = { const_cast<double*>(chr_diff_relaxed.rho[sr]) };
                    pot_hxc_gs->cal_v_eff(rho_in, this->ucell_, v_lin[sl], { sl, sr });
                }
            }
            std::vector<const double*> vr_eff(nspin_dm);
            for (int is = 0; is < nspin_dm; ++is) { vr_eff[is] = v_lin[is].c; }
            ModuleGint::cal_gint_fvl(nspin_dm, vr_eff, dm_gs.get_DMR_vector(), true, false, &fhxc_dvhxc, &stress_tmp);
        }
        else
        {
            ModuleBase::matrix v_lin(1, this->rhopw_.nrxx);   // zero-initialized
            double* rho_in[1] = { const_cast<double*>(chr_diff_relaxed.rho[0]) };
            pot_hxc_gs->cal_v_eff(rho_in, this->ucell_, v_lin);
            std::vector<const double*> vr_eff = { v_lin.c };
            ModuleGint::cal_gint_fvl(1, vr_eff, dm_gs.get_DMR_vector(), true, false, &fhxc_dvhxc, &stress_tmp);
            fhxc_dvhxc *= 2;   // for the two channels of the ground-state dm.
            fhxc_dvhxc *= gs_dm_channel_factor();
        }

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
        // `dm_trans` (D^X) must be SYMMETRIZED before entering here: `cal_pulay_fs` builds v from
        // rho[D^X] (which only sees the symmetric part) but contracts with D^X as passed, so an
        // un-symmetrized D^X makes the two slots of the bilinear form Tr[D^X d(K_H)[D^X]] disagree.
        //
        // `cal_pulay_fs` returns 2 * sum_{mn} D_{mn} \int (d phi_m) v phi_n, where the factor 2 is
        // `cal_gint_fvl`'s internal m<->n doubling, i.e. it is exactly the *bra-pair* derivative (Pulay).
        // The *ket-pair* derivative (Hellmann-Feynman) is equal to it (both slots hold the same D^X), 
        // so the total needs one more factor 2 (Pulay -> Pulay + Hellmann-Feynman).
        // Verified on H2/SZ against the analytic 4-center derivative: 4*sum D^sym P = 16.8653548
        // vs 2*d(ai|ia)/dz = 16.865355 eV/Ang (7 digits).
        const double pulay_to_total_sym = 2.0;
        // Open shell: the kernel couples the channels, so the density has to be built per channel
        // and the potential accumulated over the summed spin before contracting with $D^X_\sigma$.
        // `pulay_to_total_sym` is a Pulay -> Pulay+Hellmann-Feynman factor and is spin-independent.
        if (is_openshell_dm(dm_trans))
        {
            return PulayForceStress::cal_pulay_fs_openshell(dm_trans, this->ucell_, &pot_hxc) * pulay_to_total_sym;
        }
        return PulayForceStress::cal_pulay_fs(dm_trans, this->ucell_, &pot_hxc) * pulay_to_total_sym;
    }

    template<typename TK>
    ModuleBase::matrix LR_Force<TK>::cal_force_gxc_dmtrans(const elecstate::DensityMatrix<TK, double>& dm_trans,
        const elecstate::DensityMatrix<TK, double>& dm_gs, const PotGradXCLR& pot_grad)
    {
        // The third source of position dependence in
        //     $f^{xc}_{\kappa\lambda,\alpha\beta}
        //        =\int\phi_\kappa\phi_\lambda\,f_{xc}[\rho^\text{gs}(r)]\,\phi_\alpha\phi_\beta$.
        // `cal_force_hxc_dmtrans` above differentiates the two basis pairs of $D^X$, which for the
        // Hartree kernel $(\kappa\lambda|\alpha\beta)$ is everything. The xc kernel additionally
        // depends on the nuclear positions through $\rho^\text{gs}$ itself, and that derivative is
        //     $\int\rho^X g^{xc}\rho^X\,\partial_x\rho^\text{gs}|_\text{basis}
        //      =\int v^{(2)}[\rho^X,\rho^X]\,\partial_x\rho^\text{gs}|_\text{basis}$,
        // i.e. the Pulay derivative of the *ground-state* density against the $g^{xc}$ potential.
        //
        // No spin factor: this is the sibling of `cal_force_hxc_dmtrans` above, which likewise has
        // none, because `PotGradXCLR` (like `pot[ispin]` there) already carries the S2_singlet /
        // S2_triplet spin combination. The superficially similar Hellmann-Feynman half of
        // `cal_force_hamilt_gs_dm_relaxed_diff` DOES need a factor 2, but only because it is built
        // from `pot_hxc_gs`, which is normalized as S2_gs = S2_singlet/2.
        // Confirmed numerically on H2/SZ TDLDA (see `cal_multiplier_w_from_z.h`).
        const Charge chr_x = dm_to_charge(dm_trans);
        ModuleBase::matrix v2(1, this->rhopw_.nrxx);   // zero-initialized
        double* rho_in[1] = { const_cast<double*>(chr_x.rho[0]) };
        pot_grad.cal_v_eff(rho_in, this->ucell_, v2);

        ModuleBase::matrix f(this->ucell_.nat, 3);
        ModuleBase::matrix stress_tmp;
        std::vector<const double*> vr_eff = { v2.c };
        ModuleGint::cal_gint_fvl(1, vr_eff, dm_gs.get_DMR_vector(), true, false, &f, &stress_tmp);
        f *= gs_dm_channel_factor();
        return f;
    }

    template<typename TK>
    ModuleBase::matrix LR_Force<TK>::cal_force_gxc_dmtrans_openshell(
        const elecstate::DensityMatrix<TK, double>& dm_trans,
        const elecstate::DensityMatrix<TK, double>& dm_gs, const PotGradXCLR& pot_grad)
    {
        // Same term as `cal_force_gxc_dmtrans`, spin-resolved: the free index $\tau$ (spin channel)
        // of $v^{(2)}_\tau$ is contracted with $D^\text{gs}_\tau$, and `cal_gint_fvl` does the
        // $\sum_\tau$. No `gs_dm_channel_factor` here -- both ground-state channels are summed
        // explicitly and each carries occupation 1.
        constexpr int nspin_dm = 2;
        assert(dm_trans.get_DMR_vector().size() == nspin_dm);
        const Charge chr_x = dm_to_charge(dm_trans);
        const double* rho_in[nspin_dm] = { chr_x.rho[0], chr_x.rho[1] };

        std::vector<ModuleBase::matrix> v2(nspin_dm, ModuleBase::matrix(1, this->rhopw_.nrxx));
        for (int tau = 0; tau < nspin_dm; ++tau) { pot_grad.cal_v_eff_openshell(rho_in, this->ucell_, v2[tau], tau); }

        ModuleBase::matrix f(this->ucell_.nat, 3);
        ModuleBase::matrix stress_tmp;
        std::vector<const double*> vr_eff(nspin_dm);
        for (int is = 0; is < nspin_dm; ++is) { vr_eff[is] = v2[is].c; }
        ModuleGint::cal_gint_fvl(nspin_dm, vr_eff, dm_gs.get_DMR_vector(), true, false, &f, &stress_tmp);
        return f;
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