#include "gradient_inputs.h"
#include "gradient_output.h"
#include "gradient_checks.h"
#include "cal_edm.h"
#include "dm_trans/dm_diff.h"
#include "grad_degen.h"
#include "source_base/timer.h"
#include "source_io/module_output/output_log.h"
namespace LR
{
template <typename T>
std::vector<ModuleBase::matrix> evaluate_open_shell_force(
    const GradientInputs<T>& inputs, LR_Force<T>& lr_force,
    const module_dm::DensityMatrix<T, double>& dm_gs,
    const ct::Tensor& Xz, const ct::Tensor& Z,
    const std::vector<double>& omega, const int label_begin)
{
    ModuleBase::timer::start("LR", "evaluate_open_shell_force");
    const int nst = static_cast<int>(omega.size());
    assert(static_cast<int>(Xz.shape().dim_size(0)) == nst);
    const int ist_begin_ = label_begin;
    const int ist_end_ = label_begin + nst;
    const std::vector<int>& nvirt_g = inputs.nvirt;
    const std::vector<Parallel_2D>& paraX_g = inputs.px;
    const int nloc_g = inputs.nloc;


    const std::vector<int> ld_x = { static_cast<int>(inputs.nk * paraX_g[0].get_local_size()),
                                    static_cast<int>(inputs.nk * paraX_g[1].get_local_size()) };
    const std::vector<int> off_x = { 0, ld_x[0] };
    std::vector<psi::Psi<T>> c_spin;
    for (int is : {0, 1}) { c_spin.push_back(LR_Util::get_psi_spin(inputs.psi_ks, is, inputs.nk)); }

    inputs.ofs << "Start to calculate excited-state force of updown (open shell)" << std::endl;

    const int ist_begin = ist_begin_;
    const int ist_end = ist_end_;
    std::vector<ModuleBase::matrix> forces(ist_end - ist_begin);
    for (int istate = ist_begin;istate < ist_end;++istate)
    {
        const int offset = (istate - ist_begin) * nloc_g;   // X widened into the Z window
        const T* const X_istate = Xz.data<T>() + offset;
        const T* const Z_istate = Z.template data<T>() + offset;

        // 1. the k-space blocks of each spin channel
        std::vector<std::vector<ct::Tensor>> dmx_k(2), dmdiff_k(2), relaxed_k(2);
        for (int is : {0, 1})
        {
#ifdef __MPI
            dmx_k[is] = cal_dm_trans_pblas(X_istate + off_x[is], paraX_g[is], c_spin[is], inputs.pc,
                inputs.nbasis, inputs.nocc[is], nvirt_g[is], inputs.pmat);
            dmdiff_k[is] = cal_dm_diff_pblas(X_istate + off_x[is], paraX_g[is], c_spin[is], inputs.pc,
                inputs.nbasis, inputs.nocc[is], nvirt_g[is], inputs.pmat);
            std::vector<ct::Tensor> dmz_k = cal_dm_trans_pblas(Z_istate + off_x[is], paraX_g[is], c_spin[is],
                inputs.pc, inputs.nbasis, inputs.nocc[is], nvirt_g[is], inputs.pmat);
            for (auto& d : dmz_k) { LR_Util::matsym(d.template data<T>(), inputs.nbasis, inputs.pmat); }
#else
            dmx_k[is] = cal_dm_trans_blas(X_istate + off_x[is], c_spin[is], inputs.nocc[is], nvirt_g[is]);
            dmdiff_k[is] = cal_dm_diff_blas(X_istate + off_x[is], c_spin[is], inputs.nbasis, inputs.nocc[is], nvirt_g[is]);
            std::vector<ct::Tensor> dmz_k = cal_dm_trans_blas(Z_istate + off_x[is], c_spin[is], inputs.nocc[is], nvirt_g[is]);
            for (auto& d : dmz_k) { LR_Util::matsym(d.template data<T>(), inputs.nbasis); }
#endif
            relaxed_k[is] = dmdiff_k[is] + dmz_k;
        }

        // 2. $D^X$. Complex and UN-symmetrized first (the EXX kernel needs the full
        //    non-symmetric $D^X$), then the real symmetrized copy for the grid Hxc force --
        //    `build_dm_from_dmk_spin` symmetrizes IN PLACE, hence the ordering.
        auto dm_trans = LR_Util::build_dm_from_dmk_spin<T, T>(dmx_k,
            inputs.pmat, inputs.nk, inputs.kv.kvec_d, inputs.ucell, inputs.gd, inputs.orb_cutoff);
        LR_Util::transpose_DMR(dm_trans, inputs.pmat);
        auto dm_trans_real = LR_Util::build_dm_from_dmk_spin<T, double>(dmx_k,
            inputs.pmat, inputs.nk, inputs.kv.kvec_d, inputs.ucell, inputs.gd, inputs.orb_cutoff,
            /*symmetrize=*/true);
        LR_Util::transpose_DMR(dm_trans_real, inputs.pmat);

        // 3. the relaxed difference density matrix $T+D^Z$
        const module_dm::DensityMatrix<T, T>& relaxed_diff_dm =
            LR_Util::build_dm_from_dmk_spin<T, T>(relaxed_k,
                inputs.pmat, inputs.nk, inputs.kv.kvec_d, inputs.ucell, inputs.gd, inputs.orb_cutoff);
        module_dm::DensityMatrix<T, double> relaxed_diff_dm_real(&inputs.pmat, 2, inputs.kv.kvec_d, inputs.nk);
        LR_Util::initialize_DMR(relaxed_diff_dm_real, inputs.pmat, inputs.ucell, inputs.gd, inputs.orb_cutoff);
        LR_Util::get_DMR_real_imag_part(relaxed_diff_dm, relaxed_diff_dm_real, 'R');

        // 4. the energy-weighted density matrix
        std::weak_ptr<PotHxcLR> pot_weak = inputs.pot[0];
        const double omega_istate = omega[istate - ist_begin];
        const std::vector<std::vector<ct::Tensor>>& edm_k = cal_edm_from_XZ_istate_openshell(
            inputs, X_istate, Z_istate, omega_istate, inputs.eig_ks.c, dm_trans, pot_weak);
        module_dm::DensityMatrix<T, double> edm_real = LR_Util::build_dm_from_dmk_spin<T, double>(edm_k,
            inputs.pmat, inputs.nk, inputs.kv.kvec_d, inputs.ucell, inputs.gd, inputs.orb_cutoff,
            /*symmetrize=*/true);

        // 5. the force terms
        ModuleBase::matrix force_hxc_dmtrans = lr_force.cal_force_hxc_dmtrans(dm_trans_real, *inputs.pot[0]);
        if (inputs.test_force)
            ModuleIO::print_force(inputs.ofs, inputs.ucell, "HXC DMTRANS FORCE (eV/Angstrom)", force_hxc_dmtrans, false);


        // the $g^{xc}$ half of $\partial_x K[D^X]D^X$, i.e. the derivative of the xc kernel
        // through the ground-state density. Only for local kernels.
        if (LR_Util::has_local_xc(inputs.xc_kernel))
        {
            // The density dependence of K[D^X]D^X belongs to the LR functional.
            const std::shared_ptr<PotHxcLR>& pot_lr = inputs.pot[0];
            PotGradXCLR pot_grad(pot_lr->xc_kernel_components(), pot_lr->get_rho_basis(),
                inputs.ucell, pot_lr->nrxx, /*triplet=*/false);
            ModuleBase::matrix force_gxc_dmtrans =
                lr_force.cal_force_gxc_dmtrans_openshell(dm_trans_real, dm_gs, pot_grad);
            if (inputs.test_force)
                ModuleIO::print_force(inputs.ofs, inputs.ucell, "GXC DMTRANS FORCE (eV/Angstrom)", force_gxc_dmtrans, false);
            force_hxc_dmtrans += force_gxc_dmtrans;
        }

        ModuleBase::matrix force_hamiltgs_relaxed_diff = lr_force.cal_force_hamilt_gs_dm_relaxed_diff(
            relaxed_diff_dm_real, dm_gs, false, inputs.pot_hxc_gs.get());
        if (inputs.test_force)
            ModuleIO::print_force(inputs.ofs, inputs.ucell, "H_GS-(T+Z) FORCE (without EXX) (eV/Angstrom)", force_hamiltgs_relaxed_diff, false);

        ModuleBase::matrix force_overlap_edm = lr_force.cal_force_overlap_edm(edm_real);
        if (inputs.test_force)
            ModuleIO::print_force(inputs.ofs, inputs.ucell, "OVERLAP-EDM FORCE (eV/Angstrom)", force_overlap_edm, false);

#ifdef __EXX
        const double& alpha = inputs.hybrid_alpha;
        // Exchange is spin-diagonal, so each channel is done independently and summed.
        // `get_exx_Ds_gs` returns the channels unscaled (SPIN_multiple = 1 at nspin=2), unlike
        // the closed-shell `get_exx_Ds_spin1` which returns 0.5*D. Counting the closed-shell
        // `alpha*4.0` back down for both terms:
        //   $D^XD^X$: each slot goes 0.5*D^X_tot -> D^X_is, i.e. x2 each, and the explicit
        //     sum over is adds another x2 -- but $D^X_\text{tot}=\sqrt2 D^X_\sigma$ eats one,
        //     so 4/(2*2) * ... = `alpha`.
        //   $D^\text{gs}(T{+}D^Z)$: the left slot goes 0.5*D_up -> D_is (x2); the right slot
        //     goes 0.5*(T+Z)_tot = (T+Z)_up -> (T+Z)_is (x1, no sqrt2 here); the explicit sum
        //     over is adds x2. So 4/(2*1*2) = `alpha` as well -- NOT `2*alpha`: the earlier
        //     comment forgot that the closed-shell right slot is already the spin SUM, which is
        //     exactly what the `for (is)` loop below now supplies.
        if (LR::exx_kernel_list().count(inputs.xc_kernel))
        {
            const auto& Ds_trans = LR_Util::get_exx_Ds_gs(dm_trans, inputs.ucell, inputs.kv, inputs.pmat);
            ModuleBase::matrix force_exx_dmtrans(inputs.ucell.nat, 3);
            for (int is : {0, 1})
            {
                force_exx_dmtrans += lr_force.cal_force_exx_dm_trans(Ds_trans[is], alpha, std::to_string(is));
            }
            if (inputs.test_force)
                ModuleIO::print_force(inputs.ofs, inputs.ucell, "EXX DMTRANS FORCE (eV/Angstrom)", force_exx_dmtrans, false);
            force_hxc_dmtrans += force_exx_dmtrans;
        }
        if (LR::gs_is_hybrid(inputs.dft_functional))
        {
            const auto& Ds_gs = LR_Util::get_exx_Ds_gs(dm_gs, inputs.ucell, inputs.kv, inputs.pmat);
            const auto& Ds_relaxed_diff = LR_Util::get_exx_Ds_gs(relaxed_diff_dm, inputs.ucell, inputs.kv, inputs.pmat);
            ModuleBase::matrix force_exx_gs_relaxed_diff(inputs.ucell.nat, 3);
            for (int is : {0, 1})
            {
                force_exx_gs_relaxed_diff += lr_force.cal_force_exx_gs_dm_relaxed_diff(
                    Ds_gs[is], Ds_relaxed_diff[is], alpha, std::to_string(is));
            }
            if (inputs.test_force)
                ModuleIO::print_force(inputs.ofs, inputs.ucell, "EXX GS-(T+Z) FORCE (eV/Angstrom)", force_exx_gs_relaxed_diff, false);
            force_hamiltgs_relaxed_diff += force_exx_gs_relaxed_diff;
        }
#endif
        forces[istate - ist_begin] = force_hxc_dmtrans + force_hamiltgs_relaxed_diff + force_overlap_edm;
    }
    print_lr_force(forces, std::cout, ist_begin);
    print_lr_force(forces, inputs.ofs, ist_begin);
    ModuleBase::timer::end("LR", "evaluate_open_shell_force");
    return forces;
}


template std::vector<ModuleBase::matrix> evaluate_open_shell_force<double>(const GradientInputs<double>&, LR_Force<double>&, const module_dm::DensityMatrix<double, double>&, const ct::Tensor&, const ct::Tensor&, const std::vector<double>&, int);
template std::vector<ModuleBase::matrix> evaluate_open_shell_force<std::complex<double>>(const GradientInputs<std::complex<double>>&, LR_Force<std::complex<double>>&, const module_dm::DensityMatrix<std::complex<double>, double>&, const ct::Tensor&, const ct::Tensor&, const std::vector<double>&, int);
}
