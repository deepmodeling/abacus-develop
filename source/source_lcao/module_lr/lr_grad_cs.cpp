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
std::vector<ModuleBase::matrix> evaluate_closed_shell_force(
    const GradientInputs<T>& inputs, LR_Force<T>& lr_force,
    const module_dm::DensityMatrix<T, double>& dm_gs,
    const ct::Tensor& Xz, const ct::Tensor& Z,
    const std::vector<double>& omega, const int label_begin, const int ispin)
{
    ModuleBase::timer::start("LR", "evaluate_closed_shell_force");
    // for each block, calculate dm_trans, dm_relaxed_diff, edm and force
    const int nst = static_cast<int>(omega.size());
    assert(static_cast<int>(Xz.shape().dim_size(0)) == nst);
    // `ist_begin`/`ist_end` label the blocks in the output only; the gradient itself never looks
    // up a state, it only uses `omega[i]`. That is what lets a caller pass excitation vectors
    // that are not the stored eigenvectors (see `cal_grad_matrix_degenerate`).
    const int ist_begin = label_begin;
    const int ist_end = label_begin + nst;
    const std::vector<int>& nvirt_g = inputs.nvirt;
    const std::vector<Parallel_2D>& paraX_g = inputs.px;
    const int nloc_g = inputs.nloc;


    ModuleBase::TITLE("ESolver_LR", "cal_force");
    // Spin channel 0, NOT `ispin`. `ispin` indexes `spin_types` = {singlet, triplet}.
    // Closed shell always use spin-up channel of psi_ks, i.e. psi_ks(0).
    const auto& c = LR_Util::get_psi_spin(inputs.psi_ks, 0, inputs.nk);   // wavefunction coefficients of ground state

    // calculate the force (the partial gradient of Lagrangian)
    inputs.ofs << "Start to calculate excited-state force of " << inputs.spin_types[ispin] << std::endl;
    // ground state dm for currrent spin (only for test the correctness of the force)
    // module_dm::DensityMatrix<T, T> dm_gs(inputs.pmat, 1, inputs.kv.kvec_d, inputs.nk);

    std::vector<ModuleBase::matrix> forces(ist_end - ist_begin);
    for (int istate = ist_begin;istate < ist_end;++istate)
    {
        const int offset = (istate - ist_begin) * nloc_g;   // block of X (widened into the Z window)
        const int zoffset = offset;                         // block of Z
        // The imag part will be cancelled in the force calculation, so we use double DM(R) to calculate force.
        // But complex transition DM(R) is still used in energy density matrix calculation.
#ifdef __MPI
        const auto& dm_trans_k = cal_dm_trans_pblas(Xz.data<T>() + offset, paraX_g[ispin], c, inputs.pc, inputs.nbasis, inputs.nocc[ispin], nvirt_g[ispin], inputs.pmat);
#else
        const auto& dm_trans_k = cal_dm_trans_blas(Xz.data<T>() + offset, c, inputs.nocc[ispin], nvirt_g[ispin]);
#endif
        // D(X) complex, for the EXX (LibRI) force. Built FIRST and left UN-symmetrized:
        // the exchange kernel (mu kappa | nu lambda) puts the two indices of one D^X into
        // different electron coordinates, so Tr[D^X D^X K_exx] = (aa|ii) requires the full
        // non-symmetric D^X. Symmetrizing would give 1/2[(aa|ii)+(ai|ia)], which is wrong.
        // (`cal_force_exx_dm_trans` feeds the same tensor to both slots, so it is consistent.)
        auto dm_trans =   // D(X) complex
            LR_Util::build_dm_from_dmk<T, T>(dm_trans_k,
                inputs.pmat, inputs.nk, inputs.kv.kvec_d, inputs.ucell, inputs.gd, inputs.orb_cutoff);
        LR_Util::transpose_DMR(dm_trans, inputs.pmat);
        // D(X) real, for the grid Hxc force. The Coulomb kernel (mu nu | kappa lambda) is
        // symmetric within each index pair, so it only ever sees the symmetric part of D^X.
        // In `PulayForceStress::cal_pulay_fs`, `cal_gint_rho` (which builds v) symmetrizes
        // implicitly, while `cal_gint_fvl`'s internal factor 2 assumes D_{mu nu} = D_{nu mu}.
        // Passing an un-symmetrized D^X makes the two slots of the bilinear form disagree.
        // NOTE: `build_dm_from_dmk` symmetrizes `dm_trans_k` IN PLACE, hence the ordering.
        auto dm_trans_real =   // D(X), double (FIXME: not enough for periodic system!)
            LR_Util::build_dm_from_dmk<T, double>(dm_trans_k,
                inputs.pmat, inputs.nk, inputs.kv.kvec_d, inputs.ucell, inputs.gd, inputs.orb_cutoff,
                /*symmetrize=*/true);
        LR_Util::transpose_DMR(dm_trans_real, inputs.pmat);
        // LR_Util::print_DMR(dm_trans, "dm_trans of istate " + std::to_string(istate));
        // difference density matrix
#ifdef __MPI
        std::vector<ct::Tensor> dm_diff_k = cal_dm_diff_pblas(Xz.data<T>() + offset, paraX_g[ispin], c, inputs.pc, inputs.nbasis, inputs.nocc[ispin], nvirt_g[ispin], inputs.pmat);
#else
        std::vector<ct::Tensor> dm_diff_k = cal_dm_diff_blas(Xz.data<T>() + offset, c, inputs.nbasis, inputs.nocc[ispin], nvirt_g[ispin]);
#endif
        // std::cout << "dm_diff_k T(k) before symmetrization, istate " + std::to_string(istate) << std::endl;
        // LR_Util::print_value(dm_diff_k[0].data<T>(), inputs.pmat.get_col_size(), inputs.pmat.get_row_size());
        // for (auto& d : dm_diff_k) { LR_Util::matsym(d.data<T>(), inputs.nbasis, inputs.pmat); }   // symmetrize
        // std::cout << "dm_diff_k T(k) after symmetrization, istate " + std::to_string(istate) << std::endl;
        // LR_Util::print_value(dm_diff_k[0].data<T>(), inputs.pmat.get_col_size(), inputs.pmat.get_row_size());

#ifdef __MPI
        const std::vector<ct::Tensor>& dm_relaxed_k = cal_dm_trans_pblas(Z.template data<T>() + zoffset, paraX_g[ispin], c, inputs.pc, inputs.nbasis, inputs.nocc[ispin], nvirt_g[ispin], inputs.pmat);
#else
        const std::vector<ct::Tensor>& dm_relaxed_k = cal_dm_trans_blas(Z.template data<T>() + zoffset, c, inputs.nocc[ispin], nvirt_g[ispin]);
#endif
        // std::cout << "dm_relaxed_k Z(k) before symmetrization, istate " + std::to_string(istate) << std::endl;
        // LR_Util::print_value(dm_relaxed_k[0].data<T>(), inputs.pmat.get_col_size(), inputs.pmat.get_row_size());
#ifdef __MPI
        for (auto& d : dm_relaxed_k) { LR_Util::matsym(d.data<T>(), inputs.nbasis, inputs.pmat); }    // symmetrize
#else
        for (auto& d : dm_relaxed_k) { LR_Util::matsym(d.data<T>(), inputs.nbasis); }    // symmetrize
#endif
        // std::cout << "dm_relaxed_k Z(k) after symmetrization, istate " + std::to_string(istate) << std::endl;
        // LR_Util::print_value(dm_relaxed_k[0].data<T>(), inputs.pmat.get_col_size(), inputs.pmat.get_row_size());
        // relaxed difference density matrix
        const std::vector<ct::Tensor>& relaxed_diff_dm_k = dm_diff_k + dm_relaxed_k;
        const module_dm::DensityMatrix<T, T>& diff_dm =
            LR_Util::build_dm_from_dmk<T, T>(dm_diff_k,
                inputs.pmat, inputs.nk, inputs.kv.kvec_d, inputs.ucell, inputs.gd, inputs.orb_cutoff);
        const module_dm::DensityMatrix<T, T>& relaxed_diff_dm =
            LR_Util::build_dm_from_dmk<T, T>(relaxed_diff_dm_k,
                inputs.pmat, inputs.nk, inputs.kv.kvec_d, inputs.ucell, inputs.gd, inputs.orb_cutoff);
        // LR_Util::print_DMR(relaxed_diff_dm, "relaxed_diff_dm T+Z (Z symmetrized) of istate " + std::to_string(istate));

        // module_dm::DensityMatrix<T, T> relaxed_diff_dm =    // T+D(Z), (R) can be complex
        //     LR_Util::build_dm_from_dmk<T, T>(
        //         // LR_Util::operator+(
        //             cal_dm_diff_pb las(Xz.data<T>() + offset, paraX_g[ispin], c, inputs.pc, inputs.nbasis, inputs.nocc[ispin], nvirt_g[ispin], inputs.pmat)
        //             + cal_dm_trans_pblas(Z.template data<T>() + offset, paraX_g[ispin], c, inputs.pc, inputs.nbasis, inputs.nocc[ispin], nvirt_g[ispin], inputs.pmat)
        //             ,// ),
        //         inputs.pmat, inputs.nk, inputs.kv.kvec_d, inputs.ucell, inputs.gd, inputs.orb_cutoff);
        // LR_Util::print_DMR(relaxed_diff_dm, "relaxed_diff_dm of istate " + std::to_string(istate));
        module_dm::DensityMatrix<T, double> relaxed_diff_dm_real(&inputs.pmat, 1, inputs.kv.kvec_d, inputs.nk);
        LR_Util::initialize_DMR(relaxed_diff_dm_real, inputs.pmat, inputs.ucell, inputs.gd, inputs.orb_cutoff);
        LR_Util::get_DMR_real_imag_part(relaxed_diff_dm, relaxed_diff_dm_real, 'R');

        // get edm of type DensityMatrix
        // weak_ptr here is to avoid "could not match 'weak_ptr' against 'shared_ptr'"
        // but why there're no bug in the previous code (HamiltLR and HamiltULF)?
        // seems because those two are classes having constructors
        // but `cal_edm_from_XZ_istate` here is a functions
        std::weak_ptr<PotHxcLR> pot_weak = inputs.pot[ispin];
        const T* const x_istate = Xz.data<T>() + offset;
        const T* const z_istate = Z.template data<T>() + zoffset;
        const double omega_istate = omega[istate - ist_begin];
        const std::vector<ct::Tensor>& edm_k = cal_edm_from_XZ_istate(
            inputs, x_istate, z_istate, omega_istate, inputs.eig_ks.c,
            dm_trans, c, pot_weak, inputs.spin_types[ispin]);
        if (inputs.test_force && inputs.nocc[0] == 1 && nvirt_g[0] == 1)
        {
#ifdef __MPI
            const std::vector<ct::Tensor>& dm_diff = cal_dm_diff_pblas(Xz.data<T>() + offset, paraX_g[0], c, inputs.pc, inputs.nbasis, inputs.nocc[0], nvirt_g[0], inputs.pmat);
#else
            const std::vector<ct::Tensor>& dm_diff = cal_dm_diff_blas(Xz.data<T>() + offset, c, inputs.nbasis, inputs.nocc[0], nvirt_g[0]);
#endif
            // test_dm_diff_H2<T>(relaxed_diff_dm.get_dmk_ptr(0), c, inputs.nbasis);
            test_dm_diff_H2<T>(dm_diff[0].data<T>(), c, inputs.nbasis);
            test_edm_H2<T>(edm_k[0].data<T>(), inputs.eig_ks.c, c, inputs.nbasis);
        }
        module_dm::DensityMatrix<T, double> edm_real = LR_Util::build_dm_from_dmk<T, double>(edm_k,
            inputs.pmat, inputs.nk, inputs.kv.kvec_d, inputs.ucell, inputs.gd, inputs.orb_cutoff, /*symmetrize=*/true);
        // print edm_real (R)
        if (inputs.test_force)
        {
            LR_Util::save_DMR(edm_real, "data-EDMR-sparse" + std::string(inputs.excited_relax ? "_state" + std::to_string(istate) : ""), inputs.pmat, inputs.out_dir, inputs.nbasis, inputs.my_rank);
            // LR_Util::print_DMR(edm_real, "edm_real (R) of istate " + std::to_string(istate));
        }

        ModuleBase::matrix force_hxc_dmtrans = lr_force.cal_force_hxc_dmtrans(dm_trans_real, *inputs.pot[ispin]);
        if (inputs.test_force)
            ModuleIO::print_force(inputs.ofs, inputs.ucell, "HXC DMTRANS FORCE (eV/Angstrom)", force_hxc_dmtrans, false);


        // the $g^{xc}$ half of $\partial_x K[D^X]D^X$, i.e. the derivative of the xc kernel through
        // the ground-state density (see `cal_force_gxc_dmtrans`). Only for local kernels.
        if (LR_Util::has_local_xc(inputs.xc_kernel))
        {
            // The density dependence of K[D^X]D^X belongs to the LR functional.
            const std::shared_ptr<PotHxcLR>& pot_lr = inputs.pot[ispin];
            const bool triplet = inputs.spin_types[ispin] == "triplet";
            PotGradXCLR pot_grad(pot_lr->xc_kernel_components(), pot_lr->get_rho_basis(),
                inputs.ucell, pot_lr->nrxx, triplet);
            ModuleBase::matrix force_gxc_dmtrans = lr_force.cal_force_gxc_dmtrans(dm_trans_real, dm_gs, pot_grad);
            if (inputs.test_force)
                ModuleIO::print_force(inputs.ofs, inputs.ucell, "GXC DMTRANS FORCE (eV/Angstrom)", force_gxc_dmtrans, false);
            force_hxc_dmtrans += force_gxc_dmtrans;
        }
        ModuleBase::matrix force_hamiltgs_relaxed_diff = lr_force.cal_force_hamilt_gs_dm_relaxed_diff(relaxed_diff_dm_real, dm_gs, false, inputs.pot_hxc_gs.get());
        if (inputs.test_force)
            ModuleIO::print_force(inputs.ofs, inputs.ucell, "H_GS-(T+Z) FORCE (without EXX) (eV/Angstrom)", force_hamiltgs_relaxed_diff, false);

        ModuleBase::matrix force_overlap_edm = lr_force.cal_force_overlap_edm(edm_real);    // "-" sign has been included in the force factor
        if (inputs.test_force)
            ModuleIO::print_force(inputs.ofs, inputs.ucell, "OVERLAP-EDM FORCE (eV/Angstrom)", force_overlap_edm, false);

        if (inputs.test_force)
        {
            // test H[T] force (Z=0), non-EXX part
            module_dm::DensityMatrix<T, double> diff_dm_real(&inputs.pmat, 1, inputs.kv.kvec_d, inputs.nk);
            LR_Util::initialize_DMR(diff_dm_real, inputs.pmat, inputs.ucell, inputs.gd, inputs.orb_cutoff);
            LR_Util::get_DMR_real_imag_part(diff_dm, diff_dm_real, 'R');

            inputs.ofs << "========== [TEST H_GS-(T) force (Z=0), non-EXX part] ===========" << std::endl;
            ModuleBase::matrix force_hamiltgs_diff = lr_force.cal_force_hamilt_gs_dm_relaxed_diff(diff_dm_real, dm_gs);
            ModuleIO::print_force(inputs.ofs, inputs.ucell, "H_GS-T FORCE (without EXX) (eV/Angstrom)", force_hamiltgs_diff, false);
            inputs.ofs << "========== [\\TEST H_GS-(T) force (Z=0), non-EXX part] ===========" << std::endl;
        }


#ifdef __EXX
        const double& alpha = inputs.hybrid_alpha;

        if (LR::exx_kernel_list().count(inputs.xc_kernel))
        {
            const auto& Ds_trans = LR_Util::get_exx_Ds_spin1(dm_trans, inputs.ucell, inputs.kv, inputs.pmat);
            ModuleBase::matrix force_exx_dmtrans = lr_force.cal_force_exx_dm_trans(Ds_trans, alpha * 4.0);  // cancel the two 0.5s in Ds
            if (inputs.test_force)
                ModuleIO::print_force(inputs.ofs, inputs.ucell, "EXX DMTRANS FORCE (eV/Angstrom)", force_exx_dmtrans, false);
            force_hxc_dmtrans += force_exx_dmtrans;

        }

        if (LR::gs_is_hybrid(inputs.dft_functional))
        {
            const auto& Ds_gs = LR_Util::get_exx_Ds_spin1(dm_gs, inputs.ucell, inputs.kv, inputs.pmat);    // returns 0.5*D[0]
            const auto& Ds_relaxed_diff = LR_Util::get_exx_Ds_spin1(relaxed_diff_dm, inputs.ucell, inputs.kv, inputs.pmat);   // returns 0.5*D[0]
            // LR_Util::print_CV(Ds_relaxed_diff, "Ds_relaxed_diff for EXX force");
            // `get_exx_Ds_spin1` feeds `split_m2D_ktoR(..., nspin=1)`, which reads only channel 0
            // with a 0.5 prefactor. For `dm_gs` that channel is $D^\text{gs}_\uparrow$ at nspin=2
            // but the spin-summed $D^\text{gs}$ at nspin=1, i.e. twice as large.
            ModuleBase::matrix force_exx_gs_relaxed_diff = lr_force.cal_force_exx_gs_dm_relaxed_diff(Ds_gs, Ds_relaxed_diff, alpha * 4.0) * gs_dm_channel_factor(inputs.nspin);  // cancel the two 0.5s in Ds
            if (inputs.test_force)
                ModuleIO::print_force(inputs.ofs, inputs.ucell, "EXX GS-(T+Z) FORCE (eV/Angstrom)", force_exx_gs_relaxed_diff, false);
            force_hamiltgs_relaxed_diff += force_exx_gs_relaxed_diff;

            if (inputs.test_force)
            {
                // test H[T] force (Z=0), EXX part
                const auto& Ds_diff = LR_Util::get_exx_Ds_spin1(diff_dm, inputs.ucell, inputs.kv, inputs.pmat);   // returns 0.5*D[0]
                inputs.ofs << "========== [TEST H_GS-(T) force (Z=0), EXX part] ===========" << std::endl;
                ModuleBase::matrix force_exx_gs_diff = lr_force.cal_force_exx_gs_dm_relaxed_diff(Ds_gs, Ds_diff, alpha * 4.0) * gs_dm_channel_factor(inputs.nspin);  // cancel the two 0.5s in Ds
                ModuleIO::print_force(inputs.ofs, inputs.ucell, "H_GS-T EXX FORCE (Z=0) (eV/Angstrom)", force_exx_gs_diff, false);
                inputs.ofs << "========== [\\TEST H_GS-(T) force (Z=0), EXX part] ===========" << std::endl;
            }
        }
#endif
        forces[istate - ist_begin] = force_hxc_dmtrans + force_hamiltgs_relaxed_diff + force_overlap_edm;
    }
    // total force
    print_lr_force(forces, std::cout, ist_begin);
    print_lr_force(forces, inputs.ofs, ist_begin);
    ModuleBase::timer::end("LR", "evaluate_closed_shell_force");
    return forces;
}


template std::vector<ModuleBase::matrix> evaluate_closed_shell_force<double>(const GradientInputs<double>&, LR_Force<double>&, const module_dm::DensityMatrix<double, double>&, const ct::Tensor&, const ct::Tensor&, const std::vector<double>&, int, int);
template std::vector<ModuleBase::matrix> evaluate_closed_shell_force<std::complex<double>>(const GradientInputs<std::complex<double>>&, LR_Force<std::complex<double>>&, const module_dm::DensityMatrix<std::complex<double>, double>&, const ct::Tensor&, const ct::Tensor&, const std::vector<double>&, int, int);
}
