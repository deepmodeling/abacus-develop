#include "source_lcao/module_lr/lr_core.h"

#include "source_base/module_external/scalapack_connector.h"
#include "source_base/timer.h"
#include "source_base/tool_quit.h"
#include "source_base/tool_title.h"
#include "source_basis/module_nao/two_center_bundle.h"
#include "source_basis/module_nao/two_center_integrator.h"
#include "source_cell/klist.h"
#include "source_cell/module_neighbor/sltk_grid_driver.h"
#include "source_cell/unitcell.h"
#include "source_estate/module_dm/dm_from_psi.h"
#include "source_io/module_output/output_log.h"
#include "source_io/module_parameter/input_parameter.h"
#include "source_lcao/module_lr/gradient_inputs.h"
#include "source_lcao/module_lr/gradient_output.h"
#include "source_lcao/module_lr/lr_amp.h"
#include "source_lcao/module_lr/lr_dm.h"
#include "source_lcao/module_lr/lr_force.h"
#include "source_lcao/module_lr/lr_gs.h"
#include "source_lcao/module_lr/lr_sforce.h"
#include "source_lcao/module_lr/lr_zeq.h"
#include "source_lcao/module_lr/utils/lr_util.h"
#include "source_lcao/module_lr/utils/lr_util_hcontainer.h"

#include <algorithm>
#include <cassert>
#include <complex>

namespace LR
{

template <typename T>
void Core<T>::read_wfc(const std::string& in_dir, const Input_para& inp, const K_Vectors& kv, std::ostream& log)
{
    assert(this->gs.psi != nullptr);
    GsWfc<T> wfc = read_gs_wfc<T>(*this->gs.psi, this->gs.eig, this->gs.wg, this->gs.para_mat, this->gs.para_mat_all,
                                  in_dir, inp, this->dim.nocc_in, this->dim.nvirt_in, this->dim.nocc_max,
                                  this->dim.nbands, this->dim.nbasis, kv, inp.cal_force, log);
    this->gs.psi_all_own = std::move(wfc.all_own);
    this->gs.psi_all = wfc.all;
    this->gs.eig_all = std::move(wfc.eig_all);
    this->gs.wg_all = std::move(wfc.wg_all);
}

template <typename T>
void Core<T>::fill_z(const int* desc_src, const K_Vectors& kv, const bool cal_force, const bool first,
                     std::ostream& log)
{
    ModuleBase::TITLE("LR::Core", "fill_z");
    ZWin<T> z = build_z_win<T>(this->gs.psi_all, this->gs.eig_all, this->dim.nocc, this->dim.nvirt, this->dim.nocc_max,
                               this->dim.nbands, this->dim.nbasis, this->dim.nk, this->dim.openshell, kv,
                               this->amp.para_c, this->gs.para_mat, desc_src, cal_force, first, log);
    this->zw.psi = std::move(z.psi);
    this->zw.eig = std::move(z.eig);
    this->zw.para_c = std::move(z.para_c);
    this->zw.para_x = std::move(z.para_x);
    this->zw.nvirt = std::move(z.nvirt);
    this->zw.nbands = z.nbands;
    this->zw.nloc = z.nloc;
}

template <typename T>
void Core<T>::refresh_gs(const KsBands<T>& bands, const K_Vectors& kv, const Input_para& inp, const bool first,
                         std::ostream& log)
{
    // The ground-state solver owns `psi` and reuses it across ionic steps, so it cannot be stolen.
    // `eig_all` / `wg_all` are nspin x nbands matrices -- a few kB, copied rather than aliased so
    // that they survive the KS solver overwriting its `pelec` on the next step.
    this->gs.psi_all = bands.psi;
    this->gs.eig_all = *bands.ekb;
    this->gs.wg_all = *bands.wg;
    const int max_occupied_bands = *std::max_element(this->dim.nocc.begin(), this->dim.nocc.end());
    const int start_band = this->dim.nocc_max - max_occupied_bands;

    for (int ik = 0; ik < kv.get_nks(); ++ik)
    {
        // copy the KS orbitals in the [nocc, nvirt] window
#ifdef __MPI
        this->gs.psi_all->fix_k(ik);
        this->gs.psi->fix_k(ik);
        const int first_band = start_band + 1;
        Cpxgemr2d(this->dim.nbasis, this->dim.nbands, this->gs.psi_all->get_pointer(), 1, first_band,
                  const_cast<int*>(bands.pv->desc_wfc), this->gs.psi->get_pointer(), 1, 1,
                  const_cast<int*>(this->amp.para_c.desc()), this->amp.para_c.ictxt());
#else
        // serial: each band is `nbasis` contiguous coefficients, so the window is a plain
        // band-by-band copy
        for (int ib = 0; ib < this->dim.nbands; ++ib)
        {
            const auto* start = &(*this->gs.psi_all)(ik, start_band + ib, 0);
            auto* to = &(*this->gs.psi)(ik, ib, 0);
            std::copy(start, start + this->dim.nbasis, to);
        }
#endif
        // copy the KS bands in the [nocc, nvirt] window
        for (int ib = 0; ib < this->dim.nbands; ++ib)
        {
            this->gs.eig(ik, ib) = this->gs.eig_all(ik, start_band + ib);
        }
    }

    if (this->dim.nspin == 2)
    {
        const int nupdown_now = LR_Util::spin_population_difference(*bands.wg, kv.para_k.kpar, kv.para_k.rank_in_pool);
        if (first)
        {
            this->dim.nupdown = nupdown_now;
            apply_spin_pol(this->dim, inp); // shifts nocc/nvirt between the spin channels: exactly once
        }
        else if (nupdown_now != this->dim.nupdown)
        {
            ModuleBase::WARNING_QUIT(
                "LR::Core::refresh_gs",
                "the ground-state spin population changed between ionic steps, but nocc/nvirt and"
                " the distributions built from them were fixed at the first step.");
        }
    }

    // the Z-vector window: after `apply_spin_pol`, so nocc/nvirt/openshell are final
#ifdef __MPI
    this->fill_z(bands.pv->desc_wfc, kv, inp.cal_force, first, log);
#else
    this->fill_z(nullptr, kv, inp.cal_force, first, log);
#endif
}

template <typename T>
void Core<T>::setup_ks_layout(const Parallel_Orbitals& pv, const K_Vectors& kv, const UnitCell& ucell,
                              const Input_para& inp)
{
    // KS and LR must use the same AO distribution, including its BLACS grid.
    const int ks_block_size = pv.para_2d().nb();
    LR_Util::setup_2d_division(this->gs.para_mat.para_2d(), ks_block_size, this->dim.nbasis, this->dim.nbasis
#ifdef __MPI
                               ,
                               pv.para_2d().ictxt()
#endif
    );
    LR::set_band_layout(this->gs.para_mat, this->dim.nbasis, this->dim.nbands, ucell,
                        inp.ri_hartree_benchmark == "aims");
    if (inp.cal_force)
    {
        LR_Util::setup_2d_division(this->gs.para_mat_all.para_2d(), ks_block_size, this->dim.nbasis, this->dim.nbasis
#ifdef __MPI
                                   ,
                                   pv.para_2d().ictxt()
#endif
        );
        LR::set_band_layout(this->gs.para_mat_all, this->dim.nbasis, inp.nbands, ucell,
                            inp.ri_hartree_benchmark == "aims");
    }

    this->gs.para_mat.atom_begin_row = pv.atom_begin_row;
    this->gs.para_mat.atom_begin_col = pv.atom_begin_col;
    this->gs.para_mat.iat2iwt_ = ucell.get_iat2iwt();

    LR_Util::setup_2d_division(this->amp.para_c, ks_block_size, this->dim.nbasis, this->dim.nbands
#ifdef __MPI
                               ,
                               this->gs.para_mat.para_2d().ictxt()
#endif
    );

    // allocate psi_ks and eig_ks in the [nocc, nvirt] window
#ifdef __MPI
    this->gs.psi.reset(
        new psi::Psi<T>(kv.get_nks(), this->amp.para_c.ncloc(), this->amp.para_c.nrloc(), kv.ngk, true));
#else
    this->gs.psi.reset(new psi::Psi<T>(kv.get_nks(), this->dim.nbands, this->dim.nbasis, kv.ngk, true));
#endif
    this->gs.eig.create(kv.get_nks(), this->dim.nbands);
}

template <typename T>
void Core<T>::setup_file_layout(const K_Vectors& kv, const UnitCell& ucell, const Input_para& inp)
{
    LR_Util::setup_2d_division(this->gs.para_mat.para_2d(), 1, this->dim.nbasis, this->dim.nbasis);
    LR::set_band_layout(this->gs.para_mat, this->dim.nbasis, this->dim.nbands, ucell,
                        inp.ri_hartree_benchmark == "aims");
    if (inp.cal_force)
    {
        LR_Util::setup_2d_division(this->gs.para_mat_all.para_2d(), 1, this->dim.nbasis, this->dim.nbasis);
        LR::set_band_layout(this->gs.para_mat_all, this->dim.nbasis, inp.nbands, ucell,
                            inp.ri_hartree_benchmark == "aims");
    }

    // the ground state is read from file, so the band window spans the stored bands directly
    this->gs.psi.reset(new psi::Psi<T>(kv.get_nks(), this->gs.para_mat.ncol_bands,
                                       this->gs.para_mat.para_2d().nrloc(), kv.ngk, true));

    LR_Util::setup_2d_division(this->amp.para_c, 1, this->dim.nbasis, this->dim.nbands
#ifdef __MPI
                               ,
                               this->gs.para_mat.para_2d().ictxt()
#endif
    );
}

template <typename T>
GradientInputs<T> Core<T>::gradient_inputs() const
{
    return {*this->refs.ucell, *this->refs.kv, *this->refs.gd, *this->refs.orb_cutoff, this->gs.para_mat,
            this->zw.para_c, this->zw.para_x, *this->zw.psi, this->zw.eig, this->dim.nocc, this->zw.nvirt,
            this->dim.nspin, this->dim.nk, this->dim.nbasis, this->zw.nloc, *this->refs.xc_kernel,
            this->refs.inp->dft_functional, this->refs.inp->ks_solver, this->refs.inp->test_force,
            this->relax.target.active(), *this->refs.out_dir, this->refs.my_rank, this->amp.spin_types, this->pot.pot,
            this->pot.pot_hxc_gs, *this->refs.log
#ifdef __EXX
            ,
            *this->refs.exx_lri, this->refs.hybrid_alpha
#endif
    };
}

template <typename T>
GsForceCtx<T> Core<T>::gs_force_ctx() const
{
    const LR::ForceSettings settings{this->dim.nspin, this->refs.inp->test_force, this->refs.inp->vl_in_h,
                                     this->refs.inp->vh_in_h, this->refs.inp->dft_functional};
    return {*this->gs.psi_all, this->gs.wg_all, *this->refs.pv_all, *this->refs.rho, *this->refs.rhod,
            *this->refs.vloc, *this->refs.sfac, *this->refs.tcb, settings, this->dim.openshell,
            this->refs.inp->lr_grad_solver, *this->refs.in_dir};
}

template <typename T>
StateCtx<T> Core<T>::state_ctx(const ModuleBase::matrix& ekb) const
{
    return {this->amp.x, this->amp.para_x, this->zw.para_x, this->dim.nocc, this->dim.nvirt, this->dim.nk,
            this->amp.nloc, this->zw.nloc, this->dim.nstates, this->dim.openshell, ekb.c};
}

template <typename T>
std::vector<ModuleBase::matrix> Core<T>::cal_force(const int ispin, const ModuleBase::matrix& ekb,
                                                   const int istate_only)
{
    if (this->refs.inp->test_force && ispin == 0)
    {
        this->test_force();
    }
    return this->dim.openshell
               ? LR::state_forces_os<T>(this->gradient_inputs(), this->gs_force_ctx(), this->state_ctx(ekb),
                                        istate_only)
               : LR::state_forces<T>(this->gradient_inputs(), this->gs_force_ctx(), this->state_ctx(ekb), ispin,
                                     istate_only);
}

template <typename T>
std::vector<ModuleBase::matrix> Core<T>::cal_force_xz(const int ispin, const ct::Tensor& xz,
                                                      const std::vector<double>& omega, const int label_begin)
{
    ModuleBase::timer::start("LR::Core", "cal_force_xz");
    const auto forces = LR::force_xz<T>(this->gradient_inputs(), this->gs_force_ctx(), ispin, xz, omega, label_begin);
    ModuleBase::timer::end("LR::Core", "cal_force_xz");
    return forces;
}

template <typename T>
std::vector<ModuleBase::matrix> Core<T>::cal_force_os(const ct::Tensor& xz, const std::vector<double>& omega,
                                                      const int label_begin)
{
    ModuleBase::timer::start("LR::Core", "cal_force_os");
    const auto forces = LR::force_xz_os<T>(this->gradient_inputs(), this->gs_force_ctx(), xz, omega, label_begin);
    ModuleBase::timer::end("LR::Core", "cal_force_os");
    return forces;
}

template <typename T>
void Core<T>::cal_force_grad(const int ispin, const ModuleBase::matrix& ekb, std::ofstream& ofs)
{
    const std::vector<ModuleBase::matrix> forces = this->cal_force(ispin, ekb, -1);
    if (this->refs.inp->lr_degen_thr <= 0.0)
    {
        return;
    }
    // Rebuilt from scratch for this channel: a relaxation calls this once per ionic step, and a
    // stale multiplet from the previous geometry must not survive into the next one.
    this->relax.mpl.erase(
        std::remove_if(
            this->relax.mpl.begin(), this->relax.mpl.end(), [ispin](const LR::MplLvc& m) { return m.ispin == ispin; }),
        this->relax.mpl.end());
    // The per-state gradients above are the diagonal of the degenerate-subspace gradient matrix in
    // whatever basis the eigensolver returned, so inside a multiplet only their trace means
    // anything. `lr_degen_thr` asks for the off-diagonal part too.
    for (LR::MplLvc& m : LR::multiplet_grad<T>(this->gradient_inputs(), this->gs_force_ctx(), this->state_ctx(ekb),
                                               ispin, this->refs.inp->lr_degen_thr, forces, ofs))
    {
        this->relax.mpl.push_back(std::move(m));
    }
}

template <typename T>
module_dm::DensityMatrix<T, double> Core<T>::cal_dm_gs() const
{
    return LR::make_gs_dm<T>(*this->gs.psi_all, this->gs.wg_all, *this->refs.pv_all, this->gs.para_mat,
                             *this->refs.ucell, *this->refs.gd, *this->refs.orb_cutoff, this->dim.nspin,
                             this->refs.kv->kvec_d, this->dim.nk);
}

template <typename T>
void Core<T>::test_force()
{
    const LR::ForceSettings force_settings{this->dim.nspin, this->refs.inp->test_force, this->refs.inp->vl_in_h,
                                           this->refs.inp->vh_in_h, this->refs.inp->dft_functional};
    LR_Force<T> lr_force((*this->refs.ucell), this->refs.kv->kvec_d, this->gs.para_mat, *this->refs.rhod,
                         *this->refs.rho, *this->refs.vloc, *this->refs.sfac, *this->refs.gd, *this->refs.tcb,
                         force_settings, *this->refs.log
#ifdef __EXX
                         ,
                         std::weak_ptr<Exx_LRI<T>>(*this->refs.exx_lri), this->refs.hybrid_alpha
#endif
    );

    const module_dm::DensityMatrix<T, double>& dm_gs = this->cal_dm_gs();
    ///========================== test 1: reproduce the force of ground state =========================
    // energy density matrix of the ground state
    module_dm::DensityMatrix<T, double> edm_gs(&this->gs.para_mat, this->dim.nspin, this->refs.kv->kvec_d,
                                               this->dim.nk); // DX
    ModuleBase::matrix wg_ekb_ks_all(this->dim.nspin, this->refs.inp->nbands);
    const int occupation_count = this->dim.nspin * this->refs.inp->nbands;
    for (int index = 0; index < occupation_count; ++index)
    {
        wg_ekb_ks_all.c[index] = this->gs.wg_all.c[index] * this->gs.eig_all.c[index];
    }
    module_dm::dm_from_psi(this->refs.pv_all, wg_ekb_ks_all, *this->gs.psi_all, edm_gs);
    LR_Util::initialize_DMR(edm_gs, this->gs.para_mat, (*this->refs.ucell), *this->refs.gd, *this->refs.orb_cutoff);
    edm_gs.cal_dmr(-1);
    // ground-state force
    ModuleBase::matrix force_gs = lr_force.reproduce_force_gs(*this->refs.kv, dm_gs, edm_gs);
    ModuleIO::print_force(*this->refs.log, (*this->refs.ucell), "Ground State FORCE (eV/Angstrom)", force_gs, false);
    /// ======================================= END test 1 =========================================
    ///========================== test 2: reproduce the DX Hartree term =========================
    ModuleBase::matrix f_hxc_potgs = lr_force.reproduce_force_gs_loc(dm_gs, *this->pot.gs_pot.pot_hartree());
    ModuleIO::print_force(*this->refs.log, (*this->refs.ucell),
                          "GS Hartree force calculated by 'cal_pulay_fs' from potential (eV/Angstrom)", f_hxc_potgs,
                          false);
    ModuleBase::matrix f_hxc_potlr = lr_force.cal_force_hxc_dmtrans(dm_gs, *this->pot.pot[0]);
    ModuleIO::print_force(*this->refs.log, (*this->refs.ucell),
                          "GS Hxc force calculated by 'LR_Force' from kernel (eV/Angstrom)", f_hxc_potlr, false);
    /// ======================================= END test 2 =========================================
}

template <typename T>
void Core<T>::setup_relax()
{
    this->relax.target.setup(*this->refs.inp, this->dim.openshell, this->dim.nspin, this->dim.nstates, *this->refs.log);
    if (this->relax.target.active())
    {
        this->relax.force_gs.create(this->refs.ucell->nat, 3);
        this->relax.lr_force.create(this->refs.ucell->nat, 3);
    }
}

template <typename T>
void Core<T>::follow_target(const ModuleBase::matrix& ekb, std::ofstream& ofs)
{
    ModuleBase::timer::start("LR::Core", "follow_target");
    if (!this->relax.target.active())
    {
        ModuleBase::timer::end("LR::Core", "follow_target");
        return;
    }
    const int channel = this->dim.openshell ? 0 : this->relax.target.ispin();
    const T* const amplitudes = this->amp.x[channel].template data<T>();
    const TwoCenterIntegrator& overlap_integrator = *this->refs.tcb->overlap_orb;
    std::vector<double> energies(this->dim.nstates);
    std::copy_n(ekb.c + this->relax.target.ekb_off(this->dim.nstates), this->dim.nstates, energies.begin());
    this->relax.target.follow(*this->refs.inp, *this->refs.ucell, *this->refs.orb_cutoff, overlap_integrator,
                              *this->gs.psi, this->amp.para_c, this->gs.para_mat, this->amp.para_x, this->dim.nocc,
                              this->dim.nvirt, amplitudes, this->amp.nloc, this->dim.nstates, energies, ofs);
    ModuleBase::timer::end("LR::Core", "follow_target");
}

template <typename T>
ct::Tensor Core<T>::pad_group_z(const int ispin, const std::vector<int>& group) const
{
    return LR::pad_group<T>(this->amp.x, this->dim.openshell, this->amp.para_x, this->zw.para_x, this->dim.nocc,
                            this->dim.nvirt, this->dim.nk, this->amp.nloc, this->zw.nloc, ispin, group);
}

template <typename T>
ModuleBase::matrix Core<T>::cal_lr_force(const ModuleBase::matrix& ekb, std::ofstream& ofs)
{
    this->relax.target.resolve_group(*this->refs.inp, ekb.c, this->dim.nstates, *this->refs.log);
    const std::vector<int>& group = this->relax.target.group();
    if (group.empty())
    {
        return this->cal_force(this->relax.target.ispin(), ekb, this->relax.target.state())[0];
    }
    const bool jt_mode = (LR_Util::tolower(this->refs.inp->lr_degen_mode) == "jt");
    // `average` needs only the DIAGONAL of the gradient matrix: the average is basis-independent by
    // construction, so the off-diagonal part is not involved. The Jahn-Teller direction is
    // orthogonal to the average and does need it.
    const int d = static_cast<int>(group.size());
    const int ekb_off = this->relax.target.ekb_off(this->dim.nstates);
    const ct::Tensor xz = this->pad_group_z(this->relax.target.ispin(), group);
    std::vector<double> omega(d);
    for (int k = 0; k < d; ++k)
    {
        omega[k] = ekb.c[ekb_off + group[k]];
    }
    const std::vector<ModuleBase::matrix> forces =
        this->dim.openshell ? this->cal_force_os(xz, omega, group.front())
                            : this->cal_force_xz(this->relax.target.ispin(), xz, omega, group.front());
    ofs << " Followed state " << this->relax.target.state() << " is degenerate with";
    for (const int ist : group)
    {
        if (ist != this->relax.target.state())
        {
            ofs << " " << ist;
        }
    }
    ofs << " (Omega_bar = " << this->relax.target.omega(ekb.c, this->dim.nstates) << " Ry)." << std::endl;
    if (!jt_mode)
    {
        ofs << " lr_degen_mode=average: following the multiplet average, which keeps the"
               " geometry on the current average surface; split groups are not recombined."
            << std::endl;
        return average_forces(forces);
    }
    return this->cal_jt_force(forces, ekb, ofs);
}

template <typename T>
ModuleBase::matrix Core<T>::cal_jt_force(const std::vector<ModuleBase::matrix>& diag, const ModuleBase::matrix& ekb,
                                         std::ofstream& ofs)
{
    // Regime (c): descend the Jahn-Teller branch. This needs the whole gradient matrix, so the
    // off-diagonal elements are assembled here.
    const std::vector<int>& group = this->relax.target.group();
    const LR::MplLvc mpl = LR::build_mpl_state<T>(this->gradient_inputs(), this->gs_force_ctx(),
                                                  this->state_ctx(ekb), this->relax.target.ispin(), group, diag, ofs);
    LR::JtForce jt = LR::eval_jt_force(mpl.g, ofs);
    const int channel = this->dim.openshell ? 0 : this->relax.target.ispin();
    this->relax.target.save_mixed(this->amp.x[channel].template data<T>(), this->amp.nloc, jt.mixing);
    // The returned force is the descending branch's own gradient; once a step splits the
    // multiplet `resolve_group` finds no group and the single-state path takes over.
    return jt.force;
}

template void Core<double>::read_wfc(const std::string&, const Input_para&, const K_Vectors&, std::ostream&);
template void Core<std::complex<double>>::read_wfc(const std::string&, const Input_para&, const K_Vectors&,
                                                   std::ostream&);
template void Core<double>::fill_z(const int*, const K_Vectors&, bool, bool, std::ostream&);
template void Core<std::complex<double>>::fill_z(const int*, const K_Vectors&, bool, bool, std::ostream&);
template void Core<double>::refresh_gs(const KsBands<double>&, const K_Vectors&, const Input_para&, bool,
                                       std::ostream&);
template void Core<std::complex<double>>::refresh_gs(const KsBands<std::complex<double>>&, const K_Vectors&,
                                                     const Input_para&, bool, std::ostream&);
template void Core<double>::setup_ks_layout(const Parallel_Orbitals&, const K_Vectors&, const UnitCell&,
                                            const Input_para&);
template void Core<std::complex<double>>::setup_ks_layout(const Parallel_Orbitals&, const K_Vectors&, const UnitCell&,
                                                          const Input_para&);
template void Core<double>::setup_file_layout(const K_Vectors&, const UnitCell&, const Input_para&);
template void Core<std::complex<double>>::setup_file_layout(const K_Vectors&, const UnitCell&, const Input_para&);

template std::vector<ModuleBase::matrix> Core<double>::cal_force(int, const ModuleBase::matrix&, int);
template std::vector<ModuleBase::matrix> Core<std::complex<double>>::cal_force(int, const ModuleBase::matrix&, int);
template std::vector<ModuleBase::matrix> Core<double>::cal_force_xz(int, const ct::Tensor&, const std::vector<double>&,
                                                                    int);
template std::vector<ModuleBase::matrix> Core<std::complex<double>>::cal_force_xz(int, const ct::Tensor&,
                                                                                  const std::vector<double>&, int);
template std::vector<ModuleBase::matrix> Core<double>::cal_force_os(const ct::Tensor&, const std::vector<double>&, int);
template std::vector<ModuleBase::matrix> Core<std::complex<double>>::cal_force_os(const ct::Tensor&,
                                                                                  const std::vector<double>&, int);
template void Core<double>::cal_force_grad(int, const ModuleBase::matrix&, std::ofstream&);
template void Core<std::complex<double>>::cal_force_grad(int, const ModuleBase::matrix&, std::ofstream&);
template module_dm::DensityMatrix<double, double> Core<double>::cal_dm_gs() const;
template module_dm::DensityMatrix<std::complex<double>, double> Core<std::complex<double>>::cal_dm_gs() const;
template void Core<double>::test_force();
template void Core<std::complex<double>>::test_force();
template void Core<double>::setup_relax();
template void Core<std::complex<double>>::setup_relax();
template void Core<double>::follow_target(const ModuleBase::matrix&, std::ofstream&);
template void Core<std::complex<double>>::follow_target(const ModuleBase::matrix&, std::ofstream&);
template ModuleBase::matrix Core<double>::cal_lr_force(const ModuleBase::matrix&, std::ofstream&);
template ModuleBase::matrix Core<std::complex<double>>::cal_lr_force(const ModuleBase::matrix&, std::ofstream&);
template ModuleBase::matrix Core<double>::cal_jt_force(const std::vector<ModuleBase::matrix>&,
                                                       const ModuleBase::matrix&, std::ofstream&);
template ModuleBase::matrix Core<std::complex<double>>::cal_jt_force(const std::vector<ModuleBase::matrix>&,
                                                                     const ModuleBase::matrix&, std::ofstream&);
template ct::Tensor Core<double>::pad_group_z(int, const std::vector<int>&) const;
template ct::Tensor Core<std::complex<double>>::pad_group_z(int, const std::vector<int>&) const;

} // namespace LR
