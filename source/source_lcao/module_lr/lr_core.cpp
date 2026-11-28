#include "source_lcao/module_lr/lr_core.h"

#include "source_base/module_external/scalapack_connector.h"
#include "source_base/parallel_grid.h"
#include "source_base/timer.h"
#include "source_base/tool_quit.h"
#include "source_base/tool_title.h"
#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_basis/module_nao/two_center_bundle.h"
#include "source_basis/module_nao/two_center_integrator.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_cell/klist.h"
#include "source_cell/module_neighbor/sltk_grid_driver.h"
#include "source_cell/unitcell.h"
#include "source_estate/module_dm/dm_from_psi.h"
#include "source_io/module_output/output_log.h"
#include "source_io/module_parameter/input_parameter.h"
#include "source_lcao/module_lr/gradient_inputs.h"
#include "source_lcao/module_lr/gradient_output.h"
#include "source_lcao/module_lr/hamilt_casida.h"
#include "source_lcao/module_lr/hamilt_ulr.h"
#include "source_lcao/module_lr/hsolver_lrtd.h"
#include "source_lcao/module_lr/kernel_pots.h"
#include "source_lcao/module_lr/lr_amp.h"
#include "source_lcao/module_lr/lr_density.h"
#include "source_lcao/module_lr/lr_dm.h"
#include "source_lcao/module_lr/lr_force.h"
#include "source_lcao/module_lr/lr_gs.h"
#include "source_lcao/module_lr/lr_sforce.h"
#include "source_lcao/module_lr/lr_spectrum.h"
#include "source_lcao/module_lr/lr_zeq.h"
#include "source_lcao/module_lr/operator_casida/operator_lr_diag.h"
#include "source_lcao/module_lr/utils/lr_util.h"
#include "source_lcao/module_lr/utils/lr_util_hcontainer.h"
#include "source_lcao/module_lr/utils/lr_util_print.h"
#include "source_lcao/module_lr/utils/spectrum_weights.h"

#include <algorithm>
#include <cassert>
#include <complex>
#include <iostream>

namespace LR
{

template <typename T>
void Core<T>::read_wfc()
{
    assert(this->gs.psi != nullptr);
    const Input_para& inp = *this->refs.inp;
    GsWfc<T> wfc = read_gs_wfc<T>(*this->gs.psi, this->gs.eig, this->gs.wg, this->gs.para_mat, this->gs.para_mat_all,
                                  *this->refs.in_dir, inp, this->dim.nocc_in, this->dim.nvirt_in, this->dim.nocc_max,
                                  this->dim.nbands, this->dim.nbasis, *this->refs.kv, inp.cal_force, *this->refs.log);
    this->gs.psi_all_own = std::move(wfc.all_own);
    this->gs.psi_all = wfc.all;
    this->gs.eig_all = std::move(wfc.eig_all);
    this->gs.wg_all = std::move(wfc.wg_all);
}

template <typename T>
void Core<T>::read_chg(Charge& chg)
{
    const Input_para& inp = *this->refs.inp;
    LR::read_gs_chg(*this->refs.in_dir, inp, this->dim.nspin, *this->refs.ucell, *this->refs.pgrid, this->refs.rho, chg,
                    this->refs.my_rank, *this->refs.log);
}

template <typename T>
void Core<T>::fill_z_win(const int* desc_src, const bool first)
{
    this->fill_z(desc_src, *this->refs.kv, this->refs.inp->cal_force, first, *this->refs.log);
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
void Core<T>::solve(double* ekb, const int istep)
{
    const Input_para& inp = *this->refs.inp;
    // set once here rather than inside the solver branches: the open-shell `spectrum` branch used
    // to leave it empty, and both the output path and the gradient index it
    this->amp.spin_types =
        this->dim.openshell ? std::vector<std::string>({"updown"}) : std::vector<std::string>({"singlet", "triplet"});
    if (this->dim.nspin == 1)
    {
        this->amp.spin_types = {inp.lr_target_spin};
    }

    // a relaxation writes these once per ionic step; without the suffix every step would overwrite
    // the last, and the trajectory would be impossible to inspect afterwards
    const std::string step_suffix = this->relax.target.active() ? "_step" + std::to_string(istep) : "";
    auto efile_out = [this, &step_suffix](const std::string& label) -> std::string
    { return *this->refs.out_dir + "Excitation_Energy_" + label + step_suffix + ".dat"; };
    auto vfile_out = [this, &step_suffix](const std::string& label) -> std::string
    {
        return *this->refs.out_dir + "Excitation_Amplitude_" + label + step_suffix + "_" +
               std::to_string(this->refs.my_rank + 1) + ".dat";
    };

    if (inp.lr_solver == "elpa")
    {
        ModuleBase::WARNING_QUIT("LR::Core::solve", "ESolver_LR doesn't support elpa now.");
    }
    else if (inp.lr_solver != "spectrum")
    {
        auto write_states = [this, &efile_out, &vfile_out](
                                const std::string& label, const Real<T>* e, const T* v, const int& dim, const int& nst,
                                const int& prec = 8) -> void
        {
            if (this->refs.my_rank == 0)
            {
                assert(nst == LR_Util::write_value(efile_out(label), prec, e, nst));
            }
            assert(nst * dim == LR_Util::write_value(vfile_out(label), prec, v, nst, dim));
        };
        std::vector<double> precondition(inp.lr_solver == "lapack" ? 0 : this->amp.nloc, 1.0);
        // allocate and initialize A matrix and density matrix
        if (this->dim.openshell)
        {
            for (int is : {0, 1})
            {
                if (inp.lr_solver != "lapack")
                {
                    const int offset_is = is * this->amp.para_x[0].size();
                    OperatorLRDiag<double> pre_op(
                        this->gs.eig.c + is * this->dim.nk * (this->dim.nocc[0] + this->dim.nvirt[0]), this->amp.para_x[is],
                        this->dim.nk, this->dim.nocc[is], this->dim.nvirt[is], false);
                    pre_op.act(1, offset_is, 1, precondition.data() + offset_is, precondition.data() + offset_is);
                }
            }
            std::cout << "Solving spin-conserving excitation for open-shell system." << std::endl;
            HamiltULR<T> hulr(
                *this->refs.xc_kernel, this->dim.nspin, this->dim.nbasis, this->dim.nocc, this->dim.nvirt,
                *this->refs.ucell, *this->refs.orb_cutoff, *this->refs.gd, *this->gs.psi, this->gs.eig,
#ifdef __EXX
                *this->refs.exx_lri, this->refs.hybrid_alpha,
#endif
                this->pot.pot, *this->refs.kv, this->amp.para_x, this->amp.para_c, this->gs.para_mat);
            LR::HSolver::solve(hulr, this->amp.x[0].template data<T>(), this->amp.nloc, this->dim.nstates,
                               this->dim.nk, this->dim.nocc, this->dim.nvirt, this->amp.para_x, ekb, inp.lr_solver,
                               inp.lr_thr, precondition, true);
            if (inp.out_wfc_lr)
            {
                write_states("openshell", ekb, this->amp.x[0].template data<T>(), this->amp.nloc, this->dim.nstates);
            }
        }
        else
        {
            if (inp.lr_solver != "lapack")
            {
                OperatorLRDiag<double> pre_op(this->gs.eig.c, this->amp.para_x[0], this->dim.nk, this->dim.nocc[0],
                                              this->dim.nvirt[0], false);
                pre_op.act(1, this->amp.nloc, 1, precondition.data(), precondition.data());
            }
            const std::vector<std::string>& spin_types = this->amp.spin_types;
            for (int is = 0; is < this->dim.nspin; ++is)
            {
                std::cout << " Calculating " << spin_types[is] << " excitations" << std::endl;
                HamiltLR<T> hlr(
                    *this->refs.xc_kernel, this->dim.nspin, this->dim.nbasis, this->dim.nocc, this->dim.nvirt,
                    *this->refs.ucell, *this->refs.orb_cutoff, *this->refs.gd, *this->gs.psi, this->gs.eig,
#ifdef __EXX
                    *this->refs.exx_lri, this->refs.hybrid_alpha,
#endif
                    this->pot.pot[is], *this->refs.kv, this->amp.para_x, this->amp.para_c, this->gs.para_mat,
                    spin_types[is], *this->refs.in_dir, *this->refs.out_dir, inp.ri_hartree_benchmark,
                    (inp.ri_hartree_benchmark == "aims" ? inp.aims_nbasis : std::vector<int>({})));
                LR::HSolver::solve(hlr, this->amp.x[is].template data<T>(), this->amp.nloc, this->dim.nstates,
                                   this->dim.nk, this->dim.nocc, this->dim.nvirt, this->amp.para_x,
                                   ekb + is * this->dim.nstates, inp.lr_solver, inp.lr_thr, precondition, true);
                if (inp.out_wfc_lr)
                {
                    write_states(spin_types[is], ekb + is * this->dim.nstates, this->amp.x[is].template data<T>(),
                                 this->amp.nloc, this->dim.nstates);
                }
            }
        }
    }
    else // lr_solver == "spectrum", read the eigenvalues
    {
        auto efile_in = [&](const std::string& label) -> std::string
        { return *this->refs.in_dir + "Excitation_Energy_" + label + ".dat"; };
        auto vfile_in = [this](const std::string& label) -> std::string
        {
            return *this->refs.in_dir + "Excitation_Amplitude_" + label + "_" +
                   std::to_string(this->refs.my_rank + 1) + ".dat";
        };

        auto read_states = [&](const std::string& label, Real<T>* e, T* v, const int& dim, const int& nst) -> void
        {
            if (this->refs.my_rank == 0)
            {
                assert(nst == LR_Util::read_value(efile_in(label), e, nst));
                std::cout << "Rank " << this->refs.my_rank << ": finish reading " << efile_in(label) << std::endl;
            }
#ifdef __MPI
            // in velocity gauge, the eigenvalues may be used to calculate the transition dipole, so we'd better
            // broadcast them
            MPI_Bcast(e, nst, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
            assert(nst * dim == LR_Util::read_value(vfile_in(label), v, nst, dim));
            std::cout << "Rank " << this->refs.my_rank << ": finish reading " << vfile_in(label) << std::endl;
        };
        std::cout << "reading the excitation states from file: \n";
        if (this->dim.openshell)
        {
            read_states("openshell", ekb, this->amp.x[0].template data<T>(), this->amp.nloc, this->dim.nstates);
        }
        else
        {
            const std::vector<std::string>& spin_types = this->amp.spin_types;
            for (int is = 0; is < this->dim.nspin; ++is)
            {
                read_states(spin_types[is], ekb + is * this->dim.nstates, this->amp.x[is].template data<T>(),
                            this->amp.nloc, this->dim.nstates);
            }
        }
    }
}

template <typename T>
void Core<T>::setup_x()
{
    ModuleBase::TITLE("LR::Core", "setup_x");
    XAmp<T> amp = LR::setup_x<T>(this->amp.para_c, this->dim.nocc, this->dim.nvirt, this->dim.npairs, this->gs.eig,
                                 this->dim.nstates, this->dim.nk, this->dim.nspin, this->dim.openshell,
                                 this->refs.inp->lr_solver == "spectrum", *this->refs.log);
    this->amp.x = std::move(amp.x);
    this->amp.para_x = std::move(amp.para_x);
    this->amp.nloc = amp.nloc;
}

template <typename T>
void Core<T>::init_pot(const Charge& chg_gs)
{
    this->pot.pot.assign(this->dim.nspin, nullptr);
    if (this->refs.inp->ri_hartree_benchmark != "none")
    {
        return;
    }
    const Parallel_Grid& grid = *this->refs.pgrid;
    LR::KernelPotentials prepared = LR::prepare_lr_potentials(*this->refs.inp, *this->refs.xc_kernel, *this->refs.rho,
                                                               *this->refs.ucell, chg_gs, grid, this->dim.nspin,
                                                               this->dim.openshell);
    this->pot.pot = std::move(prepared.pot);
    if (this->refs.inp->cal_force)
    {
        ModuleBase::TITLE("LR::Core", "init_pot_gs");
        this->pot.gs_pot.refresh(*this->refs.inp, *this->refs.xc_kernel, *this->refs.ucell, this->refs.rhod,
                                 this->refs.rho, grid, chg_gs, *this->refs.sfac, *this->refs.vloc,
                                 *this->refs.solvent, this->refs.refresh_sfac, this->refs.has_float_data);
        this->pot.pot_hxc_gs = LR::prepare_gs_potential(*this->refs.inp, *this->refs.xc_kernel, *this->refs.rho,
                                                        *this->refs.ucell, chg_gs, grid, this->dim.nspin,
                                                        this->dim.openshell, prepared);
    }
}

template <typename T>
void Core<T>::output(const ModuleBase::matrix& ekb)
{
    const Input_para& inp = *this->refs.inp;
    std::ofstream& ofs = *this->refs.log;
    if (inp.ri_hartree_benchmark != "none")
    {
        return; // no need to calculate the spectrum in the benchmark routine
    }

    // cal electron-hole density
    if (inp.out_chg[0])
    {
        LR_Density<T> lr_density(*this->refs.ucell, *this->refs.kv, *this->refs.gd, *this->gs.psi, *this->refs.orb_cutoff,
                                 *this->refs.pgrid, this->dim.nspin, this->dim.nocc, this->dim.nvirt, this->dim.nbasis,
                                 this->amp.para_x, this->amp.para_c, this->gs.para_mat, this->dim.openshell,
                                 inp.out_chg[1], *this->refs.out_dir);

        if (this->dim.openshell)
        {
            for (int is = 0; is < this->dim.nspin; ++is)
            {
                const T* const amplitudes = this->amp.x[0].template data<T>();
                const std::string density_label = is == 0 ? "up" : "down";
                lr_density.output_eh_density_all_states(amplitudes, is, this->dim.nstates, density_label, ofs);
            }
        }
        else
        {
            for (int is = 0; is < this->amp.x.size(); ++is)
            {
                const T* const amplitudes = this->amp.x[is].template data<T>();
                lr_density.output_eh_density_all_states(amplitudes, is, this->dim.nstates, this->amp.spin_types[is], ofs);
            }
        }
    }

    const std::string gauge = LR_Util::tolower(inp.abs_gauge);
    const bool calculate_spectrum = gauge != "none";
    const std::vector<double> spatial_weights = LR_Util::spatial_k_weights(this->refs.kv->wk, this->dim.nspin);
    // cal spectrum
    std::vector<std::complex<double>> velocity_mo;
    if (gauge == "velocity")
    {
        const int nspin_tmp = inp.nspin == 2 ? 2 : 1;
        velocity_mo = LR_Util::cal_velocity_mo(*this->refs.ucell, *this->refs.gd, *this->refs.tcb, this->gs.para_mat,
                                               this->amp.para_c, *this->refs.kv, *this->gs.psi, this->dim.nk, nspin_tmp,
                                               this->dim.nbasis, this->dim.nocc, this->dim.nvirt);
    }
    std::vector<double> freq;
    if (calculate_spectrum)
    {
        freq.resize(100);
        std::vector<double> abs_wavelen_range({20, 200}); // default range
        if (inp.abs_wavelen_range.size() >= 2 &&
            std::abs(inp.abs_wavelen_range[1] - inp.abs_wavelen_range[0]) > 0.02)
        {
            abs_wavelen_range = inp.abs_wavelen_range;
        }
        double lambda_diff = std::abs(abs_wavelen_range[1] - abs_wavelen_range[0]);
        double lambda_min = std::min(abs_wavelen_range[1], abs_wavelen_range[0]);
        for (int i = 0; i < freq.size(); ++i)
        {
            freq[i] = 91.126664 / (lambda_min + 0.01 * static_cast<double>(i + 1) * lambda_diff);
        }
    }
    for (int is = 0; is < this->amp.x.size(); ++is)
    {
        if (!calculate_spectrum)
        {
            const double* const energies = &ekb.c[is * this->dim.nstates];
            LR_Util::print_excitation_energies(ofs, energies, this->dim.nstates, this->amp.spin_types[is]);
        }
        else
        {
            LR_Spectrum<T> spectrum(
                this->dim.nspin, this->dim.nbasis, this->dim.nocc, this->dim.nvirt, *this->refs.rho, *this->gs.psi,
                *this->refs.ucell, *this->refs.kv, spatial_weights, *this->refs.gd, *this->refs.orb_cutoff,
                *this->refs.tcb, this->amp.para_x, this->amp.para_c, this->gs.para_mat, &ekb.c[is * this->dim.nstates],
                this->gs.eig.c, this->amp.x[is].template data<T>(), this->dim.nstates, this->dim.openshell, gauge,
                this->refs.my_rank, *this->refs.out_dir);
            if (gauge == "velocity")
            {
                spectrum.set_vmo(velocity_mo.data());
            }
            spectrum.cal_spectrum();
            spectrum.transition_analysis(this->amp.spin_types[is] + "_tda");
            if (this->amp.spin_types[is] != "triplet") // triplets has no transition dipole and no contribution to the spectrum
            {
                spectrum.optical_absorption_method1(freq, inp.abs_broadening);
                spectrum.write_transition_dipole(*this->refs.out_dir + "trans_dipole_" + this->amp.spin_types[is] + "_tda.dat");
            }
        }
        if (inp.cal_force && !this->relax.target.active())
        {
            this->cal_force_grad(is, ekb, ofs);
        }
    }
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
    const auto forces = LR::force_xz<T>(this->gradient_inputs(), this->gs_force_ctx(), this->dim.nvirt[0], ispin, xz,
                                        omega, label_begin);
    ModuleBase::timer::end("LR::Core", "cal_force_xz");
    return forces;
}

template <typename T>
std::vector<ModuleBase::matrix> Core<T>::cal_force_os(const ct::Tensor& xz, const std::vector<double>& omega,
                                                      const int label_begin)
{
    ModuleBase::timer::start("LR::Core", "cal_force_os");
    const auto forces =
        LR::force_xz_os<T>(this->gradient_inputs(), this->gs_force_ctx(), this->dim.nvirt, xz, omega, label_begin);
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

template void Core<double>::read_wfc();
template void Core<std::complex<double>>::read_wfc();
template void Core<double>::read_chg(Charge&);
template void Core<std::complex<double>>::read_chg(Charge&);
template void Core<double>::fill_z(const int*, const K_Vectors&, bool, bool, std::ostream&);
template void Core<std::complex<double>>::fill_z(const int*, const K_Vectors&, bool, bool, std::ostream&);
template void Core<double>::fill_z_win(const int*, bool);
template void Core<std::complex<double>>::fill_z_win(const int*, bool);
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
template void Core<double>::solve(double*, int);
template void Core<std::complex<double>>::solve(double*, int);
template void Core<double>::output(const ModuleBase::matrix&);
template void Core<std::complex<double>>::output(const ModuleBase::matrix&);
template void Core<double>::setup_x();
template void Core<std::complex<double>>::setup_x();
template void Core<double>::init_pot(const Charge&);
template void Core<std::complex<double>>::init_pot(const Charge&);

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
