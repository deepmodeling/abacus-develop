// LR ground-state distributions, excitation amplitudes, and file input.
#include "source_lcao/module_lr/lr_gs.h"

#include "source_base/module_external/scalapack_connector.h"
#include "source_base/module_parallel/para_bridge.h"
#include "source_base/tool_quit.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_cell/klist.h"
#include "source_cell/unitcell.h"
#include "source_estate/module_charge/charge.h"
#include "source_estate/module_charge/chg_rhog_io.h"
#include "source_hamilt/module_xc/xc_functional.h"
#include "source_io/module_output/cube_io.h"
#include "source_io/module_parameter/input_parameter.h"
#include "source_io/module_wf/read_wfc_nao.h"
#include "source_lcao/module_lr/ri_benchmark/ri_benchmark.h"

#include <algorithm>
#include <complex>
#include <sstream>

namespace LR
{

void set_band_layout(Parallel_Orbitals& pmat, int nbasis, int nbands, const UnitCell& ucell, bool aims)
{
#ifdef __MPI
    pmat.set_wfc_and_mo_layouts(nbasis, nbands);
    // The KS helper rejects empty band columns. LR windows may legitimately
    // have fewer blocks than the inherited KS process grid has columns.
    Parallel_Blkcyclic bands;
    const int block_size = pmat.para_2d().nb();
    LR_Util::setup_2d_division(bands, block_size, nbasis, nbands, pmat.para_2d().ictxt());
    pmat.nbands = nbands;
    pmat.nrow_bands = bands.nrloc();
    pmat.ncol_bands = bands.ncloc();
    pmat.nloc_wfc = bands.size();
    // Skipped for the aims benchmark: with `aims_nbasis` the per-atom orbital counts behind
    // `iat2iwt` do not match `nbasis`, so the atomic trace would be wrong. The guard came from
    // d4fe3fe84 ("Support different basis number from aims"), and was silently undone by
    // 1e4c1c6af (BSE, #7718) re-adding an unconditional call above it.
    if (!aims)
    {
        pmat.set_atomic_trace(ucell.get_iat2iwt(), ucell.nat, nbasis);
    }
#else
    pmat.nrow_bands = nbasis;
    pmat.ncol_bands = nbands;
#endif
}

AmpLayout make_amp_layout(const Parallel_Blkcyclic& para_c, const std::vector<int>& nocc,
                          const std::vector<int>& nvirt, int nk, bool openshell)
{
    AmpLayout layout;
    const int block_size = para_c.nb();
    layout.para_x.reserve(nocc.size());
    for (int is = 0; is < static_cast<int>(nocc.size()); ++is)
    {
        Parallel_Blkcyclic px;
        LR_Util::setup_2d_division(px, block_size, nvirt[is], nocc[is]
#ifdef __MPI
            ,
            para_c.ictxt()
#endif
        ); // nvirt - row, nocc - col
        layout.para_x.emplace_back(std::move(px));
    }
    layout.nloc =
        nk * (openshell ? layout.para_x[0].size() + layout.para_x[1].size() : layout.para_x[0].size());
    return layout;
}

template <typename T>
void set_x_init(std::vector<ct::Tensor>& X, const ModuleBase::matrix& eig_ks, const std::vector<int>& nocc,
                const std::vector<int>& nvirt, const std::vector<int>& npairs,
                const std::vector<Parallel_Blkcyclic>& para_x, int nstates, int nloc, int nk, bool openshell,
                std::ostream& log)
{
    // set the initial guess of X
    for (int is = 0; is < static_cast<int>(nocc.size()); ++is)
    {
        const int& no = nocc[is];
        const int& nv = nvirt[is];
        const int& np = npairs[is];
        const Parallel_Blkcyclic& px = para_x[is];

        // if (E_{lumo}-E_{homo-1} < E_{lumo+1}-E{homo}), mode = 0, else 1(smaller first)
        bool ix_mode = false; // default
        if (eig_ks.nc > no + 1 && no >= 2
            && eig_ks(is, no) - eig_ks(is, no - 2) - 1e-5 > eig_ks(is, no + 1) - eig_ks(is, no - 1))
        {
            ix_mode = true;
        }
        log << "setting the initial guess of X of spin" << is << std::endl;
        if (no >= 2 && eig_ks.nc > no)
        {
            log << "E_{lumo}-E_{homo-1}=" << eig_ks(is, no) - eig_ks(is, no - 2) << std::endl;
        }
        if (no >= 1 && eig_ks.nc > no + 1)
        {
            log << "E_{lumo+1}-E{homo}=" << eig_ks(is, no + 1) - eig_ks(is, no - 1) << std::endl;
        }
        log << "mode of X-index: " << ix_mode << std::endl;

        /// global index map between (i,c) and ix
        const auto indexmap = LR_Util::set_ix_map_diagonal(ix_mode, no, nv);
        const std::vector<std::pair<int, int>>& ix2ioiv = std::get<1>(indexmap);

        for (int ib = 0; ib < nstates; ++ib)
        {
            const int ipair = ib % np;
            const int occ_global = ix2ioiv[ipair].first;
            const int virt_global = ix2ioiv[ipair].second;
            const int ik = ib / np;
            const int xstart_b = ib * nloc; // start index of band ib
            const int xstart_bs =
                (openshell && is == 1) ? xstart_b + nk * para_x[0].size() : xstart_b; // start index of band ib, spin is
            const int is_in_x = openshell ? 0 : is; // if openshell, spin-up and spin-down are put together
            if (px.in_this_processor(virt_global, occ_global))
            {
                const int xstart_pair = ik * px.size();
                const int ipair_loc = px.icloc(occ_global) * px.nrloc() + px.irloc(virt_global);
                X[is_in_x].data<T>()[xstart_bs + xstart_pair + ipair_loc] = (static_cast<T>(1.0) / static_cast<T>(nk));
            }
        }
    }
}

template <typename T>
void read_wfc_aims(psi::Psi<T>& psi, ModuleBase::matrix& eig, const std::string& in_dir, int nocc_in, int nvirt_in,
                   int nbands, int nbasis)
{
#ifdef __EXX
    int ncore = 0;
    std::vector<double> eig_ks_vec =
        RI_Benchmark::read_aims_ebands<double>(in_dir + "band_out", nocc_in, nvirt_in, ncore);
    std::cout << "ncore=" << ncore << ", nocc=" << nocc_in << ", nvirt=" << nvirt_in << ", nbands=" << nbands
              << std::endl;
    std::cout << "eig_ks_vec.size()=" << eig_ks_vec.size() << std::endl;
    if (eig_ks_vec.size() != static_cast<std::size_t>(nbands))
    {
        ModuleBase::WARNING_QUIT("ESolver_LR", "read_aims_ebands failed.");
    }
    for (int i = 0; i < nbands; ++i)
    {
        eig(0, i) = eig_ks_vec[i];
    }
    RI_Benchmark::read_aims_eigenvectors<T>(psi, in_dir + "KS_eigenvectors.out", ncore, nbands, nbasis);
#else
    ModuleBase::WARNING_QUIT("ESolver_LR", "RI benchmark is only supported when compile with LibRI.");
#endif
}

template <typename T>
bool read_wfc_win(const std::string& in_dir, const Parallel_Orbitals& pmat, psi::Psi<T>& psi, ModuleBase::matrix& eig,
                  ModuleBase::matrix& wg, const K_Vectors& kv, int nspin, bool binary, int skip_bands)
{
    return ModuleIO::read_wfc_nao<T>(in_dir, pmat, psi, eig, wg, kv.ik2iktot, kv.get_nkstot(), nspin, binary,
                                     skip_bands);
}

template <typename T>
bool read_wfc_all(const std::string& in_dir, const Parallel_Orbitals& pmat, psi::Psi<T>& psi, ModuleBase::matrix& eig,
                  ModuleBase::matrix& wg, const K_Vectors& kv, int nspin, bool binary)
{
    return ModuleIO::read_wfc_nao<T>(in_dir, pmat, psi, eig, wg, kv.ik2iktot, kv.get_nkstot(), nspin, binary,
                                     /*skip_bands=*/0);
}

template <typename T>
ZWin<T> build_z_win(psi::Psi<T>* psi_all, const ModuleBase::matrix& eig_all, const std::vector<int>& nocc,
                    const std::vector<int>& nvirt, int nocc_max, int nbands_x, int nbasis, int nk, bool openshell,
                    const K_Vectors& kv, const Parallel_Blkcyclic& para_c, const Parallel_Orbitals& pmat,
                    const int* desc_src, bool need, bool diagnose, std::ostream& log)
{
    ZWin<T> z;
    if (!need || psi_all == nullptr)
    {
        return z;
    }

    const int max_occupied_bands = *std::max_element(nocc.begin(), nocc.end());
    const int start_band = nocc_max - max_occupied_bands;
    // Every band the ground state solved, from the window start upward. `eig_all` is the
    // authority on how many there are: on the ks-lr path it is the KS solver's `ekb`, on the
    // file path it was read with `skip_bands = 0`.
    z.nbands = eig_all.nc - start_band;
    if (z.nbands <= nbands_x)
    { // nothing to widen: nbands was already at (or below) the X window
        z.nbands = nbands_x;
    }
    const int nspin = static_cast<int>(nocc.size());
    z.nvirt.assign(nspin, 0);
    for (int is = 0; is < nspin; ++is)
    {
        z.nvirt[is] = z.nbands - nocc[is];
    }

    const int block_size = para_c.nb();
    LR_Util::setup_2d_division(z.para_c, block_size, nbasis, z.nbands
#ifdef __MPI
        ,
        pmat.para_2d().ictxt()
#endif
    );
    z.para_x.clear();
    for (int is = 0; is < nspin; ++is)
    {
        Parallel_Blkcyclic px;
        LR_Util::setup_2d_division(px, block_size, z.nvirt[is], nocc[is]
#ifdef __MPI
            ,
            z.para_c.ictxt()
#endif
        );
        z.para_x.emplace_back(std::move(px));
    }
    // Open shell solves one eigenproblem whose vector is the concatenation [up | down], so a
    // state's block is as long as both channels together; the closed-shell singlet/triplet
    // algorithm carries a single channel. (Same rule as `nloc` in `make_amp_layout`, applied to
    // the widened windows.)
    z.nloc = nk * (openshell ? z.para_x[0].size() + z.para_x[1].size() : z.para_x[0].size());

#ifdef __MPI
    z.psi.reset(new psi::Psi<T>(kv.get_nks(), z.para_c.ncloc(), z.para_c.nrloc(), kv.ngk, true));
#else
    z.psi.reset(new psi::Psi<T>(kv.get_nks(), z.nbands, nbasis, kv.ngk, true));
#endif
    z.eig.create(kv.get_nks(), z.nbands);

    for (int ik = 0; ik < kv.get_nks(); ++ik)
    {
        // same redistribution the X window does, over more bands
#ifdef __MPI
        psi_all->fix_k(ik);
        z.psi->fix_k(ik);
        const int first_band = start_band + 1;
        Cpxgemr2d(nbasis, z.nbands, psi_all->get_pointer(), 1, first_band, const_cast<int*>(desc_src),
                  z.psi->get_pointer(), 1, 1, const_cast<int*>(z.para_c.desc()), z.para_c.ictxt());
#else
        for (int ib = 0; ib < z.nbands; ++ib)
        {
            const auto* start = &(*psi_all)(ik, start_band + ib, 0);
            std::copy(start, start + nbasis, &(*z.psi)(ik, ib, 0));
        }
#endif
        for (int ib = 0; ib < z.nbands; ++ib)
        {
            z.eig(ik, ib) = eig_all(ik, start_band + ib);
        }
    }
    log << "Z-vector window: nbands = " << z.nbands << " (X window: " << nbands_x << "), nvirt =";
    for (int is = 0; is < nspin; ++is)
    {
        log << " " << z.nvirt[is] << "(X window: " << nvirt[is] << ")";
    }
    log << std::endl;
    // Initialization precedes ks_initialized_; subsequent geometry refreshes keep
    // the same windows. Diagnose once, without changing the truncated eigenproblem.
    if (diagnose)
    {
        const std::string warning =
            LR_Util::force_window_warning(nbasis, eig_all.nc, start_band, nvirt, z.nvirt);
        if (!warning.empty())
        {
            log << warning << std::endl;
        }
    }
    return z;
}

void read_gs_chg(const std::string& in_dir, const Input_para& inp, int nspin, const UnitCell& ucell,
                 Parallel_Grid& pgrid, ModulePW::PW_Basis* pw_rho, Charge& chg, int rank, std::ostream& log)
{
    chg.set_rhopw(pw_rho);
    const bool kin_den = XC_Functional::get_ked_flag() || (inp.out_elf[0] > 0); // mohan add 20251202
    chg.allocate(nspin, kin_den, XC_Functional::get_ked_flag(), inp.test_charge);
    // Keep the existing cube input authoritative. Choose one density source for
    // all spin channels and broadcast the decision before collective readers.
    const Parallel::ParaWorld world = Parallel::make_pw_world();
    const bool have_cubes = module_charge::has_complete_charge_cubes(in_dir, nspin, world);
    if (!have_cubes)
    {
        const std::string density_file = in_dir + inp.suffix + "-CHARGE-DENSITY.restart";
        const bool loaded = module_charge::read_rho_from_binary(
            density_file, *pw_rho, nspin, chg.rhog, chg.rho, world, log);
        if (!loaded)
        {
            const std::string message = "No complete charge cube set or readable density restart in " + in_dir;
            ModuleBase::WARNING_QUIT("ESolver_LR::read_ks_chg", message);
        }
        return;
    }
    log << " try to read charge from file : ";
    for (int is = 0; is < nspin; ++is)
    {
        std::stringstream ssc;
        if (nspin == 1)
        {
            ssc << in_dir << "chg.cube";
        }
        else
        {
            ssc << in_dir << "chgs" << is + 1 << ".cube";
        }
        log << ssc.str() << std::endl;
        if (ModuleIO::read_vdata_palgrid(pgrid, rank, log, ssc.str(), chg.rho[is], ucell.nat))
        {
            log << " Read in the charge density: " << ssc.str() << std::endl;
        }
    }
}

template void set_x_init<double>(std::vector<ct::Tensor>&, const ModuleBase::matrix&, const std::vector<int>&,
                                 const std::vector<int>&, const std::vector<int>&,
                                 const std::vector<Parallel_Blkcyclic>&, int, int, int, bool, std::ostream&);
template void set_x_init<std::complex<double>>(std::vector<ct::Tensor>&, const ModuleBase::matrix&,
                                               const std::vector<int>&, const std::vector<int>&,
                                               const std::vector<int>&, const std::vector<Parallel_Blkcyclic>&, int,
                                               int, int, bool, std::ostream&);
template void read_wfc_aims<double>(psi::Psi<double>&, ModuleBase::matrix&, const std::string&, int, int, int, int);
template void read_wfc_aims<std::complex<double>>(psi::Psi<std::complex<double>>&, ModuleBase::matrix&,
                                                  const std::string&, int, int, int, int);
template bool read_wfc_win<double>(const std::string&, const Parallel_Orbitals&, psi::Psi<double>&,
                                   ModuleBase::matrix&, ModuleBase::matrix&, const K_Vectors&, int, bool, int);
template bool read_wfc_win<std::complex<double>>(const std::string&, const Parallel_Orbitals&,
                                                 psi::Psi<std::complex<double>>&, ModuleBase::matrix&,
                                                 ModuleBase::matrix&, const K_Vectors&, int, bool, int);
template bool read_wfc_all<double>(const std::string&, const Parallel_Orbitals&, psi::Psi<double>&,
                                   ModuleBase::matrix&, ModuleBase::matrix&, const K_Vectors&, int, bool);
template bool read_wfc_all<std::complex<double>>(const std::string&, const Parallel_Orbitals&,
                                                 psi::Psi<std::complex<double>>&, ModuleBase::matrix&,
                                                 ModuleBase::matrix&, const K_Vectors&, int, bool);
template ZWin<double> build_z_win<double>(psi::Psi<double>*, const ModuleBase::matrix&, const std::vector<int>&,
                                          const std::vector<int>&, int, int, int, int, bool, const K_Vectors&,
                                          const Parallel_Blkcyclic&, const Parallel_Orbitals&, const int*, bool, bool,
                                          std::ostream&);
template ZWin<std::complex<double>> build_z_win<std::complex<double>>(
    psi::Psi<std::complex<double>>*, const ModuleBase::matrix&, const std::vector<int>&, const std::vector<int>&, int,
    int, int, int, bool, const K_Vectors&, const Parallel_Blkcyclic&, const Parallel_Orbitals&, const int*, bool, bool,
    std::ostream&);

} // namespace LR
