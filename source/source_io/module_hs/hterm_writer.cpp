#include "hterm_writer.h"
#include "source_io/module_output/spin_tag.h"
#include "source_io/module_output/write_notice.h"

#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include "source_estate/module_pot/h_hartree_pw.h"
#include "source_hamilt/module_xc/xc_functional.h"
#include "source_io/module_hs/hsr_writer.h"
#include "source_base/module_out/filename.h"
#include "source_cell/ucell_io.h"
#include "source_hamilt/module_gint/gint_interface.h"
#include "source_hamilt/module_hcontainer/hcontainer_funcs.h"
#include "source_hamilt/module_hcontainer/output_hcontainer.h"
#include "source_lcao/module_operator_lcao/nonlocal.h"
#include "source_lcao/module_operator_lcao/operator_fs_utils.h"
#ifdef __EXX
#include "source_lcao/module_operator_lcao/op_exx_lcao.h"
#include "source_lcao/module_ri/exx_lri_interface.h"
#include "source_lcao/module_ri/ri_2d_comm.h"
#endif

#include <cmath>
#include <complex>
#include <fstream>
#include <functional>
#include <iomanip>
#include <ios>
#include <string>
#include <tuple>
#include <vector>

namespace ModuleIO
{

static void setup_veff_hcontainer(hamilt::HContainer<double>& hR,
                                  const UnitCell& ucell,
                                  const Grid_Driver& gd,
                                  const Parallel_Orbitals& pv,
                                  const std::vector<double>& orb_cutoff)
{
    const Parallel_Orbitals* paraV = hR.get_paraV();
    for (int iat1 = 0; iat1 < ucell.nat; iat1++)
    {
        auto tau1 = ucell.get_tau(iat1);
        int T1 = 0, I1 = 0;
        ucell.iat2iait(iat1, &I1, &T1);

        AdjacentAtomInfo adjs;
        gd.Find_atom(ucell, tau1, T1, I1, &adjs);

        for (int ad = 0; ad < adjs.adj_num + 1; ++ad)
        {
            const int T2 = adjs.ntype[ad];
            const int I2 = adjs.natom[ad];
            const int iat2 = ucell.itia2iat(T2, I2);
            if (paraV->is_invalid_atom_pair(iat1, iat2))
            {
                continue;
            }
            const ModuleBase::Vector3<int>& R_index = adjs.box[ad];
            if (ucell.cal_dtau(iat1, iat2, R_index).norm() * ucell.lat0 < orb_cutoff[T1] + orb_cutoff[T2])
            {
                hamilt::AtomPair<double> tmp(iat1, iat2, R_index, paraV);
                hR.insert_pair(tmp);
            }
        }
    }
    hR.allocate(nullptr, true);
}

static void gather_and_write(const std::string& prefix,
                             const std::string& label,
                             const std::string& term_name,
                             hamilt::HContainer<double>& hR,
                             const UnitCell& ucell,
                             const Parallel_Orbitals& pv,
                             const int nspin,
                             const int ispin,
                             const int istep,
                             const int out_app_flag,
                             const int* iat2iwt,
                             const int nat,
                             const std::string& global_matrix_dir,
                             std::ofstream* ofs_running)
{
    // Mode 1 appends one frame per ionic step; the first step (istep == 0)
    // truncates the fixed file so a rerun does not keep stale records. Only
    // mode 2 embeds the step in the name.
    const bool append = (out_app_flag == 1) && (istep > 0);
    const bool numbered = (out_app_flag == 2);
    const int nbasis = hR.get_nbasis();
#ifdef __MPI
    Parallel_Orbitals serialV;
    serialV.para_2d().set_serial(nbasis, nbasis);
    serialV.set_atomic_trace(iat2iwt, nat, nbasis);
    hamilt::HContainer<double> hr_serial(&serialV);
    hamilt::gatherParallels(hR, &hr_serial, 0);
    if (GlobalV::MY_RANK == 0)
#endif
    {
        const std::string fname = global_matrix_dir + hsr_gen_fname(prefix, ispin, numbered, istep);
        // H(R) component terms (T, V^NL, V^L, V^H, V^XC, V^EXX) are parts of
        // the Hamiltonian, not the full H(R): they do not carry the Fermi annotation
        const double no_efermi = 0.0;
#ifdef __MPI
        const bool written = write_hcontainer_csr(fname, &ucell, 8, &hr_serial, istep, append, ispin, nspin, label, "", no_efermi, false);
#else
        const bool written = write_hcontainer_csr(fname, &ucell, 8, &hR, istep, append, ispin, nspin, label, "", no_efermi, false);
#endif
        if (written && ofs_running != nullptr)
        {
            const std::string property = "H(R) (" + term_name + " term)";
            ModuleIO::write_matrix_notice(*ofs_running, property, ispin, nspin, fname);
        }
    }
}

// Write one folded Hamiltonian-term matrix H(k) using the same text layout as
// write_dmk / write_hsk: metadata block, UnitCell block, flat matrix body
// (8 real or 4 complex values per line, no "Row N" markers).
// Returns true only if the record was fully written; on an open or write
// failure the caller must not announce the file in the running log.
static bool write_hk_text(const std::string& fname,
                          const std::vector<std::complex<double>>& hk_global,
                          const K_Vectors& kv,
                          const int ik,
                          const int ibz_ik,
                          const int nlocal,
                          const int precision,
                          const bool gamma_only,
                          const int nspin,
                          const int ispin,
                          const int istep,
                          const UnitCell& ucell,
                          const std::string& term_name,
                          const bool append)
{
    std::ios_base::openmode mode = std::ios::out;
    if (append)
    {
        mode |= std::ios::app;
    }
    std::ofstream ofs(fname.c_str(), mode);
    if (!ofs.is_open())
    {
        ModuleBase::WARNING("ModuleIO::write_h_term", "Cannot open file < " + fname + " >.");
        return false;
    }

    // metadata block, mirrors write_dmk / write_hsk
    ofs << "Ionic_Step " << istep + 1 << std::endl;
    ofs << " " << nspin << " # number of spin directions" << std::endl;
    ofs << " " << ispin + 1 << " # spin index" << std::endl;
    ofs << " " << kv.get_nkstot_nospin() << " # total k points " << std::endl;
    ofs << " " << kv.get_nkstot() << " # total k points after symmetrized (if open) " << std::endl;
    ofs << " " << ibz_ik + 1 << " # k-point index " << std::endl;
    ofs << " " << kv.kvec_c[ik].x << " " << kv.kvec_c[ik].y << " " << kv.kvec_c[ik].z
        << " # k point coordinate (Cartesian) " << std::endl;
    ofs << " " << kv.kvec_d[ik].x << " " << kv.kvec_d[ik].y << " " << kv.kvec_d[ik].z
        << " # k point coordinate (direct) " << std::endl;
    ofs << " " << kv.wk[ik] << " # weight of this k point" << std::endl;
    // A single Hamiltonian term has no Fermi-energy concept; keep a zero
    // placeholder with an explicit label so the column stays parseable.
    ofs << " 0.0 # Fermi energy in Ry (not applicable to " << term_name << ")" << std::endl;
    ofs << " " << nlocal << " # number of localized basis " << std::endl;
    ofs << " " << nlocal << " " << nlocal << " # size of this matrix " << std::endl;
    ofs << std::endl;

    // UnitCell lattice block
    ModuleIO::UcellIO::write_ucell(ofs, &ucell);

    // matrix body: flat, 8 real / 4 complex values per line
    const double dm_thr = 1.0e-16;
    ofs << std::fixed;
    ofs << std::scientific;
    ofs << std::setprecision(precision);
    ofs << std::right;
    for (int i = 0; i < nlocal; ++i)
    {
        const size_t ii = i * nlocal;
        for (int j = 0; j < nlocal; ++j)
        {
            if (gamma_only)
            {
                if (j % 8 == 0)
                {
                    ofs << "\n";
                }
                ofs << " " << hk_global[ii + j].real();
            }
            else
            {
                if (j % 4 == 0)
                {
                    ofs << "\n";
                }
                double real_v = std::real(hk_global[ii + j]);
                if (std::abs(real_v) < dm_thr)
                {
                    real_v = 0.0;
                }
                double imag_v = std::imag(hk_global[ii + j]);
                if (std::abs(imag_v) < dm_thr)
                {
                    imag_v = 0.0;
                }
                ofs << " (" << real_v << "," << imag_v << ")";
            }
        }
    }
    // Terminate the body so that an appended record (out_app_flag == 1) starts
    // on a fresh line instead of being glued to the last value.
    ofs << std::endl;
    // Flush explicitly so that a mid-write failure (e.g. a full disk) is
    // detected here instead of being dropped silently by the destructor.
    ofs.flush();
    return ofs.good();
}

// Canonical running-log name of a Hamiltonian-term H(k) output prefix, so that
// the notice matches the quantity name used by the other matrix writers.
static std::string hk_property_name(const std::string& prefix)
{
    if (prefix == "vnlk")
    {
        return "Vnl(k)";
    }
    if (prefix == "vlk")
    {
        return "Vl(k)";
    }
    if (prefix == "vhk")
    {
        return "Vh(k)";
    }
    if (prefix == "vxck")
    {
        return "Vxc(k)";
    }
    return prefix;
}

// Fold one Hamiltonian term at k point ik and write it in the dmk-compatible
// text layout. The HContainer passed in must already be complete (gathered to
// a serial container) on the writing rank.
static void write_hk_at_k(const hamilt::HContainer<double>& hR,
                          const std::string& prefix,
                          const K_Vectors& kv,
                          const UnitCell& ucell,
                          const int nspin,
                          const int nspin_name,
                          const int ispin,
                          const int istep,
                          const bool append,
                          const int ik,
                          const int nlocal,
                          const int precision,
                          const bool gamma_only,
                          const std::string& global_out_dir,
                          const int out_app_flag,
                          const std::string& term_name,
                          const bool do_write,
                          std::ofstream* ofs_running)
{
    if (!do_write)
    {
        return;
    }

    // IBZ k-point index, matching write_dmk / write_hsk: after symmetry
    // reduction both spin copies of the same IBZ point carry the same label,
    // so subtract half of nkstot for the spin-down block.
    const int half_k = (nspin == 2) ? kv.get_nkstot() / 2 : 0;
    const int ik0 = kv.ik2iktot[ik];
    const int ibz_ik = (nspin == 2 && ik0 >= half_k) ? (ik0 - half_k) : ik0;

    std::vector<std::complex<double>> hk_global(nlocal * nlocal, 0);
    hamilt::folding_HR(hR, hk_global.data(), kv.kvec_d[ik], nlocal, 0);

    const int out_label = 1; // text
    const std::string fname = ModuleIO::filename_output(global_out_dir,
                                                        prefix,
                                                        "nao",
                                                        ik,
                                                        kv.ik2iktot,
                                                        nspin_name,
                                                        kv.get_nkstot(),
                                                        out_label,
                                                        out_app_flag,
                                                        gamma_only,
                                                        istep);

    const bool written = write_hk_text(fname, hk_global, kv, ik, ibz_ik, nlocal, precision, gamma_only,
                                       nspin, ispin, istep, ucell, term_name, append);

    if (written && ofs_running != nullptr)
    {
        // Spin-independent terms are written once with nspin_name == 1, which
        // makes write_matrix_notice emit no spin label for them.
        const std::string property = hk_property_name(prefix);
        ModuleIO::write_matrix_notice(*ofs_running, property, ispin, nspin_name, fname);
    }
}

// Descriptor for one Hamiltonian term: the H(k)/H(R) filename prefixes and the
// term-specific routine that fills hR_tmp for one spin channel.
struct HTermSpec
{
    std::string timer_label;
    std::string k_prefix;
    std::string r_prefix;
    std::string r_label;
    std::string term_name;
    // Fills hR_tmp for spin channel ispin. setup_veff_hcontainer has already
    // been run when use_veff_layout is true; the builder adds the term values.
    std::function<void(hamilt::HContainer<double>&, int)> build;
    // Whether the term differs between the two spin channels when nspin == 2.
    // Spin-independent terms are built and written once and carry no spin
    // suffix in the filename, matching the S(k) convention.
    bool spin_dependent;
};

// Shared skeleton for the five non-EXX terms: loop over spin channels, build
// the term into hR_tmp, always write H(k), and write H(R) only when asked.
static void write_h_term(WriteHParams& params, const HTermSpec& spec)
{
    ModuleBase::TITLE("ModuleIO", spec.timer_label);
    ModuleBase::timer::start("ModuleIO", spec.timer_label);

    const UnitCell& ucell = *params.ucell;
    const Parallel_Orbitals& pv = *params.pv;
    const K_Vectors& kv = *params.kv;
    const int nspin = params.nspin;
    const int istep = params.istep;
    const bool append = params.append;
    const int* iat2iwt = params.iat2iwt;
    const int nat = params.nat;
    const bool also_hR = params.also_hR;
    const int nlocal = params.nlocal;
    const bool gamma_only = params.gamma_only_local;
    const std::string& global_out_dir = params.global_out_dir;
    const int out_app_flag = params.out_app_flag;

    const int nspin_out = (spec.spin_dependent && nspin == 2) ? 2 : 1;
    // Spin-independent terms pass nspin_name = 1 so that filename_output does
    // not append a meaningless spin block; spin-dependent terms keep the real
    // nspin so that the two channels land in s1 and s2 files.
    const int nspin_name = spec.spin_dependent ? nspin : 1;
    const int nks = kv.get_nks();

    for (int ispin = 0; ispin < nspin_out; ispin++)
    {
        hamilt::HContainer<double> hR_tmp(const_cast<Parallel_Orbitals*>(&pv));
        spec.build(hR_tmp, ispin);

#ifdef __MPI
        // Gather the distributed HContainer onto the writing rank so that the
        // folded H(k) is complete; folding the rank-local pieces directly
        // would scatter the basis blocks into wrong matrix positions.
        Parallel_Orbitals serialV;
        serialV.para_2d().set_serial(nlocal, nlocal);
        serialV.set_atomic_trace(iat2iwt, nat, nlocal);
        hamilt::HContainer<double> hR_serial(&serialV);
        hamilt::gatherParallels(hR_tmp, &hR_serial, 0);
        const hamilt::HContainer<double>& hk_source = hR_serial;
        const bool do_write = (GlobalV::MY_RANK == 0);
#else
        const hamilt::HContainer<double>& hk_source = hR_tmp;
        const bool do_write = true;
#endif

        // Enumerate the local k points owned by this spin channel, exactly as
        // write_hsk does, so that filename_output derives the correct spin and
        // k index from ik2iktot.
        for (int ik = 0; ik < nks; ++ik)
        {
            const int kspin = kv.isk[ik];
            if (kspin != ispin)
            {
                continue;
            }
            write_hk_at_k(hk_source, spec.k_prefix, kv, ucell, nspin, nspin_name, ispin, istep,
                          append, ik, nlocal, params.precision, gamma_only, global_out_dir,
                          out_app_flag, spec.term_name, do_write, params.ofs_running);
        }

        if (also_hR)
        {
            gather_and_write(spec.r_prefix, spec.r_label, spec.term_name, hR_tmp, ucell, pv, nspin, ispin,
                             istep, out_app_flag, iat2iwt, nat,
                             params.global_matrix_dir, params.ofs_running);
        }
    }

    ModuleBase::timer::end("ModuleIO", spec.timer_label);
}

void write_h_vnl(WriteHParams& params)
{
    const UnitCell& ucell = *params.ucell;
    const Grid_Driver& gd = *params.gd;
    const TwoCenterBundle& two_center_bundle = *params.two_center_bundle;
    const LCAO_Orbitals& orb = *params.orb;
    const K_Vectors& kv = *params.kv;
    const std::vector<double>& orb_cutoff = orb.cutoffs();

    HTermSpec spec;
    spec.timer_label = "write_h_vnl";
    spec.k_prefix = "vnlk";
    spec.r_prefix = "vnl";
    spec.r_label = "V^NL";
    spec.term_name = "nonlocal";
    spec.spin_dependent = false;
    spec.build = [&](hamilt::HContainer<double>& hR_tmp, int /*ispin*/) {
        hamilt::Nonlocal<hamilt::OperatorLCAO<double, double>>
            tmp_nonlocal(nullptr, kv.kvec_d, &hR_tmp, &ucell, orb_cutoff, &gd,
                         two_center_bundle.overlap_orb_beta.get());
        tmp_nonlocal.contributeHR();
    };
    write_h_term(params, spec);
}

void write_h_vl(WriteHParams& params)
{
    const UnitCell& ucell = *params.ucell;
    const Grid_Driver& gd = *params.gd;
    const Parallel_Orbitals& pv = *params.pv;
    const LCAO_Orbitals& orb = *params.orb;
    const elecstate::Potential* pot = params.pot;
    const std::vector<double>& orb_cutoff = orb.cutoffs();

    const double* v_local = pot->get_fixed_v(); // local pp, no Hxc

    HTermSpec spec;
    spec.timer_label = "write_h_vl";
    spec.k_prefix = "vlk";
    spec.r_prefix = "vl";
    spec.r_label = "V^L";
    spec.term_name = "local";
    spec.spin_dependent = false;
    spec.build = [&](hamilt::HContainer<double>& hR_tmp, int /*ispin*/) {
        setup_veff_hcontainer(hR_tmp, ucell, gd, pv, orb_cutoff);
        // All ranks must enter cal_gint_vl so the Alltoall inside
        // hr_gint_to_hR stays matched. A grid-parallel rank with no local
        // real-space grid points (nrxx == 0) holds v_local == nullptr;
        // cal_gint_vl then has no local biggrid and never dereferences it.
        ModuleGint::cal_gint_vl(v_local, &hR_tmp);
    };
    write_h_term(params, spec);
}

void write_h_vh(WriteHParams& params)
{
    const UnitCell& ucell = *params.ucell;
    const Grid_Driver& gd = *params.gd;
    const Parallel_Orbitals& pv = *params.pv;
    const LCAO_Orbitals& orb = *params.orb;
    const Charge* chg = params.chg;
    const ModulePW::PW_Basis* rho_basis = params.rho_basis;
    const std::vector<double>& orb_cutoff = orb.cutoffs();

    ModuleBase::matrix v_h = elecstate::H_Hartree_pw::v_hartree(
        ucell, const_cast<ModulePW::PW_Basis*>(rho_basis), params.nspin, chg->rho);

    HTermSpec spec;
    spec.timer_label = "write_h_vh";
    spec.k_prefix = "vhk";
    spec.r_prefix = "vh";
    spec.r_label = "V^H";
    spec.term_name = "hartree";
    // The Hartree potential depends only on the total density, so both spin
    // rows of v_h are identical.
    spec.spin_dependent = false;
    spec.build = [&](hamilt::HContainer<double>& hR_tmp, int ispin) {
        setup_veff_hcontainer(hR_tmp, ucell, gd, pv, orb_cutoff);
        // All ranks must enter cal_gint_vl so the Alltoall inside
        // hr_gint_to_hR stays matched. Under grid-parallel distribution a
        // rank may own no real-space grid points (nrxx == 0), leaving v_h an
        // empty matrix; pass nullptr there (cal_gint_vl then has no local
        // biggrid and never dereferences it) instead of &v_h(ispin, 0), whose
        // matrix::operator() would trip the ic < nc assertion.
        const double* vh_ptr = (v_h.c != nullptr) ? &v_h(ispin, 0) : nullptr;
        ModuleGint::cal_gint_vl(vh_ptr, &hR_tmp);
    };
    write_h_term(params, spec);
}

void write_h_vxc(WriteHParams& params)
{
    const UnitCell& ucell = *params.ucell;
    const Grid_Driver& gd = *params.gd;
    const Parallel_Orbitals& pv = *params.pv;
    const LCAO_Orbitals& orb = *params.orb;
    const Charge* chg = params.chg;
    const int nrxx = params.nrxx;
    const std::vector<double>& orb_cutoff = orb.cutoffs();

    ModuleBase::matrix v_xc;
    double etxc;
    double vtxc;
    const double hybrid_alpha = XC_Functional::get_hybrid_alpha();
#ifdef __EXX
    const double hse_omega = XC_Functional::get_hse_omega();
#else
    const double hse_omega = 0.0;
#endif
    std::tie(etxc, vtxc, v_xc) = XC_Functional::v_xc(nrxx, chg, &ucell, params.nspin,
                                                     params.domag, params.domag_z,
                                                     params.gga_grad,
                                                     hybrid_alpha, hse_omega);

    HTermSpec spec;
    spec.timer_label = "write_h_vxc";
    spec.k_prefix = "vxck";
    spec.r_prefix = "vxc";
    spec.r_label = "V^XC";
    spec.term_name = "xc";
    spec.spin_dependent = true;
    spec.build = [&](hamilt::HContainer<double>& hR_tmp, int ispin) {
        setup_veff_hcontainer(hR_tmp, ucell, gd, pv, orb_cutoff);
        // All ranks must enter cal_gint_vl so the Alltoall inside
        // hr_gint_to_hR stays matched. A grid-parallel rank with no local
        // real-space grid points (nrxx == 0) holds an empty v_xc matrix; pass
        // nullptr there (never dereferenced) instead of &v_xc(ispin, 0),
        // whose matrix::operator() would trip the ic < nc assertion.
        const double* vxc_ptr = (v_xc.c != nullptr) ? &v_xc(ispin, 0) : nullptr;
        ModuleGint::cal_gint_vl(vxc_ptr, &hR_tmp);
    };
    write_h_term(params, spec);
}

#ifdef __EXX
// Build V^EXX(R) for one interface (real or complex Hexx) into real HContainers and write them.
template <typename Tdata>
static void write_h_exx_impl(const UnitCell& ucell,
                             const Parallel_Orbitals& pv,
                             Exx_LRI_Interface<double, Tdata>* ex,
                             const K_Vectors& kv,
                             const int nspin,
                             const int istep,
                             const bool append,
                             const int* iat2iwt,
                             const int nat,
                             const bool also_hR,
                             const int npol,
                             const int nlocal,
                             const int precision,
                             const bool gamma_only,
                             const std::string& global_out_dir,
                             const std::string& global_matrix_dir,
                             const std::string& calculation,
                             const int out_app_flag,
                             const Exx_Info& exx_info,
                             std::ofstream* ofs_running)
{
    const auto& Hexxs = ex->get_Hexxs(); // vector over spin of map<iat, map<(jat,R), Tensor>>
    const int nspin_out = (nspin == 2 ? 2 : 1);
    const int nks = kv.get_nks();
    const double alpha = exx_info.info_global.hybrid_alpha;

    for (int ispin = 0; ispin < nspin_out; ispin++)
    {
        hamilt::HContainer<double> hR_tmp(const_cast<Parallel_Orbitals*>(&pv));
        // add_HexxR only fills existing matrices, so first allocate the atom-pair structure
        // from the exx-form data (native cells, consistent with the nullptr cell_nearest below).
        hamilt::reallocate_hcontainer(Hexxs, &hR_tmp);
        RI_2D_Comm::add_HexxR(ispin, alpha, Hexxs, pv, npol, hR_tmp, nullptr);

#ifdef __MPI
        // Same as write_h_term: gather so that the folded H(k) is complete on
        // the writing rank.
        Parallel_Orbitals serialV;
        serialV.para_2d().set_serial(nlocal, nlocal);
        serialV.set_atomic_trace(iat2iwt, nat, nlocal);
        hamilt::HContainer<double> hR_serial(&serialV);
        hamilt::gatherParallels(hR_tmp, &hR_serial, 0);
        const hamilt::HContainer<double>& hk_source = hR_serial;
        const bool do_write = (GlobalV::MY_RANK == 0);
#else
        const hamilt::HContainer<double>& hk_source = hR_tmp;
        const bool do_write = true;
#endif

        for (int ik = 0; ik < nks; ++ik)
        {
            const int kspin = kv.isk[ik];
            if (kspin != ispin)
            {
                continue;
            }
            write_hk_at_k(hk_source, "vexxk", kv, ucell, nspin, nspin, ispin, istep, append, ik,
                          nlocal, precision, gamma_only, global_out_dir, out_app_flag, "exx",
                          do_write, ofs_running);
        }

        if (also_hR)
        {
            gather_and_write("vexx", "V^EXX", "exx", hR_tmp, ucell, pv, nspin, ispin, istep, out_app_flag,
                             iat2iwt, nat, global_matrix_dir,
                             ofs_running);
        }
    }
}

void write_h_exx(WriteHParams& params, const Exx_Info& exx_info)
{
    ModuleBase::TITLE("ModuleIO", "write_h_exx");
    ModuleBase::timer::start("ModuleIO", "write_h_exx");

    // Multi-k out_mat_exxk is rejected upstream at the call site (setup_exx_h_params in
    // ctrl_scf_lcao.cpp); this function is only reached on the gamma-only path.
    const UnitCell& ucell = *params.ucell;
    const Parallel_Orbitals& pv = *params.pv;
    const K_Vectors& kv = *params.kv;
    const int nspin = params.nspin;
    const int istep = params.istep;
    const bool append = params.append;
    const int* iat2iwt = params.iat2iwt;
    const int nat = params.nat;
    const bool also_hR = params.also_hR;

    // exd (real Hexx) and exc (complex Hexx) are mutually exclusive; pick by real_number.
    if (exx_info.info_ri.real_number)
    {
        if (params.exd != nullptr)
        {
            write_h_exx_impl(ucell, pv, params.exd, kv, nspin, istep, append, iat2iwt, nat, also_hR,
                             params.npol, params.nlocal, params.precision, params.gamma_only_local,
                             params.global_out_dir, params.global_matrix_dir,
                             params.calculation, params.out_app_flag, exx_info,
                             params.ofs_running);
        }
    }
    else
    {
        if (params.exc != nullptr)
        {
            write_h_exx_impl(ucell, pv, params.exc, kv, nspin, istep, append, iat2iwt, nat, also_hR,
                             params.npol, params.nlocal, params.precision, params.gamma_only_local,
                             params.global_out_dir, params.global_matrix_dir,
                             params.calculation, params.out_app_flag, exx_info,
                             params.ofs_running);
        }
    }

    ModuleBase::timer::end("ModuleIO", "write_h_exx");
}
#endif

} // namespace ModuleIO
