#include "hsk_writer.h"
#include "hs_dense_io.h"
#include "source_base/module_external/scalapack_connector.h"
#include "source_base/module_out/filename.h" // use filename_output function
#include "source_base/global_function.h"
#include "source_base/parallel_common.h"
#include "source_base/timer.h"
#include "source_cell/klist.h"
#include "source_cell/ucell_io.h"
#include "source_cell/unitcell.h"
#include "source_io/module_output/spin_tag.h"
#include "source_io/module_output/write_notice.h"

#include <complex>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>

namespace
{
// Threshold below which a complex matrix element is flushed to zero in the
// text output so that tiny numerical noise does not clutter the file.
constexpr double k_hk_zero_thr = 1.0e-16;
} // namespace

template <typename T>
void ModuleIO::write_hsk(
        const std::string &global_out_dir,
        const int nspin,
        const int nks,
        const int nkstot,
        const std::vector<int> &ik2iktot,
        const std::vector<int> &isk,
        hamilt::Hamilt<T>* p_hamilt,
        const Parallel_Orbitals &pv,
        const K_Vectors &kv,
        const UnitCell* ucell,
        const std::vector<double> &efs,
        const bool gamma_only,
        const int out_app_flag,
        const int istep,
        const int out_type,
        const int precision,
        const int nlocal,
        const std::string &ks_solver,
        const int drank,
        std::ofstream &ofs_running)
{
    ModuleBase::TITLE("ModuleIO", "write_hsk");
    ModuleBase::timer::start("ModuleIO", "write_hsk");

    // Binary output keeps the save_mat layout unchanged so callers that
    // only consume .dat files (e.g. post-processing scripts) are not
    // affected by the text-format alignment.
    const bool binary = (out_type == 2);

    // The text body is documented as row-major (i, j). A column-major solver
    // stores its local block as mat[i + j * nrloc], and the Cpxgemr2d serial
    // destination keeps that layout, so reading mat[i * nlocal + j] would emit
    // the transpose (H(k)^T = conj(H(k)) for a Hermitian matrix). Select the
    // traversal order from the solver storage convention instead.
    const bool column_major_source
        = ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER(ks_solver);
    auto dense_index = [nlocal, column_major_source](const int i, const int j) -> size_t {
        return column_major_source ? static_cast<size_t>(j) * nlocal + i
                                   : static_cast<size_t>(i) * nlocal + j;
    };

    ofs_running << "\n WRITE H(k) OR S(k)" << std::endl;

    // text mode: align with write_dmk filename convention so that hk and dmk
    // files share the same naming scheme (no .dat suffix).  Binary mode
    // still uses filename_output which appends the .dat/.bin suffix.
    //
    // The filename's k index is the IBZ k-point index (0..nkstot_nospin-1),
    // not the global k-spin pair index (0..nkstot-1).  That way hk1s1_nao.txt
    // and hk1s2_nao.txt refer to the same IBZ k-point, matching dmk naming.
    auto make_fname = [&](const std::string& property, const int ik_global,
                          const int ibz_ik, const int ispin,
                          const bool spin_suffix) -> std::string {
        if (binary)
        {
            return ModuleIO::filename_output(global_out_dir, property, "nao",
                    ik_global, ik2iktot, nspin, nkstot, out_type, out_app_flag,
                    gamma_only, istep);
        }

        // Text path mirrors dmk_gen_fname so that hk1s1_nao.txt and
        // dmk1s1_nao.txt stay side-by-side and easy to compare.
        std::string fn = global_out_dir + property;
        if (!gamma_only)
        {
            fn += "k" + std::to_string(ibz_ik + 1);
        }
        if (spin_suffix)
        {
            fn += "s" + std::to_string(ispin + 1);
        }
        if (out_app_flag == 2)
        {
            fn += "g" + std::to_string(istep + 1);
        }
        fn += "_nao.txt";
        return fn;
    };

    // Each multi-step run (MD or relax) appends a new record to the fixed
    // file, mirroring the open_sparse_file convention; single-point SCF
    // passes istep = -1 and starts a fresh file. The first written step
    // (istep == 0) truncates so a fresh run drops records left by a previous
    // run; only later steps (istep > 0) append.
    const bool append_on_restart = (out_app_flag == 1) && (istep > 0);

    const Parallel_Blkcyclic& pv_dist = pv.para_2d();
    const int my_rank = drank;
    Parallel_Blkcyclic pv_glb;

    for (int ik = 0; ik < nks; ++ik)
    {
        p_hamilt->updateHk(ik);

        hamilt::MatrixBlock<T> h_mat;
        hamilt::MatrixBlock<T> s_mat;
        p_hamilt->matrix(h_mat, s_mat);

        // ------------------------------------------------------------------
        // H(k) — always written for every k-point
        // ------------------------------------------------------------------
        const int ispin = isk[ik];
        // IBZ k-point index: after symmetry reduction nkstot has already
        // been shrunk to IBZ size and then doubled for spin-polarised
        // calculations, so IBZ_k_count = nkstot / nspin.  ispin == 0 owns
        // ik in [0, IBZ_k_count), ispin == 1 owns [IBZ_k_count, nkstot),
        // hence ibz_ik = ik - ispin * IBZ_k_count recovers the irreducible
        // label that matches dmk/s1, dmk/s2, etc.
        const int ibz_k_count = kv.get_nkstot() / nspin;
        const int ibz_ik = ik - ispin * ibz_k_count;
        // Text path property is the letter before "k" so that make_fname
        // produces hk1s1 and sk1 matching dmk1s1.  Binary path keeps the
        // full "hk"/"sk" property since filename_output handles its own
        // suffix scheme.
        const std::string h_property = binary ? "hk" : "h";
        const bool hk_spin_suffix = (nspin == 2);
        const std::string h_fn = make_fname(h_property, ik, ibz_ik, ispin,
                hk_spin_suffix);

        if (binary)
        {
            // save_mat takes a Boolean append flag; modes 0 and 2 overwrite.
            const bool append_mat = append_on_restart;
            ModuleIO::save_mat(istep,
                    h_mat.p,
                    nlocal,
                    binary,
                    precision,
                    true, // triangular
                    append_mat,
                    h_fn,
                    pv_dist,
                    drank,
                    ks_solver);
        }
        else // text mode
        {
            // Gather the distributed H(k) matrix to rank 0 so we can write
            // the full matrix body without MPI-aware row-by-row logic.
            std::vector<T> h_global(my_rank == 0 ? nlocal * nlocal : 0);
#ifdef __MPI
            pv_glb.set(nlocal, nlocal, nlocal, pv_dist.ictxt());
            Cpxgemr2d(
            nlocal, nlocal,
            h_mat.p, 1, 1,
            const_cast<int*>(pv_dist.desc()),
            h_global.data(), 1, 1,
            const_cast<int*>(pv_glb.desc()),
            pv_glb.ictxt());
#else
            h_global.assign(h_mat.p, h_mat.p + nlocal * nlocal);
#endif

            if (my_rank == 0)
            {
                std::ios_base::openmode mode = std::ios::out;
                if (append_on_restart)
                {
                    mode |= std::ios::app;
                }
                std::ofstream ofs(h_fn.c_str(), mode);
                if (!ofs.is_open())
                {
                    ModuleBase::WARNING("ModuleIO::write_hsk",
                            "Cannot open H(k) file < " + h_fn + " >.");
                }
                else
                {
                    ModuleIO::write_matrix_notice(ofs_running, "H(k)", ispin, nspin, h_fn);

                    // --- metadata block (mirrors write_dmk exactly) ---
                    ofs << "Ionic_Step " << istep + 1 << std::endl;
                    ofs << " " << nspin << " # number of spin directions" << std::endl;
                    ofs << " " << ispin + 1 << " # spin index" << std::endl;
                    ofs << " " << kv.get_nkstot_nospin() << " # total k points " << std::endl;
                    ofs << " " << kv.get_nkstot() << " # total k points after symmetrized (if open) " << std::endl;
                    ofs << " " << ibz_ik + 1 << " # k-point index " << std::endl;
                    ofs << " " << kv.kvec_c[ibz_ik].x << " " << kv.kvec_c[ibz_ik].y << " " << kv.kvec_c[ibz_ik].z
                        << " # k point coordinate (Cartesian) " << std::endl;
                    ofs << " " << kv.kvec_d[ibz_ik].x << " " << kv.kvec_d[ibz_ik].y << " " << kv.kvec_d[ibz_ik].z
                        << " # k point coordinate (direct) " << std::endl;
                    ofs << " " << kv.wk[ibz_ik] << " # weight of this k point" << std::endl;
                    ofs << " " << efs[ispin] << " # Fermi energy in Ry " << std::endl;
                    ofs << " " << nlocal << " # number of localized basis " << std::endl;
                    ofs << " " << nlocal << " " << nlocal << " # size of this matrix " << std::endl;
                    ofs << std::endl;

                    // --- UnitCell lattice block ---
                    if (ucell != nullptr)
                    {
                        ModuleIO::UcellIO::write_ucell(ofs, ucell);
                    }
                    ofs << std::endl;

                    // --- matrix body: flat, 4 complex / 8 double per line ---
                    ofs << std::fixed;
                    ofs << std::scientific;
                    ofs << std::setprecision(precision);
                    ofs << std::right;
                    for (int i = 0; i < nlocal; ++i)
                    {
                        for (int j = 0; j < nlocal; ++j)
                        {
                            if (std::is_same<T, double>::value)
                            {
                                if (j % 8 == 0)
                                {
                                    ofs << "\n";
                                }
                                ofs << " " << h_global[dense_index(i, j)];
                            }
                            else if (std::is_same<T, std::complex<double>>::value)
                            {
                                if (j % 4 == 0)
                                {
                                    ofs << "\n";
                                }

                                double real_v = std::real(h_global[dense_index(i, j)]);
                                if (std::abs(real_v) < k_hk_zero_thr)
                                {
                                    real_v = 0.0;
                                }
                                double imag_v = std::imag(h_global[dense_index(i, j)]);
                                if (std::abs(imag_v) < k_hk_zero_thr)
                                {
                                    imag_v = 0.0;
                                }

                                ofs << " (" << real_v << "," << imag_v << ")";
                            }
                        }
                    }
                    // Terminate the matrix body with a newline so that the next
                    // appended record (out_app_flag == 1) starts on a fresh line.
                    ofs << std::endl;
                }
            } // rank0
        } // text mode

        // ------------------------------------------------------------------
        // S(k) — skip spin-down channels because both spins share the same S
        // ------------------------------------------------------------------
        const int current_spin = isk[ik];
        if (current_spin == 1)
        {
            continue;
        }
        const std::string s_property = binary ? "sk" : "s";
        // S(k) is spin-independent so we never append a spin suffix.
        const std::string s_fn = make_fname(s_property, ik, ibz_ik, 0, false);

        if (binary)
        {
            const bool append_mat = append_on_restart;
            ofs_running << " The output filename is " << s_fn << std::endl;
            ModuleIO::save_mat(istep,
                    s_mat.p,
                    nlocal,
                    binary,
                    precision,
                    true,
                    append_mat,
                    s_fn,
                    pv_dist,
                    drank,
                    ks_solver);
        }
        else // text mode
        {
            std::vector<T> s_global(my_rank == 0 ? nlocal * nlocal : 0);
#ifdef __MPI
            Cpxgemr2d(
            nlocal, nlocal,
            s_mat.p, 1, 1,
            const_cast<int*>(pv_dist.desc()),
            s_global.data(), 1, 1,
            const_cast<int*>(pv_glb.desc()),
            pv_glb.ictxt());
#else
            s_global.assign(s_mat.p, s_mat.p + nlocal * nlocal);
#endif

            if (my_rank == 0)
            {
                std::ios_base::openmode mode = std::ios::out;
                if (append_on_restart)
                {
                    mode |= std::ios::app;
                }
                std::ofstream ofs(s_fn.c_str(), mode);
                if (!ofs.is_open())
                {
                    ModuleBase::WARNING("ModuleIO::write_hsk",
                            "Cannot open S(k) file < " + s_fn + " >.");
                    continue;
                }
                ModuleIO::write_matrix_notice(ofs_running, "S(k)", s_fn);

                // --- metadata block ---
                ofs << "Ionic_Step " << istep + 1 << std::endl;
                ofs << " " << nspin << " # number of spin directions" << std::endl;
                ofs << " 1 # spin index" << std::endl;
                ofs << " " << kv.get_nkstot_nospin() << " # total k points " << std::endl;
                ofs << " " << kv.get_nkstot() << " # total k points after symmetrized (if open) " << std::endl;
                ofs << " " << ibz_ik + 1 << " # k-point index " << std::endl;
                ofs << " " << kv.kvec_c[ibz_ik].x << " " << kv.kvec_c[ibz_ik].y << " " << kv.kvec_c[ibz_ik].z
                    << " # k point coordinate (Cartesian) " << std::endl;
                ofs << " " << kv.kvec_d[ibz_ik].x << " " << kv.kvec_d[ibz_ik].y << " " << kv.kvec_d[ibz_ik].z
                    << " # k point coordinate (direct) " << std::endl;
                ofs << " " << kv.wk[ibz_ik] << " # weight of this k point" << std::endl;
                // S(k) has no Fermi-energy concept; write a zero placeholder
                // and label it clearly so readers are not confused.
                ofs << " 0.0 # Fermi energy in Ry (not applicable to S(k))" << std::endl;
                ofs << " " << nlocal << " # number of localized basis " << std::endl;
                ofs << " " << nlocal << " " << nlocal << " # size of this matrix " << std::endl;
                ofs << std::endl;

                // --- UnitCell lattice block ---
                if (ucell != nullptr)
                {
                    ModuleIO::UcellIO::write_ucell(ofs, ucell);
                }
                ofs << std::endl;

                // --- matrix body ---
                ofs << std::fixed;
                ofs << std::scientific;
                ofs << std::setprecision(precision);
                ofs << std::right;
                for (int i = 0; i < nlocal; ++i)
                {
                    for (int j = 0; j < nlocal; ++j)
                    {
                        if (std::is_same<T, double>::value)
                        {
                            if (j % 8 == 0)
                            {
                                ofs << "\n";
                            }
                            ofs << " " << s_global[dense_index(i, j)];
                        }
                        else if (std::is_same<T, std::complex<double>>::value)
                        {
                            if (j % 4 == 0)
                            {
                                ofs << "\n";
                            }

                            double real_v = std::real(s_global[dense_index(i, j)]);
                            if (std::abs(real_v) < k_hk_zero_thr)
                            {
                                real_v = 0.0;
                            }
                            double imag_v = std::imag(s_global[dense_index(i, j)]);
                            if (std::abs(imag_v) < k_hk_zero_thr)
                            {
                                imag_v = 0.0;
                            }

                            ofs << " (" << real_v << "," << imag_v << ")";
                        }
                    }
                }
                ofs << std::endl;
            } // rank0
        } // text mode
    } // end ik

    ModuleBase::timer::end("ModuleIO", "write_hsk");
}

// Explicit instantiations
template void ModuleIO::write_hsk<double>(
    const std::string&, const int, const int, const int,
    const std::vector<int>&, const std::vector<int>&,
    hamilt::Hamilt<double>*, const Parallel_Orbitals&,
    const K_Vectors&, const UnitCell*, const std::vector<double>&,
    const bool, const int, const int, const int, const int,
    const int, const std::string&, const int, std::ofstream&);
template void ModuleIO::write_hsk<std::complex<double>>(
    const std::string&, const int, const int, const int,
    const std::vector<int>&, const std::vector<int>&,
    hamilt::Hamilt<std::complex<double>>*, const Parallel_Orbitals&,
    const K_Vectors&, const UnitCell*, const std::vector<double>&,
    const bool, const int, const int, const int, const int,
    const int, const std::string&, const int, std::ofstream&);
