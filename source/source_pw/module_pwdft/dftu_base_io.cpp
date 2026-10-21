#include "source_pw/module_pwdft/dftu_base_io.h"

#include "source_cell/unitcell.h"
#include "source_pw/module_pwdft/dftu_base.h"
#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "source_base/parallel_common.h"
#include "source_base/parallel_global.h"
#include "source_base/timer.h"
#include "source_main/version.h"

#include <ctime>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>

namespace
{

// Write the STRU-style provenance header shared by occ_mat.txt and the
// per-ionic-step snapshot files: ABACUS version, local timestamp and the
// 1-based relaxation step this file belongs to.
void write_provenance_header(std::ostream& os, const int istep)
{
    std::time_t now = std::time(nullptr);
    std::tm tm_local;
    localtime_r(&now, &tm_local);
    char time_buf[64];
    std::strftime(time_buf, sizeof(time_buf), "%Y-%m-%d %H:%M:%S", &tm_local);
    os << "# ABACUS version: " << VERSION << "\n";
    os << "# Written at " << time_buf << "\n";
    os << "# RELAX STEP " << istep + 1 << "\n";
}

} // namespace


namespace DFTU_BASE
{

bool is_ion_step_output_step(const int istep, const OccmatOutputCfg& cfg)
{
    if (cfg.out_freq_ion <= 0)
    {
        return false;
    }
    return (istep % cfg.out_freq_ion == 0);
}

bool is_elec_snapshot_trigger(const int iter,
                              const bool conv_esolver,
                              const OccmatOutputCfg& cfg)
{
    const bool periodic = (iter % cfg.out_freq_elec == 0);
    const bool last_step = (iter == cfg.scf_nmax);
    return periodic || last_step || conv_esolver;
}

std::string gen_ion_step_occ_mat_filename(const std::string& out_dir, const int istep)
{
    std::stringstream ss;
    ss << out_dir << "occ_matg" << (istep + 1) << ".txt";
    return ss.str();
}

std::string find_first_existing_file(const std::string& dir,
                                     const std::vector<std::string>& candidates)
{
    for (const std::string& name : candidates)
    {
        const std::string full = dir + name;
        std::ifstream probe(full.c_str());
        if (probe.is_open())
        {
            return full;
        }
    }
    return std::string();
}

#ifdef __MPI
/// Broadcast the local occupation number matrices from rank 0 to all ranks.
///
/// Each occupation matrix is broadcast as one contiguous block
/// (matrix::c stores nr * nc consecutive doubles) instead of element by
/// element.
void local_occup_bcast(const UnitCell& ucell,
                       OccupationMatrix& occ,
                       const std::vector<int>& l_channel,
                       int nspin,
                       int npol)
{
    ModuleBase::TITLE("DFTU_BASE", "local_occup_bcast");

    for (int T = 0; T < ucell.ntype; T++)
    {
        if (l_channel[T] == -1)
        {
            continue;
        }

        for (int I = 0; I < ucell.atoms[T].na; I++)
        {
            const int iat = ucell.itia2iat(T, I);
            const int L = l_channel[T];

            for (int l = 0; l <= ucell.atoms[T].nwl; l++)
            {
                if (l != l_channel[T])
                {
                    continue;
                }

                if (nspin == 1 || nspin == 2)
                {
                    for (int spin = 0; spin < 2; spin++)
                    {
                        Parallel_Common::bcast_double(occ.mat(iat, l, spin).c,
                                                      occ.mat(iat, l, spin).nr * occ.mat(iat, l, spin).nc);
                    }
                }
                else if (nspin == 4) // SOC
                {
                    Parallel_Common::bcast_double(occ.mat(iat, l, 0).c,
                                                  occ.mat(iat, l, 0).nr * occ.mat(iat, l, 0).nc);
                }
            }
        }
    }
    return;
}
#endif


void prepare_ion_step_file(const std::string& global_out_dir,
                           const int istep,
                           const OccmatOutputCfg& cfg)
{
    // Only the root process creates the per-ionic-step file.
    if (GlobalV::MY_RANK != 0)
    {
        return;
    }

    // trunc makes a rerun start from a clean file instead of appending
    // snapshots left over from a previous calculation.
    const std::string ion_step_fn = gen_ion_step_occ_mat_filename(global_out_dir, istep);
    std::ofstream ofs_ion_step;
    ofs_ion_step.open(ion_step_fn, std::ios::out | std::ios::trunc);
    if (!ofs_ion_step)
    {
        ModuleBase::WARNING_QUIT("DFTU_BASE::prepare_ion_step_file",
                                 "Can't create per-ionic-step occupation-matrix file");
    }
    // Provenance header shared with occ_mat.txt.
    write_provenance_header(ofs_ion_step, istep);
    ofs_ion_step.close();
    return;
}


void output(const Plus_U_Base& dftu,
            const UnitCell& ucell,
            const std::string& global_out_dir,
            int nspin,
            int npol,
            int istep,
            int iter,
            const OccmatOutputCfg& cfg,
            OccmatSocLayout soc_layout)
{
    ModuleBase::TITLE("DFTU_BASE", "output");

    if (istep < 0)
    {
        ModuleBase::WARNING_QUIT("DFTU_BASE::output", "istep must be >= 0");
    }
    if (iter < 1)
    {
        ModuleBase::WARNING_QUIT("DFTU_BASE::output", "iter must be >= 1");
    }
    if (cfg.out_freq_ion < 0 || cfg.out_freq_elec < 1 || cfg.scf_nmax < 1)
    {
        ModuleBase::WARNING_QUIT("DFTU_BASE::output",
                                "invalid occupation-matrix output frequency configuration");
    }

    // The per-type U values and the per-atom occupation matrices are not
    // dumped into the running log here: for large cells (thousands of
    // atoms) this floods the log at every electronic step. The dedicated
    // occupation-matrix files (occ_mat.txt, occ_matg{#}.txt) gated by
    // out_occ_mat carry the same information instead.

    // occ_mat.txt is not written here: its section records drho, which
    // is only available after the electronic solve. write_latest_occmat()
    // overwrites it at the iter_finish stage.

    // At the first electronic step of an output ionic step, (re)create the
    // per-ionic-step file and write its header. The LCAO path enters here
    // directly; the PW path additionally prepares the file from
    // iter_init_dftu_pw(), because it skips output() at istep 0 / iter 1.
    // Electronic-step sections are appended later by append_ion_step_snapshot().
    const bool ion_step_output = is_ion_step_output_step(istep, cfg);
    if (cfg.out_occ_mat && ion_step_output && iter == 1)
    {
        prepare_ion_step_file(global_out_dir, istep, cfg);
    }

    return;
}

// Write one electronic-step section shared by occ_mat.txt and the
// per-ionic-step snapshot files: the step marker, the configured
// charge-density convergence threshold, the actual charge-density
// residual and the full occupation matrices. When occmat_ready is false,
// the matrix body is an "N/A" placeholder.
void write_snapshot_section(std::ostream& os,
                            const Plus_U_Base& dftu,
                            const UnitCell& ucell,
                            const int iter,
                            const double scf_thr,
                            const double drho,
                            const bool occmat_ready,
                            const int nspin,
                            const int npol,
                            const OccmatSocLayout soc_layout)
{
    os << "# Electronic step " << iter << "\n";
    os << "# scf_thr " << std::scientific << std::setprecision(8) << scf_thr << "\n";
    os << "# drho " << std::scientific << std::setprecision(8) << drho << "\n";

    // The PW path has no occupation matrix at istep 0 / iter 1 unless it was
    // loaded from file; record an explicit placeholder instead of a silent
    // zero matrix. Such an incomplete snapshot is diagnostic output and is
    // intentionally rejected by the strict restart reader.
    if (!occmat_ready)
    {
        os << "\n N/A\n";
        return;
    }

    // Full occupation matrices together with per-atom magnetism.
    write_occup_m(dftu, ucell, os, nspin, soc_layout);
}

void append_ion_step_snapshot(const Plus_U_Base& dftu,
                              const UnitCell& ucell,
                              const std::string& global_out_dir,
                              int nspin,
                              int npol,
                              int istep,
                              int iter,
                              bool conv_esolver,
                              bool occmat_ready,
                              double scf_thr,
                              double drho,
                              const OccmatOutputCfg& cfg,
                              OccmatSocLayout soc_layout)
{
    ModuleBase::TITLE("DFTU_BASE", "append_ion_step_snapshot");

    // dft_plus_u <= 0 means the occupation matrix was never computed
    // (Plus_U_Base::l_channel stays empty); writing here would dereference
    // a null l_channel.data() in write_occup_m. Match the documented contract
    // that out_occ_mat only takes effect for DFT+U calculations. Check this
    // before validating iter/istep/nspin: when DFT+U is off this function is
    // a no-op, so an invalid iter (e.g. from an upstream bug) must not crash
    // here. Mirrors the guard-first pattern in write_latest_occmat.
    if (!cfg.out_occ_mat || cfg.dft_plus_u <= 0)
    {
        return;
    }

    if (nspin != 1 && nspin != 2 && nspin != 4)
    {
        ModuleBase::WARNING_QUIT("DFTU_BASE::append_ion_step_snapshot", "nspin must be 1, 2 or 4");
    }
    if (istep < 0)
    {
        ModuleBase::WARNING_QUIT("DFTU_BASE::append_ion_step_snapshot", "istep must be >= 0");
    }
    if (iter < 1)
    {
        ModuleBase::WARNING_QUIT("DFTU_BASE::append_ion_step_snapshot", "iter must be >= 1");
    }
    // out_freq_ion == 0 is the valid default (no numbered file) and is
    // handled by the gates below; only negative values are an error.
    if (cfg.out_freq_ion < 0 || cfg.out_freq_elec < 1 || cfg.scf_nmax < 1)
    {
        ModuleBase::WARNING_QUIT("DFTU_BASE::append_ion_step_snapshot",
                                "invalid occupation-matrix output frequency configuration");
    }

    const bool ion_step_output = is_ion_step_output_step(istep, cfg);
    const bool elec_trigger = is_elec_snapshot_trigger(iter, conv_esolver, cfg);
    if (!ion_step_output || !elec_trigger || GlobalV::MY_RANK != 0)
    {
        return;
    }

    const std::string ion_step_fn = gen_ion_step_occ_mat_filename(global_out_dir, istep);

    // The provenance header written by prepare_ion_step_file() must not be
    // counted as a recorded section: scan the file for an existing section
    // marker instead of testing the stream position.
    bool section_exists = false;
    {
        std::ifstream checker(ion_step_fn.c_str());
        std::string line;
        while (std::getline(checker, line))
        {
            if (line.rfind("# Electronic step", 0) == 0)
            {
                section_exists = true;
                break;
            }
        }
    }

    std::ofstream ofs_ion_step;
    ofs_ion_step.open(ion_step_fn, std::ios::app);
    if (!ofs_ion_step)
    {
        ModuleBase::WARNING_QUIT("DFTU_BASE::append_ion_step_snapshot",
                                 "Can't open per-ionic-step occupation-matrix file");
    }

    // Separate consecutive sections by a blank line; the first section
    // directly follows the provenance header.
    if (section_exists)
    {
        ofs_ion_step << "\n";
    }
    write_snapshot_section(ofs_ion_step, dftu, ucell, iter, scf_thr, drho,
                           occmat_ready, nspin, npol, soc_layout);

    ofs_ion_step.close();

    GlobalV::ofs_running << " Append local occupation matrices of electronic step " << iter
                         << " to file: " << ion_step_fn << std::endl;
    return;
}


void write_latest_occmat(const Plus_U_Base& dftu,
                         const UnitCell& ucell,
                         const std::string& global_out_dir,
                         int nspin,
                         int npol,
                         int istep,
                         int iter,
                         double scf_thr,
                         double drho,
                         const OccmatOutputCfg& cfg,
                         OccmatSocLayout soc_layout)
{
    ModuleBase::TITLE("DFTU_BASE", "write_latest_occmat");

    if (!cfg.out_occ_mat || cfg.dft_plus_u <= 0 || GlobalV::MY_RANK != 0)
    {
        return;
    }

    const std::string latest_fn = global_out_dir + "occ_mat.txt";
    std::ofstream ofdftu;
    ofdftu.open(latest_fn, std::ios::out | std::ios::trunc);
    if (!ofdftu)
    {
        ModuleBase::WARNING_QUIT("DFTU_BASE::write_latest_occmat", "Can't create file occ_mat.txt");
    }

    // occ_mat.txt is a single-section snapshot file: the provenance
    // header shared with the g files followed by the same section writer.
    write_provenance_header(ofdftu, istep);
    // Callers of write_latest_occmat() guarantee the matrix exists: the PW
    // esolver gates on latest_ready, the LCAO esolver computes it every
    // electronic iteration.
    write_snapshot_section(ofdftu, dftu, ucell, iter, scf_thr, drho,
                           true, nspin, npol, soc_layout);

    ofdftu.close();
    return;
}


} // namespace DFTU_BASE
