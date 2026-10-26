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

namespace DFTU_BASE
{

// Write the STRU-style provenance header shared by occ_mat.txt and the
// per-ionic-step snapshot files: ABACUS version, local timestamp and the
// 1-based relaxation step this file belongs to.
static void write_provenance_header(std::ostream& os, const int istep)
{
    std::time_t now = std::time(nullptr);
    std::tm tm_local;
    localtime_r(&now, &tm_local);
    char time_buf[64];
    std::strftime(time_buf, sizeof(time_buf), "%Y-%m-%d %H:%M:%S", &tm_local);
    os << "# ABACUS version: " << VERSION << "\n";
    os << "# Written at " << time_buf << "\n";
    os << "# RELAX STEP " << istep + 1 << "\n\n";
}

static void validate_output_configuration(const char* context,
                                          const int istep,
                                          const int iter,
                                          const OccmatOutputCfg& cfg)
{
    if (istep < 0)
    {
        ModuleBase::WARNING_QUIT(context, "istep must be >= 0");
    }
    if (iter < 1)
    {
        ModuleBase::WARNING_QUIT(context, "iter must be >= 1");
    }
    if (cfg.out_freq_ion < 0 || cfg.out_freq_elec < 1 || cfg.scf_nmax < 1)
    {
        ModuleBase::WARNING_QUIT(context,
                                 "invalid occupation-matrix output frequency configuration");
    }
}

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
                       const std::vector<int>& l_channel)
{
    ModuleBase::TITLE("DFTU_BASE", "local_occup_bcast");
    const int nspin = occ.nspin();

    for (int T = 0; T < ucell.ntype; T++)
    {
        if (l_channel[T] == -1)
        {
            continue;
        }
        const int l = l_channel[T];
        for (int I = 0; I < ucell.atoms[T].na; I++)
        {
            const int iat = ucell.itia2iat(T, I);
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
#endif


static void prepare_ion_step_file(const std::string& global_out_dir,
                                  const int istep)
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
}

void prepare_ion_step_output(const std::string& global_out_dir,
                             const int istep,
                             const int iter,
                             const OccmatOutputCfg& cfg)
{
    ModuleBase::TITLE("DFTU_BASE", "prepare_ion_step_output");
    validate_output_configuration("DFTU_BASE::prepare_ion_step_output",
                                  istep,
                                  iter,
                                  cfg);
    const bool ion_step_output = is_ion_step_output_step(istep, cfg);
    if (cfg.out_occ_mat && ion_step_output && iter == 1)
    {
        prepare_ion_step_file(global_out_dir, istep);
    }
}

// Write one electronic-step section shared by occ_mat.txt and the
// per-ionic-step snapshot files: the step marker, the configured
// charge-density convergence threshold, the actual charge-density
// residual and the full occupation matrices. When occmat_ready is false,
// the matrix body is an "N/A" placeholder.
static void write_snapshot_section(std::ostream& os,
                                   const Plus_U_Base& dftu,
                                   const UnitCell& ucell,
                                   const OccmatStepInfo& step,
                                   const bool occmat_ready,
                                   const OccmatSocLayout soc_layout)
{
    os << "# Electronic step " << step.iter << "\n";
    os << "# scf_thr " << std::scientific << std::setprecision(8) << step.scf_thr << "\n";
    os << "# drho " << std::scientific << std::setprecision(8) << step.drho << "\n";

    // The PW path has no occupation matrix at istep 0 / iter 1 unless it was
    // loaded from file; record an explicit placeholder instead of a silent
    // zero matrix. Such an incomplete snapshot is diagnostic output and is
    // intentionally rejected by the strict restart reader.
    if (!occmat_ready)
    {
        os << "\n N/A\n\n";
        return;
    }

    // Full occupation matrices together with per-atom magnetism.
    write_occup_m(dftu, ucell, os, soc_layout);
    os << "\n";
}

void append_ion_step_snapshot(const Plus_U_Base& dftu,
                              const UnitCell& ucell,
                              const std::string& global_out_dir,
                              const OccmatStepInfo& step,
                              const bool conv_esolver,
                              const bool occmat_ready,
                              const OccmatOutputCfg& cfg,
                              const OccmatSocLayout soc_layout)
{
    ModuleBase::TITLE("DFTU_BASE", "append_ion_step_snapshot");

    // dft_plus_u <= 0 means the occupation matrix was never computed
    // (Plus_U_Base::l_channel stays empty); writing here would dereference
    // a null l_channel.data() in write_occup_m. Match the documented contract
    // that out_occ_mat only takes effect for DFT+U calculations. Check this
    // before validating the step: when DFT+U is off this function is
    // a no-op, so an invalid iter (e.g. from an upstream bug) must not crash
    // here. Mirrors the guard-first pattern in write_latest_occmat.
    if (!cfg.out_occ_mat || cfg.dft_plus_u <= 0)
    {
        return;
    }

    validate_output_configuration("DFTU_BASE::append_ion_step_snapshot",
                                  step.istep,
                                  step.iter,
                                  cfg);

    const bool ion_step_output = is_ion_step_output_step(step.istep, cfg);
    const bool elec_trigger = is_elec_snapshot_trigger(step.iter, conv_esolver, cfg);
    if (!ion_step_output || !elec_trigger || GlobalV::MY_RANK != 0)
    {
        return;
    }

    const std::string ion_step_fn = gen_ion_step_occ_mat_filename(global_out_dir, step.istep);

    std::ofstream ofs_ion_step;
    ofs_ion_step.open(ion_step_fn, std::ios::app);
    if (!ofs_ion_step)
    {
        ModuleBase::WARNING_QUIT("DFTU_BASE::append_ion_step_snapshot",
                                 "Can't open per-ionic-step occupation-matrix file");
    }

    // A direct caller may not have prepared the file at iter 1. Add the
    // provenance header when append mode opens an empty file.
    if (ofs_ion_step.tellp() == std::streampos(0))
    {
        write_provenance_header(ofs_ion_step, step.istep);
    }
    write_snapshot_section(ofs_ion_step, dftu, ucell, step, occmat_ready, soc_layout);

    GlobalV::ofs_running << " Append local occupation matrices of electronic step " << step.iter
                         << " to file: " << ion_step_fn << std::endl;
}


void write_latest_occmat(const Plus_U_Base& dftu,
                         const UnitCell& ucell,
                         const std::string& global_out_dir,
                         const OccmatStepInfo& step,
                         const OccmatOutputCfg& cfg,
                         const OccmatSocLayout soc_layout)
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
    write_provenance_header(ofdftu, step.istep);
    // Callers of write_latest_occmat() guarantee the matrix exists: the PW
    // esolver gates on latest_ready, the LCAO esolver computes it every
    // electronic iteration.
    write_snapshot_section(ofdftu, dftu, ucell, step, true, soc_layout);
}


} // namespace DFTU_BASE
