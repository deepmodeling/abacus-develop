#include <vector>
#include <algorithm>

#include "charge.h"
#include "chg_init.h"
#include "chg_atomic.h"
#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "source_base/libm/libm.h"
#include "source_base/timer.h"
#include "source_cell/magnetism.h"
#include "source_base/parallel_grid.h"
#include "source_io/module_output/cube_io.h"
#include "chg_rhog_io.h"
#include "source_io/module_wf/read_wf2rho_pw.h"
#include "source_io/module_restart/restart.h"
#include "source_cell/klist.h"
#include "source_base/module_parallel/para_world.h"
#include "source_base/module_parallel/para_tag.h"
#include "source_base/module_parallel/para_bridge.h"

// ---------------------------------------------------------------------------
// Thin member wrapper: Charge::init_rho delegates to the free function in
// module_charge so that the charge-init workflow lives outside the class.
// ---------------------------------------------------------------------------
void Charge::init_rho(const UnitCell& ucell,
                      const Parallel_Grid& pgrid,
                      const ModuleBase::ComplexMatrix& strucFac,
                      ModuleSymmetry::Symmetry& symm,
                      const void* klist,
                      const void* wfcpw,
                      const module_charge::InitRhoCfg& cfg)
{
    module_charge::init_rho(*this, *this->rhopw, ucell, pgrid, strucFac, symm, klist, wfcpw, cfg);
}

namespace module_charge
{

namespace
{

/// Aggregated file-reading configuration for read_rho_file / read_kin_file
struct ReadCfg
{
    const std::string& suffix;
    const std::string& readin_dir;
    int rank;
    std::ostream& ofs_running;
    std::ostream& ofs_warning;
};

/**
 * @brief Read charge density from restart binary or cube files into chr.rho.
 *
 * Charge members accessed: chr.rhog, chr.rho, chr.nspin.
 *
 * @param chr [inout] Charge object supplying the rho/rhog buffers
 * @param rhopw [in] plane-wave basis for file decoding and Fourier transforms
 * @param cfg [in] file-reading configuration (suffix, dir, rank, logs)
 * @param read_error [out] whether rho reading failed
 */
void read_rho_file(Charge& chr,
                   const ModulePW::PW_Basis& rhopw,
                   const UnitCell& ucell,
                   const Parallel_Grid& pgrid,
                   const ReadCfg& cfg,
                   bool& read_error)
{
    const int nspin = chr.nspin;
    std::complex<double>** const rhog = chr.rhog;
    double** const rho = chr.rho;
    const std::string& suffix = cfg.suffix;
    const std::string& readin_dir = cfg.readin_dir;
    const int rank = cfg.rank;
    std::ostream& ofs_running = cfg.ofs_running;
    std::ostream& ofs_warning = cfg.ofs_warning;

    ofs_running << " Read electron density from file" << std::endl;

    // try to read charge from binary file first, which is the same as QE
    // liuyu 2023-12-05
    std::stringstream binary;
    binary << readin_dir << suffix + "-CHARGE-DENSITY.restart";
    // Temporary bridge: use factory until ParaCollection is wired into driver.
    Parallel::ParaWorld pw_world = Parallel::make_pw_world();
    if (module_charge::read_rhog(binary.str(), &rhopw, nspin, rhog, pw_world, &ofs_warning))
    {
        ofs_running << " Read electron density from file: " << binary.str() << std::endl;
        for (int is = 0; is < nspin; ++is)
        {
            rhopw.recip2real(rhog[is], rho[is]);
        }
    }
    else
    {
        for (int is = 0; is < nspin; ++is)
        {
            std::stringstream ssc;

            if (nspin == 1)
            {
                ssc << readin_dir << "chg.cube";
            }
            else
            {
                ssc << readin_dir << "chgs" << is + 1 << ".cube";
            }

            if (ModuleIO::read_vdata_palgrid(pgrid,
                                             rank,
                                             ofs_running,
                                             ssc.str(),
                                             rho[is],
                                             ucell.nat))
            {
                ofs_running << " Read electron density from file: " << ssc.str() << std::endl;
            }
            else if (is > 0)    // nspin=2 or 4
            {
                if (is == 1)    // failed at the second spin
                {
                    std::cout << " Incomplete electron density file." << std::endl;
                    read_error = true;
                    break;
                }
                else if (is == 2)   // read 2 files when nspin=4
                {
                    ofs_running << " Didn't read in the electron density but would rearrange it later. "
                                << std::endl;
                }
                else if (is == 3)   // read 2 files when nspin=4
                {
                    ofs_running << " rearrange electron density " << std::endl;
                    for (int ir = 0; ir < rhopw.nrxx; ir++)
                    {
                        rho[3][ir] = rho[0][ir] - rho[1][ir];
                        rho[0][ir] = rho[0][ir] + rho[1][ir];
                        rho[1][ir] = 0.0;
                        rho[2][ir] = 0.0;
                    }
                }
            }
            else
            {
                read_error = true;
                break;
            }
        }
    }
}

/**
 * @brief Read kinetic-energy density from restart binary or cube files.
 *
 * Charge members accessed: chr.kin_r, chr.nspin.
 *
 * @param chr [inout] Charge object supplying the kin_r buffer
 * @param rhopw [in] plane-wave basis for file decoding and Fourier transforms
 * @param suffix [in] restart file prefix
 * @param readin_dir [in] directory to read from
 * @param rank [in] this processor's rank for palgrid reads
 * @param ofs_running [inout] running log stream
 * @param ofs_warning [inout] warning log stream
 * @param read_kin_error [out] whether kinetic-density reading failed
 */
void read_kin_file(Charge& chr,
                   const ModulePW::PW_Basis& rhopw,
                   const UnitCell& ucell,
                   const Parallel_Grid& pgrid,
                   const ReadCfg& cfg,
                   bool& read_kin_error)
{
    const int nspin = chr.nspin;
    double** const kin_r = chr.kin_r;
    const std::string& suffix = cfg.suffix;
    const std::string& readin_dir = cfg.readin_dir;
    const int rank = cfg.rank;
    std::ostream& ofs_running = cfg.ofs_running;
    std::ostream& ofs_warning = cfg.ofs_warning;

    ofs_running << " try to read kinetic energy density from file" << std::endl;
    std::vector<std::complex<double>> kin_g_space(nspin * rhopw.npw, {0.0, 0.0});
    std::vector<std::complex<double>*> kin_g;
    for (int is = 0; is < nspin; is++)
    {
        kin_g.push_back(kin_g_space.data() + is * rhopw.npw);
    }

    Parallel::ParaWorld pw_world = Parallel::make_pw_world();
    std::stringstream binary;
    binary << readin_dir << suffix + "-TAU-DENSITY.restart";
    if (module_charge::read_rhog(binary.str(), &rhopw, nspin, kin_g.data(), pw_world, &ofs_warning))
    {
        ofs_running << " Read in the kinetic energy density: " << binary.str() << std::endl;
        for (int is = 0; is < nspin; ++is)
        {
            rhopw.recip2real(kin_g[is], kin_r[is]);
        }
    }
    else
    {
        for (int is = 0; is < nspin; is++)
        {
            std::stringstream ssc;
            ssc << readin_dir << "SPIN" << is + 1 << "_TAU.cube";
            // mohan update 2012-02-10, sunliang update 2023-03-09
            if (ModuleIO::read_vdata_palgrid(
                    pgrid,
                    rank,
                    ofs_running,
                    ssc.str(),
                    kin_r[is],
                    ucell.nat))
            {
                ofs_running << " Read in the kinetic energy density: " << ssc.str() << std::endl;
            }
            else
            {
                read_kin_error = true;
                std::cout << " WARNING: \"init_chg\" is enabled but ABACUS failed to read kinetic energy "
                             "density from file.\n"
                             " Please check if there is SPINX_TAU.cube (X=1,...) or "
                             "{suffix}-TAU-DENSITY.restart in the directory.\n"
                          << std::endl;
                break;
            }
        }
    }
}

/**
 * @brief Atomic-density fallback plus Thomas-Fermi kinetic-energy-density init.
 *
 * Charge members accessed: chr.rho, chr.kin_r, chr.nspin.
 *
 * @param chr [inout] Charge object supplying rho/kin_r buffers
 * @param rhopw [in] plane-wave basis for atomic superposition and grid size
 * @param omega [in] unit-cell volume
 * @param init_chg [in] INPUT.init_chg
 * @param meta_gga [in] whether the functional is meta-GGA (tau TF-init needed)
 * @param read_error [in] whether rho reading failed
 * @param read_kin_error [in] whether kinetic-density reading failed
 */
void init_rho_atomic_and_tau(Charge& chr,
                             const ModulePW::PW_Basis& rhopw,
                             const UnitCell& ucell,
                             const ModuleBase::ComplexMatrix& strucFac,
                             const double& omega,
                             const std::string& init_chg,
                             const bool meta_gga,
                             const bool read_error,
                             const bool read_kin_error,
                             const AtomicRhoCfg& atomic_rho_cfg)
{
    const int nspin = chr.nspin;

    if (init_chg == "atomic" || read_error)
    {
        if (read_error)
        {
            std::cout << " Charge::init_rho: use atomic initialization instead." << std::endl;
        }
        module_charge::atomic_rho(nspin, omega, chr.rho, strucFac, ucell, &rhopw, atomic_rho_cfg);
    }

    // initial tau = 3/5 rho^2/3, Thomas-Fermi
    if (meta_gga)
    {
        if (init_chg == "atomic" || read_kin_error)
        {
            if (read_kin_error)
            {
                std::cout << " Charge::init_rho: init kinetic energy density from rho." << std::endl;
            }
            const double fact = (3.0 / 5.0) * pow(3.0 * ModuleBase::PI * ModuleBase::PI, 2.0 / 3.0);
            for (int is = 0; is < nspin; ++is)
            {
                for (int ir = 0; ir < rhopw.nrxx; ++ir)
                {
                    chr.kin_r[is][ir] = fact * pow(std::abs(chr.rho[is][ir]) * nspin, 5.0 / 3.0) / nspin;
                }
            }
        }
    }
}

/**
 * @brief Load charge density from the restart disk cache if requested.
 *
 * Charge members accessed: chr.rho, chr.nspin.
 *
 * @param chr [inout] Charge object supplying rho buffer
 * @param rhopw [in] plane-wave basis supplying the local real-space grid size
 * @param restart [inout] restart manager
 * @param readin_dir [in] fallback cube-file directory
 * @param rank [in] this processor's rank for palgrid reads
 * @param ofs_running [inout] running log stream
 */
void load_rho_from_restart(Charge& chr,
                           const ModulePW::PW_Basis& rhopw,
                           const UnitCell& ucell,
                           const Parallel_Grid& pgrid,
                           Restart& restart,
                           const std::string& readin_dir,
                           const int rank,
                           std::ostream& ofs_running)
{
    const int nspin = chr.nspin;

    // Peize Lin add 2020.04.04
    if (restart.info_load.load_charge && !restart.info_load.load_charge_finish)
    {
        for (int is = 0; is < nspin; ++is)
        {
            try
            {
                restart.load_disk("charge", is, rhopw.nrxx, chr.rho[is]);
            }
            catch (const std::exception& e)
            {
                // try to load from the output of `out_chg`
                std::stringstream ssc;
                ssc << readin_dir << "chgs" << is + 1 << ".cube";
                if (ModuleIO::read_vdata_palgrid(pgrid,
                                                 rank,
                                                 ofs_running,
                                                 ssc.str(),
                                                 chr.rho[is],
                                                 ucell.nat))
                {
                    ofs_running << " Read in electron density: " << ssc.str() << std::endl;
                }
            }
        }
        restart.info_load.load_charge_finish = true;
    }
}

} // anonymous namespace

// ---------------------------------------------------------------------------
// Public orchestrator: decides which initialization path(s) to run based on
// INPUT.init_chg and dispatches to the stage helpers above.
// ---------------------------------------------------------------------------
void init_rho(Charge& chr,
              const ModulePW::PW_Basis& rhopw,
              const UnitCell& ucell,
              const Parallel_Grid& pgrid,
              const ModuleBase::ComplexMatrix& strucFac,
              ModuleSymmetry::Symmetry& symm,
              const void* klist,
              const void* wfcpw,
              const InitRhoCfg& cfg)
{
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "init_chg", cfg.init_chg);

    const int nspin = chr.nspin;
    assert(nspin > 0);

    std::string init_chg_upper = cfg.init_chg;
    std::transform(init_chg_upper.begin(), init_chg_upper.end(), init_chg_upper.begin(), ::toupper);
    std::cout << " START CHARGE         : " << init_chg_upper << std::endl;

    chr.pgrid = &pgrid;

    const std::string& init_chg = cfg.init_chg;
    const std::string& suffix = cfg.suffix;
    const std::string& readin_dir = cfg.global_readin_dir;
    const int rank = (cfg.esolver_type == "sdft" ? GlobalV::RANK_IN_BPGROUP : GlobalV::MY_RANK);

    bool read_error = false;
    bool read_kin_error = false;
    // Capture before the local ReadCfg (also named cfg) shadows the argument.
    const bool meta_gga = cfg.meta_gga;
    if (init_chg == "file" || init_chg == "auto")
    {
        ReadCfg cfg{suffix, readin_dir, rank,
                    GlobalV::ofs_running, GlobalV::ofs_warning};
        read_rho_file(chr, rhopw, ucell, pgrid, cfg, read_error);

        if (read_error)
        {
            const std::string warn_msg
                = " WARNING: \"init_chg\" is enabled but ABACUS failed to read\n charge density from file.\n"
                  " Please check if there is chg.cube (for nspin=1) or"
                  " chgsx.cube (x=1,2,etc.) or\n"
                  " {suffix}-CHARGE-DENSITY.restart in the "
                  "directory.\n";
            std::cout << warn_msg;
            if (init_chg == "file")
            {
                ModuleBase::WARNING_QUIT("Charge::init_rho",
                                         "Failed to read in charge density from file.\n For initializing atomic "
                                         "charge in calculations,\n please set init_chg to atomic in INPUT.");
            }
        }

        // If the charge density is not read in, then the kinetic energy density is not read in either
        if (meta_gga)
        {
            if (!read_error)
            {
                read_kin_file(chr, rhopw, ucell, pgrid, cfg, read_kin_error);
            }
            else
            {
                read_kin_error = true;
            }
        }
    }

    const AtomicRhoCfg atomic_rho_cfg{
        cfg.nelec,
        cfg.test_charge,
        cfg.domag,
        cfg.domag_z,
        GlobalV::ofs_warning};
    init_rho_atomic_and_tau(chr, rhopw, ucell, strucFac, ucell.omega,
                            init_chg, cfg.meta_gga, read_error, read_kin_error,
                            atomic_rho_cfg);

    load_rho_from_restart(chr, rhopw, ucell, pgrid, GlobalC::restart,
                          readin_dir, rank, GlobalV::ofs_running);

    if (init_chg == "wfc")
    {
        if (wfcpw == nullptr)
        {
            ModuleBase::WARNING_QUIT("Charge::init_rho", "wfc is only supported for PW-KSDFT.");
        }

        const ModulePW::PW_Basis_K* pw_wfc = reinterpret_cast<ModulePW::PW_Basis_K*>(const_cast<void*>(wfcpw));
        const K_Vectors* kv = reinterpret_cast<const K_Vectors*>(klist);

        ModuleIO::read_wf2rho_pw(pw_wfc, symm, chr,
                                 readin_dir,
                                 GlobalV::KPAR, GlobalV::MY_POOL, GlobalV::MY_RANK,
                                 GlobalV::NPROC_IN_POOL, GlobalV::RANK_IN_POOL,
                                 cfg.nbands, nspin, cfg.npol,
                                 kv->get_nkstot(), kv->ik2iktot, kv->isk, GlobalV::ofs_running);
    }
}

} // namespace module_charge
