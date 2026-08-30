#include "source_base/constants.h"
#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "source_base/kernels/math_kernel_op.h"
#include "source_base/module_device/device.h"
#include "source_base/module_device/memory_op.h"
#include "source_base/parallel_cell.h"
#include "source_cell/check_atomic_stru.h"
#include "source_cell/distributed_mdcell_reader.h"
#include "source_cell/md_cell.h"
#include "source_cell/module_neighbor/sltk_atom_arrange.h"
#include "source_cell/print_cell.h"
#include "source_esolver/esolver_factory.h"
#include "source_hsolver/kernels/hegvd_op.h"
#include "source_io/module_json/para_json.h"
#include "source_io/module_output/print_info.h"
#include "source_io/module_parameter/parameter.h"
#include "source_main/driver.h"
#include "source_md/run_md.h"
#include "source_relax/relax_driver.h"

#include <ATen/kernels/blas.h>
#include <ATen/kernels/lapack.h>

#ifdef __DSP
#include "source_base/kernels/dsp/dsp_connector.h"
#endif

/**
 * @brief This is the driver function which defines the workflow of ABACUS
 * calculations. It relies on the class Esolver, which is a class that organizes
 * workflows of single point calculations.
 *
 * For calculations involving change of configuration (lattice parameter & ionic
 * motion), this driver calls Esolver::Run and the configuration-changing
 * subroutine in a alternating manner.
 *
 * Information is passed between the two subroutines by class UnitCell
 *
 * Esolver::Run takes in a configuration and provides force and stress,
 * the configuration-changing subroutine takes force and stress and updates the
 * configuration
 */
void Driver::driver_run()
{
    ModuleBase::TITLE("Driver", "driver_run");

    //! 1: setup cell and atom information
    // this warning should not be here, mohan 2024-05-22
#ifndef __LCAO
    if (PARAM.inp.basis_type == "lcao_in_pw" || PARAM.inp.basis_type == "lcao")
    {
        ModuleBase::WARNING_QUIT("driver", "to use LCAO basis, compile with __LCAO");
    }
#endif

    const std::string cal = PARAM.inp.calculation;
    const Input_para& input = PARAM.inp;

    this->init_hardware();
    ModuleESolver::ESolver* p_esolver = ModuleESolver::init_esolver(PARAM.inp);

    // UnitCell is initialized only for workflows that require its full DFT state.
    UnitCell ucell;
    bool ucell_initialized = false;
    const auto initialize_ucell = [&ucell, &ucell_initialized, &input]()
    {
        if (ucell_initialized)
        {
            return;
        }

        ucell.setup_from_input(input.latname,
                               input.ntype,
                               input.lmaxmax,
                               input.init_vel,
                               input.fixed_axes);
        ucell.setup_cell(PARAM.globalv.global_in_stru,
                         GlobalV::ofs_running,
                         input.symmetry_prec,
                         input.dfthalf_type,
                         input.pseudo_dir,
                         input.nspin,
                         input.basis_type,
                         input.orbital_dir,
                         input.init_wfc,
                         input.onsite_radius,
                         PARAM.globalv.deepks_setorb,
                         input.rpa,
                         input.fixed_atoms,
                         input.noncolin,
                         input.calculation,
                         input.esolver_type,
                         std::stoi(input.symmetry));
        unitcell::check_atomic_stru(ucell, input.min_dist_coef);
        ucell_initialized = true;

#ifdef __RAPIDJSON
        Json::gen_stru_wrapper(&ucell, input);
#endif
    };

    if (cal == "md")
    {
        const ModuleBase::CommunicationDomain communication_domain = ModuleBase::world_communication_domain();
        if (p_esolver->supports_mdcell())
        {
            const Input_para& input = PARAM.inp;
            const double cutoff = p_esolver->mdcell_cutoff(input);
            if (cutoff <= 0.0)
            {
                ModuleBase::WARNING_QUIT("Driver::driver_run",
                                         "An ESolver supporting MDCell must provide a positive cutoff.");
            }
            const std::vector<int> effective_replicate = input.mdp.md_restart
                                                              ? std::vector<int>{1, 1, 1}
                                                              : input.cell_replica;
            MdStruFileMetadata stru_metadata;
            MDCell mdcell = DistributedMDCellReader::read_stru(PARAM.globalv.global_in_stru,
                                                                 effective_replicate,
                                                                 cutoff,
                                                                 input.mdp.md_neighbor_skin / ModuleBase::BOHR_TO_A,
                                                                 stru_metadata,
                                                                 communication_domain);
            GlobalV::ofs_running << std::endl;
            ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "TOTAL ATOM NUMBER", mdcell.nat());
            GlobalV::ofs_running << std::endl;
            p_esolver->before_all_runners(mdcell, input);
            Run_MD::md_line(mdcell, p_esolver, PARAM, stru_metadata);
            p_esolver->after_all_runners(mdcell);
        }
        else
        {
            initialize_ucell();
            MDCell mdcell(ucell, 0.0, 0.0, communication_domain);
            const MdStruFileMetadata stru_metadata = unitcell::make_md_stru_file_metadata(ucell);
            p_esolver->before_all_runners(ucell, PARAM.inp);
            Run_MD::md_line(mdcell, p_esolver, PARAM, stru_metadata);
            p_esolver->after_all_runners(ucell);
        }
    }
    else
    {
        initialize_ucell();
        p_esolver->before_all_runners(ucell, PARAM.inp);
        if (cal == "scf" || cal == "relax" || cal == "cell-relax" || cal == "nscf")
        {
            Relax_Driver rl_driver;
            rl_driver.relax_driver(p_esolver, ucell, PARAM.inp, GlobalV::ofs_running);
        }
        else if (cal == "get_s")
        {
            p_esolver->runner(ucell, 0);
        }
        else if (cal == "get_pchg" || cal == "get_wf" || cal == "gen_bessel" || cal == "gen_opt_abfs"
                 || cal == "test_memory" || cal == "test_neighbour")
        {
            p_esolver->others(ucell, 0);
        }
        else
        {
            ModuleBase::WARNING_QUIT("Driver::driver_run", "cannot recognize the 'calculation' command");
        }
        p_esolver->after_all_runners(ucell);
    }

    delete p_esolver;

    this->finalize_hardware();

    //! 6: output the json file
    if (ucell_initialized)
    {
        Json::create_Json(&ucell, PARAM);
    }

    return;
}

void Driver::init_hardware()
{
#if ((defined __CUDA) || (defined __ROCM))
    if (PARAM.inp.device == "gpu")
    {
        ModuleBase::createGpuBlasHandle();
        hsolver::createGpuSolverHandle();
        container::kernels::createGpuBlasHandle();
        container::kernels::createGpuSolverHandle();
    }
#endif

#ifdef __DSP
    std::cout << " ** Initializing DSP Hardware..." << std::endl;
    mtfunc::dspInitHandle(GlobalV::MY_RANK % PARAM.inp.dsp_count);
#endif
}

void Driver::finalize_hardware()
{
#if defined(__CUDA) || defined(__ROCM)
    if (PARAM.inp.device == "gpu")
    {
        ModuleBase::destoryBLAShandle();
        hsolver::destroyGpuSolverHandle();
        container::kernels::destroyGpuBlasHandle();
        container::kernels::destroyGpuSolverHandle();
    }
#endif

#ifdef __DSP
    std::cout << " ** Closing DSP Hardware..." << std::endl;
    mtfunc::dspDestoryHandle(GlobalV::MY_RANK % PARAM.inp.dsp_count);
#endif
}
