#include "driver.h"
#include "src_pw/global.h"
#include "input.h"
#include "src_io/winput.h"
#include "module_neighbor/sltk_atom_arrange.h"
#include "src_io/print_info.h"
#include "src_lcao/run_md_lcao.h"
#include "src_pw/run_md_pw.h"

// This is the driver function which defines the workflow of ABACUS calculations
// It relies on the class Esolver, which is a class that organizes workflows of single point calculations.
// For calculations involving change of configuration (lattice parameter & ionic motion),
// this driver calls Esolver::Run and the configuration-changing subroutine
// in a alternating manner.
// Information is passed between the two subroutines by class UnitCell_Pseudo
// Esolver::Run takes in a configuration and provides force and stress, 
// the configuration-changing subroutine takes force and stress and updates the configuration
void Driver::driver_run()
{
    ModuleBase::TITLE("Driver", "driver_line");
    ModuleBase::timer::tick("Driver", "driver_line");

    // 1. Initialzie type of Esolver
    ModuleESolver::ESolver *p_esolver = nullptr;
    ModuleESolver::init_esolver(p_esolver);

    // 2. Setup cell and atom information
#ifdef __LCAO
    GlobalC::ucell.setup_cell(GlobalC::ORB, GlobalV::global_pseudo_dir, GlobalV::stru_file, GlobalV::ofs_running);
#else
    if(GlobalV::BASIS_TYPE == "lcao_in_pw" || GlobalV::BASIS_TYPE == "lcao")
    {
        ModuleBase::WARNING_QUIT("driver","to use LCAO basis, compile with __LCAO");
    }
    GlobalC::ucell.setup_cell(GlobalV::global_pseudo_dir, GlobalV::stru_file, GlobalV::ofs_running);
#endif

    // 3. For these two types of calculations
    // nothing else need to be initialized
    if(GlobalV::CALCULATION == "test_neighbour" || GlobalV::CALCULATION == "test_memory")
    {
        p_esolver->Run(0, GlobalC::ucell);
        ModuleBase::QUIT();
    }

    // 4. Initialize Esolver
    p_esolver->Init(INPUT, GlobalC::ucell);

    //------------------------------------------------------------
    // This part onward needs to be refactored.
    //---------------------------MD/Relax------------------
    if(GlobalV::CALCULATION == "md" && GlobalV::BASIS_TYPE=="lcao")
    {
        Run_MD_LCAO run_md_lcao;
        run_md_lcao.opt_ions(p_esolver);
    }
    else if(GlobalV::CALCULATION == "md" || GlobalV::CALCULATION == "sto-md")
    {
        Run_MD_PW run_md_pw;
        run_md_pw.md_ions_pw(p_esolver);
    }
    else // scf; cell relaxation; nscf; etc
    {
        Ions ions;
        ions.opt_ions(p_esolver);
    }
    //---------------------------MD/Relax------------------

    // 6. clean up esolver
    p_esolver->postprocess();
    ModuleESolver::clean_esolver(p_esolver);

    ModuleBase::timer::tick("Driver", "driver_line");
    return;
}
