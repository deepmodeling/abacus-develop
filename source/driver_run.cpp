#include "driver.h"
#include "module_cell/check_atomic_stru.h"
#include "module_cell/module_neighbor/sltk_atom_arrange.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_io/input.h"
#include "module_io/para_json.h"
#include "module_io/print_info.h"
#include "module_io/winput.h"
#include "module_md/run_md.h"

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
void Driver::driver_run() {
    ModuleBase::TITLE("Driver", "driver_line");
    ModuleBase::timer::tick("Driver", "driver_line");

    //! 1: initialize the ESolver 
    ModuleESolver::ESolver *p_esolver = ModuleESolver::init_esolver();

    //! 2: setup cell and atom information

    // this warning should not be here, mohan 2024-05-22
#ifndef __LCAO
    if (GlobalV::BASIS_TYPE == "lcao_in_pw" || GlobalV::BASIS_TYPE == "lcao") {
        ModuleBase::WARNING_QUIT("driver",
                                 "to use LCAO basis, compile with __LCAO");
    }
#endif

    // the life of ucell should begin here, mohan 2024-05-12
    // delete ucell as a GlobalC in near future
    GlobalC::ucell.setup_cell(GlobalV::stru_file, GlobalV::ofs_running);
    Check_Atomic_Stru::check_atomic_stru(GlobalC::ucell,
                                         GlobalV::MIN_DIST_COEF);

    //! 3: initialize Esolver and fill json-structure
    p_esolver->before_all_runners(INPUT, GlobalC::ucell);

    // this Json part should be moved to before_all_runners, mohan 2024-05-12
#ifdef __RAPIDJSON
    Json::gen_stru_wrapper(&GlobalC::ucell);
#endif

    const std::string cal_type = GlobalV::CALCULATION;

    //! 4: different types of calculations
    if (cal_type == "md") {
        Run_MD::md_line(GlobalC::ucell, p_esolver, INPUT.mdp);
    } else if (cal_type == "scf" || cal_type == "relax"
               || cal_type == "cell-relax") {
        Relax_Driver rl_driver;
        rl_driver.relax_driver(p_esolver);
    } else {
        //! supported "other" functions:
        //! nscf(PW,LCAO),
        //! get_pchg(LCAO),
        //! test_memory(PW,LCAO),
        //! test_neighbour(LCAO),
        //! get_S(LCAO),
        //! gen_bessel(PW), et al.
        const int istep = 0;
        p_esolver->others(istep);
    }

    //! 5: clean up esolver
    p_esolver->after_all_runners();
    ModuleESolver::clean_esolver(p_esolver);

    ModuleBase::timer::tick("Driver", "driver_line");
    return;
}
