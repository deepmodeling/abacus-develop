#include "driver.h"

#include "module_base/memory.h"
#include "module_base/timer.h"
#include "module_esolver/esolver.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_io/cal_test.h"
#include "module_io/input.h"
#include "module_io/input_conv.h"
#include "module_io/print_info.h"
#include "module_io/winput.h"
#ifdef USE_PAW
#include "module_cell/module_paw/paw_cell.h"
#endif

Driver::Driver()
{
}

Driver::~Driver()
{
}

void Driver::init()
{
    ModuleBase::TITLE("Driver", "init");

    time_t time_start = std::time(NULL);
    ModuleBase::timer::start();

    // (1) read the input parameters.
    this->reading();

    // (2) welcome to the atomic world!
    this->atomic_world();

    // (3) output information
    time_t time_finish = std::time(NULL);
    Print_Info::print_time(time_start, time_finish);

    // (4) close all of the running logs
    INPUT.close_log();

    return;
}

void Driver::reading(void)
{
    ModuleBase::timer::tick("Driver", "reading");

    // (1) read INPUT
    INPUT.Init(GlobalV::global_in_card);

    // (2) copy the variables from INPUT to each class
    Input_Conv::Convert();

    // (3) define the 'DIAGONALIZATION' world in MPI
    Parallel_Global::split_diag_world(GlobalV::DIAGO_PROC);
    Parallel_Global::split_grid_world(GlobalV::DIAGO_PROC);
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "DRANK", GlobalV::DRANK + 1);
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "DSIZE", GlobalV::DSIZE);
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "DCOLOR", GlobalV::DCOLOR + 1);
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "GRANK", GlobalV::GRANK + 1);
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "GSIZE", GlobalV::GSIZE);

#ifdef __MPI
    // (4)  divide the GlobalV::NPROC processors into GlobalV::KPAR for k-points parallelization.
    Parallel_Global::init_pools();
#endif

    // (5) Read in parameters about wannier functions.
    winput::Init(GlobalV::global_wannier_card);

    // (6) Print the parameters into INPUT file.
    std::stringstream ss1;
    ss1 << GlobalV::global_out_dir << GlobalV::global_in_card;
    INPUT.Print(ss1.str());
    // ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running,"READING CARDS");

#ifdef USE_PAW
    if(GlobalV::use_paw)
    {
#ifdef __MPI
        GlobalC::paw_cell.io_redirect();
#endif
    }
#endif

    ModuleBase::timer::tick("Driver", "reading");
    return;
}

void Driver::atomic_world(void)
{
    ModuleBase::TITLE("Driver", "atomic_world");
    //--------------------------------------------------
    // choose basis sets:
    // pw: plane wave basis set
    // lcao_in_pw: LCAO expaned by plane wave basis set
    // lcao: linear combination of atomic orbitals
    //--------------------------------------------------

    // where the actual stuff is done
    this->driver_run();

    ModuleBase::timer::finish(GlobalV::ofs_running);
    ModuleBase::Memory::print_all(GlobalV::ofs_running);

    return;
}
