#include "source_main/driver.h"

#include "source_base/global_file.h"
#include "source_base/memory.h"
#include "source_base/timer.h"
#include "source_esolver/esolver.h"
#include "source_pw/hamilt_pwdft/global.h"
#include "source_io/cal_test.h"
#include "source_io/input_conv.h"
#include "source_io/para_json.h"
#include "source_io/print_info.h"
#include "source_io/read_input.h"
#include "source_io/winput.h"
#include "module_parameter/parameter.h"
#include "source_main/version.h"
Driver::Driver()
{
}

Driver::~Driver()
{
}

void Driver::init()
{
    // 1) Let's start by printing a title.
    ModuleBase::TITLE("Driver", "ABACUS_begins");

    // 2) Print the current time, since it may run a long time.
    time_t time_start = std::time(nullptr);
    ModuleBase::timer::start();

    // 3) Welcome to the atomic world! Let's do some fancy stuff here.
    this->atomic_world();

    // 4) All timers recorders are printed.
    ModuleBase::timer::finish(GlobalV::ofs_running);

    // 5) All memory recorders are printed.
    ModuleBase::Memory::print_all(GlobalV::ofs_running);

    // 6) Print the final time, hopefully it will not cost too long. 
    time_t time_finish = std::time(nullptr);
    ModuleIO::print_time(time_start, time_finish);

    // 7) Clean up: close all of the running logs
    ModuleBase::Global_File::close_all_log(GlobalV::MY_RANK, PARAM.inp.out_alllog,PARAM.inp.calculation);

}

void Driver::print_start_info()
{
    ModuleBase::TITLE("Driver", "print_start_info");
#ifdef VERSION
    const char* version = VERSION;
#else
    const char* version = "unknown";
#endif
#ifdef COMMIT_INFO
#include "commit.h"
    const char* commit = COMMIT;
#else
    const char* commit = "unknown";
#endif
    time_t time_now = time(nullptr);

    PARAM.set_start_time(time_now);
    GlobalV::ofs_running << "                                                  "
                            "                                   "
                         << std::endl;
    GlobalV::ofs_running << "                              ABACUS " << version << std::endl << std::endl;
    GlobalV::ofs_running << "               Atomic-orbital Based Ab-initio "
                            "Computation at UStc                    "
                         << std::endl
                         << std::endl;
    GlobalV::ofs_running << "                     Website: http://abacus.ustc.edu.cn/           "
                            "                  "
                         << std::endl;
    GlobalV::ofs_running << "               Documentation: https://abacus.deepmodeling.com/     "
                            "                  "
                         << std::endl;
    GlobalV::ofs_running << "                  Repository: "
                            "https://github.com/abacusmodeling/abacus-develop       "
                         << std::endl;
    GlobalV::ofs_running << "                              "
                            "https://github.com/deepmodeling/abacus-develop         "
                         << std::endl;
    GlobalV::ofs_running << "                      Commit: " << commit << std::endl << std::endl;
    GlobalV::ofs_running << std::setiosflags(std::ios::right);

#ifndef __MPI
    GlobalV::ofs_running << "    This is SERIES version." << std::endl;
#endif
    GlobalV::ofs_running << "    Start Time is " << ctime(&time_now);
    GlobalV::ofs_running << "                                                  "
                            "                                   "
                         << std::endl;
    GlobalV::ofs_running << " -------------------------------------------------"
                            "-----------------------------------"
                         << std::endl;

    GlobalV::ofs_running << std::setiosflags(std::ios::left);
    std::cout << std::setiosflags(std::ios::left);

    GlobalV::ofs_running << "\n READING GENERAL INFORMATION" << std::endl;
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "global_out_dir", PARAM.globalv.global_out_dir);
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "global_in_card",  PARAM.globalv.global_in_card);
}

void Driver::reading()
{
    ModuleBase::TITLE("Driver", "reading");
    ModuleBase::timer::tick("Driver", "reading");
    // temperarily
    GlobalV::MY_RANK = PARAM.globalv.myrank;
    GlobalV::NPROC = PARAM.globalv.nproc;

    // (1) read the input file
    ModuleIO::ReadInput read_input(PARAM.globalv.myrank);
    read_input.read_parameters(PARAM,  PARAM.globalv.global_in_card);

    // (2) create the output directory, running_*.log and print info
    read_input.create_directory(PARAM);
    this->print_start_info();

    // (3) write the input file
    std::stringstream ss1;
    ss1 << PARAM.globalv.global_out_dir <<  PARAM.globalv.global_in_card;
    read_input.write_parameters(PARAM, ss1.str());

    // (*temp*) copy the variables from INPUT to each class
    Input_Conv::Convert();

    // (4) define the 'DIAGONALIZATION' world in MPI
    Parallel_Global::split_diag_world(PARAM.inp.diago_proc,
                                      GlobalV::NPROC,
                                      GlobalV::MY_RANK,
                                      GlobalV::DRANK,
                                      GlobalV::DSIZE,
                                      GlobalV::DCOLOR);
    Parallel_Global::split_grid_world(PARAM.inp.diago_proc,
                                      GlobalV::NPROC,
                                      GlobalV::MY_RANK,
                                      GlobalV::GRANK,
                                      GlobalV::GSIZE);
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "DRANK", GlobalV::DRANK + 1);
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "DSIZE", GlobalV::DSIZE);
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "DCOLOR", GlobalV::DCOLOR + 1);
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "GRANK", GlobalV::GRANK + 1);
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "GSIZE", GlobalV::GSIZE);

#ifdef __MPI
    // (5)  divide the GlobalV::NPROC processors into GlobalV::KPAR for k-points
    // parallelization.
    Parallel_Global::init_pools(GlobalV::NPROC,
                                GlobalV::MY_RANK,
                                PARAM.inp.bndpar,
                                GlobalV::KPAR,
                                GlobalV::NPROC_IN_BNDGROUP,
                                GlobalV::RANK_IN_BPGROUP,
                                GlobalV::MY_BNDGROUP,
                                GlobalV::NPROC_IN_POOL,
                                GlobalV::RANK_IN_POOL,
                                GlobalV::MY_POOL);
#endif

    // (6) Read in parameters about wannier functions.
    winput::Init(PARAM.inp.wannier_card);

    ModuleBase::timer::tick("Driver", "reading");
}

void Driver::atomic_world()
{
    ModuleBase::TITLE("Driver", "atomic_world");
    ModuleBase::timer::tick("Driver", "atomic_world");

    // reading information 
    this->reading();

    // where the actual stuff is done
    this->driver_run();

    ModuleBase::timer::tick("Driver", "atomic_world");
}
