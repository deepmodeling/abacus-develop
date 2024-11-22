#include "driver.h"

#include "module_base/global_file.h"
#include "module_base/memory.h"
#include "module_base/timer.h"
#include "module_esolver/esolver.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_io/cal_test.h"
#include "module_io/input_conv.h"
#include "module_io/para_json.h"
#include "module_io/print_info.h"
#include "module_io/read_input.h"
#include "module_io/winput.h"
#include "module_parameter/parameter.h"
#include "version.h"
Driver::Driver()
{
}

Driver::~Driver()
{
    // Release the device memory within singleton object GlobalC::ppcell
    // before the main function exits.
    GlobalC::ppcell.release_memory();
}

void Driver::init()
{
    ModuleBase::TITLE("Driver", "init");

    time_t time_start = std::time(nullptr);
    ModuleBase::timer::start();

    // (1) read the input parameters.
    // INPUT should be initalized here and then pass to atomic world, mohan
    // 2024-05-12 INPUT should not be GlobalC, mohan 2024-05-12
    Driver::reading();

    // (2) welcome to the atomic world!
    this->atomic_world();

    // (3) output information
    time_t time_finish = std::time(nullptr);
    ModuleIO::print_time(time_start, time_finish);

    // (4) close all of the running logs
    ModuleBase::Global_File::close_all_log(GlobalV::MY_RANK, PARAM.inp.out_alllog,PARAM.inp.calculation);

    // (5) output the json file
    // Json::create_Json(&GlobalC::ucell.symm,GlobalC::ucell.atoms,&INPUT);
    Json::create_Json(&GlobalC::ucell, PARAM);
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
                                GlobalV::NPROC_IN_STOGROUP,
                                GlobalV::RANK_IN_STOGROUP,
                                GlobalV::MY_STOGROUP,
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
    //--------------------------------------------------
    // choose basis sets:
    // pw: plane wave basis set
    // lcao_in_pw: LCAO expaned by plane wave basis set
    // lcao: linear combination of atomic orbitals
    //--------------------------------------------------

    // where the actual stuff is done
    this->driver_run(GlobalC::ucell);

    ModuleBase::timer::finish(GlobalV::ofs_running);
    ModuleBase::Memory::print_all(GlobalV::ofs_running);
}
