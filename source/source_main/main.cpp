//==========================================================
// AUTHOR : mohan
// DATE : 2008-11-10
//==========================================================

#include "source_main/driver.h"
#include "fftw3.h"
#include "source_base/parallel_global.h"
#include "source_io/parse_args.h"
#include "source_io/module_parameter/parameter.h"
#ifdef _OPENMP
#include <omp.h>
#endif
// #include "source_lcao/module_operator_lcao/fssh_driver.h"

int main(int argc, char** argv)
{
    // ========================================================
    // FSSH 测试截断点 (测试完记得删掉这段或注释掉)
    // ========================================================
    // FsshDriver test_driver;
    // test_driver.test_math_1(); // 运行测试
    // test_driver.test_math_2(); // 运行测试
    // test_driver.test_math_3(); // 运行测试
    // test_driver.test_math_4(); // 运行测试
    // test_driver.test_math_5(); // 运行测试
    // test_driver.test_tully_model_1();
    // test_driver.test_tully_model_1_scan();
    // test_driver.test_tully_model_2();
    // test_driver.test_tully_model_2_scan();
    // test_driver.test_tully_model_3();
    // test_driver.test_tully_model_3_scan();
    // return 0;                // 直接退出程序，不跑后面的庞大DFT代码了
    // ========================================================
    /*
    read the arguement in the command-line,
    with "abacus -v", the program exit and returns version info,
    with no arguments, the program continues.
    */
    ModuleIO::parse_args(argc, argv);

    /*
    read the mpi parameters in the command-line,
    initialize the mpi environment.
    */
    int nproc = 1;
    int my_rank = 0;
    int nthread_per_proc = 1;
    Parallel_Global::read_pal_param(argc, argv, nproc, nthread_per_proc, my_rank);
#ifdef _OPENMP
    // ref: https://www.fftw.org/fftw3_doc/Usage-of-Multi_002dthreaded-FFTW.html
    fftw_init_threads();
    fftw_plan_with_nthreads(omp_get_max_threads());
#endif
    PARAM.set_pal_param(my_rank, nproc, nthread_per_proc);

    /*
    main program for doing electronic structure calculations.
    */
    Driver DD;
    DD.init();
    // FsshDriver lowdin_validator;
    // lowdin_validator.run_lowdin_validation_suite_8_5("/home/shane/ABACUS/test/lowdin_validation_4/lowdin_validation_8_5");
    /*
    After running mpi version of abacus, release the mpi resources.
    */
#ifdef __MPI
    Parallel_Global::finalize_mpi();
#endif
#ifdef _OPENMP
    fftw_cleanup_threads();
#endif

    return 0;
}
