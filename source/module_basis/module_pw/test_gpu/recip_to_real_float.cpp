//---------------------------------------------
// TEST for FFT
//---------------------------------------------
#include "../pw_basis.h"
#ifdef __MPI
#include "test_tool.h"
#include "module_base/parallel_global.h"
#include "mpi.h"
#endif
#include "fftw3.h"
#include "module_base/constants.h"
#include "module_base/global_function.h"
#include "pw_test.h"
#include "cuda_runtime.h"

using namespace std;
