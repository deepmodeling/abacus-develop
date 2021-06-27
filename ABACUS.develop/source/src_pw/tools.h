#ifndef TOOL_H
#define TOOL_H

using namespace std;
#include <cstdlib>
#include <new>
#include <cassert>

#include <complex>
#include <cmath>

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>

#include <vector>

#include "../src_global/constants.h"
#include "../src_global/vector3.h"
#include "../src_global/matrix.h"
#include "../module_base/matrix3.h"
#include "../src_global/realarray.h"
#include "../src_global/intarray.h"
#include "../src_global/complexarray.h"
#include "../src_global/complexmatrix.h"

#include "../src_global/memory.h"
#include "../src_global/timer.h"
#include "../src_global/lapack_connector.h"
#include "../src_global/export.h"

#include "../src_global/mathzone.h"
#include "../src_global/mathzone_add1.h"
#include "../src_global/global_function.h"
#include "../src_global/global_variable.h"
#include "../src_global/global_file.h"

#ifdef _MCD_CHECK
#include "../src_parallel/mcd.h"
#endif

#include "../src_parallel/parallel_reduce.h"
#include "../src_parallel/parallel_common.h"

#endif
