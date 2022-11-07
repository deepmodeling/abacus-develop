#include "optical.h"
#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../src_pw/global.h"
#include "../src_parallel/parallel_reduce.h"
#include "../module_base/timer.h"

bool Optical::opt_epsilon2 = false;
int  Optical::opt_nbands = 0;

Optical::Optical(){}
Optical::~Optical(){}

