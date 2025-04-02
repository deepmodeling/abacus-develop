#include "parameter.h"

Parameter PARAM;

// changed from set_rank_nproc in 2024-1018
void Parameter::set_pal_param(const int& myrank, const int& nproc, const int& nthread_per_proc)
{
    sys.myrank = myrank;
    sys.nproc = nproc;
    sys.nthread_per_proc = nthread_per_proc;
}

void Parameter::set_start_time(const std::time_t& start_time)
{
    sys.start_time = start_time;
}
