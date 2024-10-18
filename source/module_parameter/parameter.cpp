#include "parameter.h"

Parameter PARAM;

void Parameter::set_rank_nproc(const int& myrank, const int& nproc, const int& nthread_per_proc)
{
    sys.myrank = myrank;
    sys.nproc = nproc;
    sys.nthread_per_proc = nthread_per_proc;
}

void Parameter::set_start_time(const std::time_t& start_time)
{
    sys.start_time = start_time;
}