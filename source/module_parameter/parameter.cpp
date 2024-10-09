#include "parameter.h"

Parameter PARAM;

void Parameter::set_rank_nproc(const int& myrank, const int& nproc)
{
    sys.myrank = myrank;
    sys.nproc = nproc;
}

void Parameter::set_start_time(const std::time_t& start_time)
{
    sys.start_time = start_time;
}