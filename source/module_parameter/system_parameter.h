#ifndef SYSTEM_PARAMETER_H
#define SYSTEM_PARAMETER_H
#include <ctime>
#include <string>

struct System_para
{
    // ---------------------------------------------------------------
    // --------------       Other  Parameters         ----------------
    // ---------------------------------------------------------------
    int myrank = 0;
    int nproc = 1;
    int mypool = 0;
    int npool = 1;
    int nproc_in_pool = 1;
    std::time_t start_time = 0;
};
#endif