#ifndef OTHER_PARAMETER_H
#define OTHER_PARAMETER_H
#include <string>

struct Other_para
{
    // ---------------------------------------------------------------
    // --------------       Other  Parameters         ----------------
    // ---------------------------------------------------------------
    int myrank = 0;
    int nproc = 1;
    int mypool = 0;
    int npool = 1;
    int nproc_in_pool = 1;
};
#endif