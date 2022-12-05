#ifndef TOOL_THREADING_H 
#define TOOL_THREADING_H 

namespace ModuleBase
{

//==========================================================
// NAME : TASK_DIST_1D
// Distributing 1d tasks by worker id
//==========================================================
template <typename T_task, typename T_out>
static inline void TASK_DIST_1D(int nworker, int iworker, T_task ntask, T_out& start, T_out& len)
{
    if (nworker == 1)
    {
        start = 0;
        len = ntask;
    }
    else
    {
        const T_task tlen = ntask / nworker;
        const T_task trem = ntask - tlen * nworker;
        if (iworker < trem)
        {
            start = tlen * iworker + iworker;
            len = tlen + 1;
        }
        else
        {
            start = tlen * iworker + trem;
            len = tlen;
        }
    }
}

}

#endif
