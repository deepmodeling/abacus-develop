#ifndef TOOL_THREADING_H 
#define TOOL_THREADING_H 

#ifdef _OPENMP
#include <omp.h>
#endif

#include <functional>

namespace ModuleBase
{

//==========================================================
// NAME : TASK_DIST_1D
// Distributing 1d tasks by worker id
//==========================================================
template <typename T_task, typename T_out>
inline void TASK_DIST_1D(int nworker, int iworker, T_task ntask, T_out& start, T_out& len)
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

template <typename T_task, typename T_out>
inline void BLOCK_TASK_DIST_1D(int nworker, int iworker, T_task ntask, T_task block_size, T_out& start, T_out& len)
{
    if (nworker == 1)
    {
        start = 0;
        len = ntask;
    }
    else
    {
        const T_task nblock = (ntask + block_size - 1) / block_size;
        const T_task blen = nblock / nworker;
        const T_task brem = nblock - blen * nworker;
        if (iworker < brem)
        {
            start = blen * iworker + iworker;
            len = blen + 1;
        }
        else
        {
            start = blen * iworker + brem;
            len = blen;
        }
        start *= block_size;
        len *= block_size;
        if (start > ntask)
        {
            start = ntask;
            len = 0;
            return;
        }
        if (start + len > ntask)
        {
            len = ntask - start;
        }
    }
}

inline void OMP_PARALLEL(const std::function<void(int, int)> &f)
{
#ifdef _OPENMP
    #pragma omp parallel
    {
        f(omp_get_num_threads(), omp_get_thread_num());
    }
#else
    f(1, 0);
#endif
}

inline void TRY_OMP_PARALLEL(const std::function<void(int, int)> &f)
{
#ifdef _OPENMP
    if (!omp_in_parallel())
    {
        OMP_PARALLEL(f);
    }
    else
#endif
    {
        f(1, 0);
    }
}

}

#endif
