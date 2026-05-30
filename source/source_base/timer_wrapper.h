#ifndef TIMER_WRAPPER_H
#define TIMER_WRAPPER_H

#include <chrono>

#ifdef __MPI
#include <mpi.h>
#endif

namespace ModuleBase {

/**
 * @brief Time point type that works in both MPI and non-MPI environments
 */
typedef double TimePoint;

/**
 * @brief Get current time as a TimePoint
 * 
 * @return TimePoint Current time
 */
inline TimePoint get_time()
{
#ifdef __MPI
    int is_initialized = 0;
    MPI_Initialized(&is_initialized);
    if (is_initialized)
    {
        return MPI_Wtime();
    }
    else
    {
        return std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::system_clock::now().time_since_epoch()).count() / 1e6;
    }
#else
    return std::chrono::duration_cast<std::chrono::microseconds>(
        std::chrono::system_clock::now().time_since_epoch()).count() / 1e6;
#endif
}

/**
 * @brief Calculate duration between two TimePoints in seconds
 * 
 * @param start Start time point
 * @param end End time point
 * @return double Duration in seconds
 */
inline double get_duration(const TimePoint& start, const TimePoint& end)
{
    return end - start;
}

}

#endif // TIMER_WRAPPER_H