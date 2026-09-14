/**
 * @file check_atomic_stru.h
 * @brief Function for checking atomic structure.
 */
#ifndef CHECK_ATOMIC_STRU_H
#define CHECK_ATOMIC_STRU_H

#include "unitcell.h"

namespace unitcell
{
    /**
     * @brief Check atomic structure.
     *
     * @param ucell unit cell [in/out]
     * @param factor scaling factor [in]
     */
    void check_atomic_stru(UnitCell& ucell, const double& factor);
};

#endif
