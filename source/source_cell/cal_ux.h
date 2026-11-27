/**
 * @file cal_ux.h
 * @brief Functions for calculating ux and related operations.
 */
#ifndef CAL_UX_H
#define CAL_UX_H

#include "source_cell/unitcell.h"

namespace unitcell {

/**
 * @brief Calculate ux for the unit cell.
 *
 * @param ucell unit cell [in/out]
 * @param nspin number of spin components [in]
 */
void cal_ux(UnitCell& ucell, const int nspin);

/**
 * @brief Judge if two vectors are parallel.
 *
 * @param a first vector [in]
 * @param b second vector [in]
 * @return true if vectors are parallel
 */
bool judge_parallel(double a[3], ModuleBase::Vector3<double> b);

}

#endif