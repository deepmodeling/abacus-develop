#ifndef USPP_SUPPORT_H_
#define USPP_SUPPORT_H_

#include <string>

namespace pw
{

/**
 * Validate that a calculation using ultrasoft pseudopotentials stays within
 * the currently reviewed support boundary.
 */
void validate_uspp_support(bool use_uspp,
                           const std::string& basis_type,
                           const std::string& esolver_type,
                           int nspin,
                           int xc_func_type,
                           bool berry_phase,
                           bool towannier90,
                           bool cal_cond);

} // namespace pw

#endif
