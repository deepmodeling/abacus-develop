#ifndef DELTASPIN_PW_H
#define DELTASPIN_PW_H

#include "source_io/module_parameter/parameter.h"

namespace pw
{

/**
 * @brief Run the inner lambda loop for DeltaSpin method to constrain atomic magnetic moments.
 *
 * This function is used in the PW basis SCF iteration to optimize lambda parameters
 * for constraining atomic magnetic moments to target values using the DeltaSpin method.
 *
 * @param iter The current iteration number (0-indexed).
 * @param drho The current charge density difference.
 * @param inp The input parameters.
 * @return true if the solver should be skipped (lambda loop was executed),
 *         false otherwise.
 */
bool run_deltaspin_lambda_loop(const int iter,
                               const double drho,
                               const Input_para& inp);

}

#endif
