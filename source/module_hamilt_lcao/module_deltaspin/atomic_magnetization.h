#ifndef ATOMIC_MAGNETIZATION_H
#define ATOMIC_MAGNETIZATION_H

#include "module_cell/klist.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_hamilt.h"
#include "module_hamilt_lcao/hamilt_lcaodft/local_orbital_charge.h"

void calculate_MW_from_lambda(const int& step, LCAO_Hamilt& uhm, Local_Orbital_Charge& loc, const K_Vectors& kv);

#endif // ATOMIC_MAGNETIZATION_H