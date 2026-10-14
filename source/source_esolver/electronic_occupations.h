#ifndef ABACUS_ELECTRONIC_OCCUPATIONS_H
#define ABACUS_ELECTRONIC_OCCUPATIONS_H

namespace ModuleBase
{
class matrix;
}

namespace ModuleESolver
{
// weights already include k-point weights and spin degeneracy. Rows/bands are
// partitioned across npool groups and replicated on ranks_per_pool ranks.
// All ranks must participate; invalid occupations or a mismatch fail together.
double checked_occupation_electron_count(const ModuleBase::matrix& weights,
                                         double requested_nelect,
                                         int npool,
                                         int ranks_per_pool);
}

#endif
