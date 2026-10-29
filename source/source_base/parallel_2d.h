#ifndef _PARALLEL_2D_H_
#define _PARALLEL_2D_H_

#include "block_cyclic_distribution_information.h"

// DEPRECATED-ALIAS: Parallel_2D -> BlockCyclicDistributionInformation
// DELETION-COMMIT: remove in S5; trigger: Parallel_2D and parallel_2d.h include burn-down == 0.
using Parallel_2D = BlockCyclicDistributionInformation;

#endif
