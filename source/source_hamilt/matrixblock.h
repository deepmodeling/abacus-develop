#ifndef MATRIXBLOCK_H
#define MATRIXBLOCK_H

#include "source_base/matrix_block.h"

namespace hamilt
{

/// MatrixBlock only describes a memory layout, so it now lives in source_base
/// and eigensolvers can use it without including the Hamiltonian interface.
/// This alias keeps the historical hamilt::MatrixBlock spelling working.
using ModuleBase::MatrixBlock;

} // namespace hamilt
#endif
