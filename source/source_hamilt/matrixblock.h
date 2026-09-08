#ifndef MATRIXBLOCK_H
#define MATRIXBLOCK_H

#include "source_base/matrix_block.h"

namespace hamilt
{

/// MatrixBlock only describes a memory layout, so it now lives in source_base
/// and eigensolvers can use it without including the Hamiltonian interface.
/// This alias keeps the historical hamilt::MatrixBlock spelling working.
///
/// TODO: this header is a temporary compatibility shim. Once every call site
/// spells the type as ModuleBase::MatrixBlock and includes
/// source_base/matrix_block.h directly, delete this file.
using ModuleBase::MatrixBlock;

} // namespace hamilt
#endif
