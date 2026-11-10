#ifndef HSOLVER_ORB_BUFFER_H
#define HSOLVER_ORB_BUFFER_H

#include "source_base/module_container/ATen/core/tensor.h"
#include "source_base/module_device/memory_op.h"

namespace hsolver
{
/** @brief Column-major orbitals; spin components have ld/npol rows each. */
struct OrbLayout
{
    int ld;
    int npw;
    int npol;
    int bands;
};

/** @brief Clear all spin-component tails, including scalar-orbital padding. */
template <typename T, typename Device>
void zero_pad(T* data, const OrbLayout& layout);

/** @brief Borrow scalar input, or copy spinors into scratch and clear their padding.
 *  Scratch must have room for ld * bands coefficients; input remains unchanged.
 */
template <typename T, typename Device>
const T* clean_orbs(const T* input, const OrbLayout& layout, ct::Tensor* scratch);
} // namespace hsolver
#endif
