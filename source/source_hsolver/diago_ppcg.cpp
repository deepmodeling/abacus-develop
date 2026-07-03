#include "diago_ppcg.h"

#include "ppcg/diago_ppcg_reduce.hpp"
#include "ppcg/diago_ppcg_small_eigen.hpp"
#include "ppcg/diago_ppcg_ops.hpp"
#include "ppcg/diago_ppcg_subspace.hpp"
#include "ppcg/diago_ppcg_orth.hpp"
#include "ppcg/diago_ppcg_cg.hpp"
#include "ppcg/diago_ppcg_diag.hpp"

namespace hsolver {

// =============================================================================
// Explicit template instantiation (CPU only; extend for GPU as needed)
// =============================================================================
template class DiagoPPCG<std::complex<float>,  base_device::DEVICE_CPU>;
template class DiagoPPCG<std::complex<double>, base_device::DEVICE_CPU>;

} // namespace hsolver
