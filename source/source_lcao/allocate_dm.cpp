#include "source_lcao/allocate_dm.h"
#include "source_base/basis_flattened_indexing.h"
#include "source_base/parallel_blkcyc.h"
#include "source_base/timer.h"
#include "source_cell/klist.h"
#include <vector>

namespace LCAO_domain
{

// change init_dm to allocate_dm, mohan 2025-10-31
// Moved from member function to free function so that Setup_DM itself
// can live in module_dm without depending on LCAO-layer modules.
template <typename TK>
void allocate_dm(module_dm::Setup_DM<TK>& dmat,
                 const K_Vectors* kv,
                 const BasisFlattenedIndexing& basis_index,
                 const Parallel_Blkcyclic& dist,
                 const int nspin)
{
    const int nspin_dm = nspin == 2 ? 2 : 1;
    // pass the global physical nspin so that cal_dmr can select the spin-resolved (Pauli)
    // branch for SOC/noncollinear (nspin==4), where nspin_dm itself collapses to 1.
    dmat.dm = new module_dm::DensityMatrix<TK, double>(
        basis_index, dist, nspin_dm, kv->kvec_d, kv->get_nks() / nspin_dm, nspin);
}

template void allocate_dm<double>(module_dm::Setup_DM<double>&,
                                  const K_Vectors*,
                                  const BasisFlattenedIndexing&,
                                  const Parallel_Blkcyclic&,
                                  const int);
template void allocate_dm<std::complex<double>>(module_dm::Setup_DM<std::complex<double>>&,
                                                const K_Vectors*,
                                                const BasisFlattenedIndexing&,
                                                const Parallel_Blkcyclic&,
                                                const int);

} // namespace LCAO_domain
