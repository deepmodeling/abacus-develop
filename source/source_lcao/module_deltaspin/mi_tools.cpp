#include "mi_tools.h"

#include "source_base/tool_quit.h"
#include "source_estate/occ_comput.h"

namespace spinconstrain
{

void accumulate_Mi_from_becp(const std::complex<double>* becp,
                            int nkb,
                            int nbands,
                            int npol,
                            int spin_sign,
                            const double* wg_ik,
                            const int* nh_iat,
                            std::vector<ModuleBase::Vector3<double>>& mi)
{
    if (becp == nullptr)
    {
        ModuleBase::WARNING_QUIT("accumulate_Mi_from_becp", "becp is nullptr");
    }
    if (wg_ik == nullptr)
    {
        ModuleBase::WARNING_QUIT("accumulate_Mi_from_becp", "wg_ik is nullptr");
    }
    if (nh_iat == nullptr)
    {
        ModuleBase::WARNING_QUIT("accumulate_Mi_from_becp", "nh_iat is nullptr");
    }
    if (nkb <= 0)
    {
        ModuleBase::WARNING_QUIT("accumulate_Mi_from_becp", "nkb must be positive");
    }
    if (nbands <= 0)
    {
        ModuleBase::WARNING_QUIT("accumulate_Mi_from_becp", "nbands must be positive");
    }
    if (npol != 1 && npol != 2)
    {
        ModuleBase::WARNING_QUIT("accumulate_Mi_from_becp", "npol must be 1 or 2");
    }
    if (spin_sign != -1 && spin_sign != 1)
    {
        ModuleBase::WARNING_QUIT("accumulate_Mi_from_becp", "spin_sign must be -1 or 1");
    }

    // Compute the per-projector 2x2 occupation blocks with the shared core.
    // npol=2 -> nspin=4 (full spin density matrix); npol=1 -> nspin=2 with the
    // spin channel selected by isk (spin_sign=+1 -> isk=0, -1 -> isk=1).
    // nspin=1 is never used here: DeltaSpin only runs with nspin=2 or 4.
    const int nat = static_cast<int>(mi.size());
    const int nspin = (npol == 2) ? 4 : 2;
    const int isk = (spin_sign == 1) ? 0 : 1;
    std::vector<std::complex<double>> occ_block(nkb * 4, std::complex<double>(0.0, 0.0));
    elecstate::occ_from_proj(
        becp,
        wg_ik,
        nbands,
        npol,
        nkb,
        nspin,
        isk,
        nh_iat,
        nat,
        occ_block.data());

    // Aggregate the per-projector blocks into per-atom magnetic moments.
    // The blocks are weighted per band inside occ_from_proj, so the aggregate
    // weight is 1 here; pauli_to_moment and the z-difference read the already
    // weighted block sums.
    const double unit_weight = 1.0;
    int begin_iprj = 0;
    for (int iat = 0; iat < nat; iat++)
    {
        const int nprj = nh_iat[iat];
        if (npol == 2)
        {
            // Mi = sum_iprj pauli_to_moment(block_iprj)
            for (int iprj = 0; iprj < nprj; iprj++)
            {
                const int occ_index = (begin_iprj + iprj) * 4;
                mi[iat] += pauli_to_moment(&occ_block[occ_index], unit_weight);
            }
        }
        else
        {
            // Mz = sum_iprj (occ[0] - occ[3]) == weight * occ * spin_sign
            for (int iprj = 0; iprj < nprj; iprj++)
            {
                const int occ_index = (begin_iprj + iprj) * 4;
                mi[iat].z += (occ_block[occ_index] - occ_block[occ_index + 3]).real();
            }
        }
        begin_iprj += nprj;
    }
}

} // namespace spinconstrain
