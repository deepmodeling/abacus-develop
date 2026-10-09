#include "density_matrix.h"

#include <numeric>

#include "source_base/libm/libm.h"
#include "source_base/tool_title.h"
#include "source_base/tool_quit.h"
#include "source_base/constants.h"
#include "source_base/timer.h"
#include "source_cell/klist.h"

namespace module_dm
{

template <typename TK, typename TR_in, typename TR_out>
void cal_dmr_full(
    const DensityMatrix<TK, TR_in>& dm,
    hamilt::HContainer<TR_out>* dmR_out,
    const int ik_in)
{
    ModuleBase::TITLE("DensityMatrix", "cal_dmr_full");

    // validate ik_in: either -1 (all k-points) or a valid index
    if (ik_in < -1 || ik_in >= dm._nk)
    {
        ModuleBase::WARNING_QUIT("module_dm::cal_dmr_full",
                                 "ik_in out of range: must be -1 (all k) or 0 <= ik_in < nk");
    }

    ModuleBase::timer::start("DensityMatrix", "cal_dmr_full");
    const int ld_hk = dm.pv->nrow;
    hamilt::HContainer<TR_out>* const dmr_full = dmR_out;
    dmr_full->set_zero();
    const std::map<ModuleBase::Vector3<int>, std::complex<double>> no_hybrid_phase;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (int i = 0; i < dmr_full->size_atom_pairs(); ++i)
    {
        hamilt::AtomPair<TR_out>& atom_pair = dmr_full->get_atom_pair(i);
        const DmrBlock block = get_dmr_block(dm.pv, atom_pair.get_atom_i(), atom_pair.get_atom_j());
        const int R_size = atom_pair.get_R_size();

        // precompute k-phase factors and collect DMR block pointers
        std::vector<std::vector<TK>> kphase_vec;
        std::vector<TR_out*> dmr_mats;
        build_kphase(atom_pair, dm._kvec_d, dm._nk, no_hybrid_phase, kphase_vec, dmr_mats);

        // transpose DMK block to row-major, then axpy into each R-vector block
        // DMR_ij(R) += e^{ik·R} * DMK_ij(k)
        // (sum over ik when ik_in < 0, single ik when ik_in >= 0)
        std::vector<TK> dmk_row(block.size());
        if (ik_in >= 0)
        {
            // single k-point
            const TK* const dmk_col = dm.dmk[ik_in].data() + block.col0 * ld_hk + block.row0;
            transpose_dmk_block(dmk_col, ld_hk, block, dmk_row.data());
            for (int iR = 0; iR < R_size; ++iR)
            {
                BlasConnector::axpy(block.size(),
                                    kphase_vec[ik_in][iR],
                                    dmk_row.data(),
                                    1,
                                    dmr_mats[iR],
                                    1);
            }
        }
        else
        {
            // all k-points
            for (int ik = 0; ik < dm._nk; ++ik)
            {
                const TK* const dmk_col = dm.dmk[ik].data() + block.col0 * ld_hk + block.row0;
                transpose_dmk_block(dmk_col, ld_hk, block, dmk_row.data());
                for (int iR = 0; iR < R_size; ++iR)
                {
                    BlasConnector::axpy(block.size(),
                                        kphase_vec[ik][iR],
                                        dmk_row.data(),
                                        1,
                                        dmr_mats[iR],
                                        1);
                }
            }
        }
    }
    ModuleBase::timer::end("DensityMatrix", "cal_dmr_full");
}

template <>
void DensityMatrix<double, double>::cal_dmr_full(
    hamilt::HContainer<std::complex<double>>* dmR_out,
    const int ik_in) const
{
}
template <>
void DensityMatrix<std::complex<double>, double>::cal_dmr_full(
    hamilt::HContainer<std::complex<double>>* dmR_out,
    const int ik_in) const
{
    module_dm::cal_dmr_full(*this, dmR_out, ik_in);
}

} // namespace module_dm
