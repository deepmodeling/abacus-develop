#include "density_matrix.h"

#include "source_base/libm/libm.h"
#include "source_base/tool_title.h"
#include "source_base/tool_quit.h"
#include "source_base/constants.h"
#include "source_base/timer.h"
#include "source_cell/klist.h"

namespace module_dm
{

// shared inner loop of cal_dmr / cal_dmr_td: accumulate kphase * DMK into DMR blocks
template <typename TK, typename TR_in, typename TR_out>
void accumulate_dmr(
    DensityMatrix<TK, TR_in>& dm,
    std::vector<hamilt::HContainer<TR_out>*>& dmR_out,
    const std::map<ModuleBase::Vector3<int>, std::complex<double>>& phase_hybrid,
    const int ik_in,
    const char* func_name)
{
    const int ld_hk = dm.pv->nrow;
    for (int is = 1; is <= dm.spin_mult; ++is)
    {
        const int ik_begin = dm._nk * (is - 1); // jump dm._nk for spin_down if nspin==2
        hamilt::HContainer<TR_out>* const dmr_spin = dmR_out[is - 1];
        // set zero since this function is called in every scf step
        dmr_spin->set_zero();

        if (dm.nspin != 1 && dm.nspin != 2 && dm.nspin != 4)
        {
            ModuleBase::WARNING_QUIT(func_name, "nspin must be 1, 2 or 4");
        }

        // accumulate kphase * DMK into DMR blocks; for nspin=4 (SOC) each orbital
        // corresponds to a 2x2 spin block, so rows/cols step by 2 (physical spin
        // dimension), and each block is transformed to Pauli components
        // (rho_0, rho_x, rho_y, rho_z)
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
        for (int i = 0; i < dmr_spin->size_atom_pairs(); ++i)
        {
            hamilt::AtomPair<TR_out>& atom_pair = dmr_spin->get_atom_pair(i);
            const DmrBlock block = get_dmr_block(dm.pv, atom_pair.get_atom_i(), atom_pair.get_atom_j());

            // precompute k-phase factors and collect DMR block pointers
            std::vector<std::vector<TK>> kphase_vec;
            std::vector<TR_out*> dmr_mats;
            build_kphase(atom_pair, dm._kvec_d, dm._nk, phase_hybrid, kphase_vec, dmr_mats);

            if (dm.nspin == 1 || dm.nspin == 2)
            {
                // nspin=1/2: DMR_ij(R) += Re[ e^{ik·R} * DMK_ij(k) ]
                // (sum over ik when ik_in < 0, single ik when ik_in >= 0)
                add_dmr_real(dm, block, ik_begin, kphase_vec, ld_hk, ik_in, dmr_mats);
            }
            else if (dm.nspin == 4)
            {
                // nspin==4 (SOC): first accumulate S_ij(R) = sum_k e^{ik·R} * DMK_ij(k),
                // then for each 2x2 spin block (upup, updown, downup, downdown) transform to
                // Pauli components via xyz_to_updown:
                //   rho_0 = rho_upup + rho_downdown
                //   rho_x = rho_updown + rho_downup
                //   rho_y = Im(rho_updown) - Im(rho_downup)   (sign for conjugated stored DM)
                //   rho_z = rho_upup - rho_downdown
                add_dmr_soc(dm, block, ik_begin, kphase_vec, ld_hk, ik_in,
                            atom_pair.get_col_size(), dmr_mats);
            }
            else
            {
                ModuleBase::WARNING_QUIT(func_name, "nspin must be 1, 2 or 4");
            }
        }
    }
}

// calculate DMR from DMK using blas for multi-k calculation
template <typename TK, typename TR_in, typename TR_out>
void cal_dmr(
    DensityMatrix<TK, TR_in>& dm,
    std::vector<hamilt::HContainer<TR_out>*>& dmR_out,
    const int ik_in)
{
    ModuleBase::TITLE("DensityMatrix", "cal_dmr");
    ModuleBase::timer::start("DensityMatrix", "cal_dmr");

    // To check whether DMR has been initialized
    if (dmR_out.size() != dm.spin_mult)
    {
        ModuleBase::WARNING_QUIT("module_dm::cal_dmr",
                                 "DMR has not been initialized: dmR_out.size() != spin_mult!");
    }

    // validate ik_in: either -1 (all k-points) or a valid index
    if (ik_in < -1 || ik_in >= dm._nk)
    {
        ModuleBase::WARNING_QUIT("module_dm::cal_dmr",
                                 "ik_in out of range: must be -1 (all k) or 0 <= ik_in < nk");
    }

    const std::map<ModuleBase::Vector3<int>, std::complex<double>> no_hybrid_phase;
    accumulate_dmr(dm, dmR_out, no_hybrid_phase, ik_in, "module_dm::cal_dmr");
    dm._dmr_ready = true;
    ModuleBase::timer::end("DensityMatrix", "cal_dmr");
}

template <>
void DensityMatrix<std::complex<double>, double>::cal_dmr(const int ik_in)
{
    module_dm::cal_dmr(*this, this->dmr, ik_in);
}

template <>
void DensityMatrix<std::complex<double>, std::complex<double>>::cal_dmr(const int ik_in)
{
    module_dm::cal_dmr(*this, this->dmr, ik_in);
}

// explicit instantiations for accumulate_dmr (used by both cal_dmr here and
// cal_dmr_td in dmr_td.cpp; without these the TD instantiations are missing
// at link time)
template void accumulate_dmr<std::complex<double>, double, double>(
    DensityMatrix<std::complex<double>, double>&,
    std::vector<hamilt::HContainer<double>*>&,
    const std::map<ModuleBase::Vector3<int>, std::complex<double>>&,
    const int,
    const char*);

template void accumulate_dmr<std::complex<double>, std::complex<double>, std::complex<double>>(
    DensityMatrix<std::complex<double>, std::complex<double>>&,
    std::vector<hamilt::HContainer<std::complex<double>>*>&,
    const std::map<ModuleBase::Vector3<int>, std::complex<double>>&,
    const int,
    const char*);

} // namespace module_dm
