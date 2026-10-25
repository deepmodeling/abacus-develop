#include "lr_force.h"
#include "source_lcao/pulay_fs.h"
#include "source_lcao/module_lr/utils/lr_util_hcontainer.h"
#include "source_io/module_output/output_log.h"
#ifdef __EXX
#include "operator_casida/operator_lr_exx.h" // gs_is_hybrid
#endif
namespace LR
{
    template<typename TK>
    ModuleBase::matrix LR_Force<TK>::reproduce_force_gs(const K_Vectors& kv,
        const module_dm::DensityMatrix<TK, double>& dm_gs,
        const module_dm::DensityMatrix<TK, double>& edm_gs)
    {
        // local + Hartree + xc term, including Hellmann-Feynman and Pulay
        ModuleBase::matrix f_gs_hf_pulay = cal_force_hamilt_gs_dm_relaxed_diff(dm_gs, dm_gs, true); // pw(vl_dvl+ewald)+vnl+t_dphi+vl_dphi
        // edm term
        ModuleBase::matrix f_nonortho = cal_force_overlap_edm(edm_gs); // overlap
#ifdef __EXX
        if (gs_is_hybrid(this->dft_functional_))
        {
            const auto& Ds_gs = LR_Util::get_exx_Ds_gs(dm_gs, ucell_, kv, pv_);
            const auto& Ds_gs_2 = LR_Util::get_exx_Ds_gs(dm_gs, ucell_, kv, pv_);
            // at nspin=1 there is only one channel, and `get_exx_Ds_gs` already returns 0.5*D
            // (= D_up), so both halves below reuse channel 0.
            const int is_2nd = (Ds_gs.size() > 1) ? 1 : 0;
            ModuleBase::matrix f_gs_exx(ucell_.nat, 3);
            // test the two function using the two spin channels respectively
            // 0.5 is from dE = 0.5 dTr[D(HD)]. No 0.5 in excited-state calculateion of dTr[(T+Z)(HD)]
            f_gs_exx += cal_force_exx_gs_dm_relaxed_diff(Ds_gs.at(0), Ds_gs_2.at(0), alpha_, std::to_string(0)) * 0.5;    // test passed， = 0.5 groud-state EXX force
            f_gs_exx += cal_force_exx_dm_trans(Ds_gs.at(is_2nd), alpha_, std::to_string(is_2nd)) * 0.5;
            if (this->test_force_)
                ModuleIO::print_force(this->ofs_running_, ucell_, "EXX GS FORCE reproduce (eV/Angstrom)", f_gs_exx, false);
            f_gs_hf_pulay += f_gs_exx;
        }
#endif
        return f_gs_hf_pulay + f_nonortho;
    }

    template<typename TK>
    ModuleBase::matrix LR_Force<TK>::reproduce_force_gs_loc(
        const module_dm::DensityMatrix<TK, double>& dm_gs,
        const elecstate::Potential& pot_gs)
    {
        Charge chr_gs;
        this->dm_to_charge(dm_gs, chr_gs);
        //  local pp (Pulay) + Hartree + xc (grid integration)
        ModuleBase::matrix fvl_dphi(this->ucell_.nat, 3);
        ModuleBase::matrix stress_tmp;  // no use now, only for passing into interfaces
        PulayForceStress::cal_pulay_fs(dm_gs.get_dmr_vec().size()/*nspin*/, fvl_dphi, stress_tmp,
            dm_gs, this->ucell_, &pot_gs, true, false);
        Parallel_Reduce::reduce_pool(fvl_dphi.c, fvl_dphi.nr * fvl_dphi.nc);   // see lr_force.cpp's `fvl_dphi`
        return fvl_dphi;
    }

}

template class LR::LR_Force<double>;
template class LR::LR_Force<std::complex<double>>;
