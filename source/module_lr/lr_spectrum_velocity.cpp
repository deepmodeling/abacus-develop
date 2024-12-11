#include "lr_spectrum.h"
#include "module_lr/dm_trans/dm_trans.h"
#include "module_hamilt_lcao/module_tddft/td_current.h"
#include "module_lr/utils/lr_util_hcontainer.h"
#include "math.h"
#include "module_parameter/parameter.h"
namespace LR
{
    /// get the velocity matrix v(R)
    inline TD_current get_velocity_matrix_R(const UnitCell& ucell,
        Grid_Driver& gd,
        const Parallel_Orbitals& pmat,
        const TwoCenterBundle& two_center_bundle)
    {
        // convert the orbital object to the old class for TD_current
        LCAO_Orbitals orb;
        const auto& inp = PARAM.inp;
        two_center_bundle.to_LCAO_Orbitals(orb, inp.lcao_ecut, inp.lcao_dk, inp.lcao_dr, inp.lcao_rmax);
        // actually this class calculates the velocity matrix v(R) at A=0
        TD_current vR(&ucell, &gd, &pmat, orb, two_center_bundle.overlap_orb.get());
        vR.calculate_vcomm_r(); // $<\mu, 0|[Vnl, r]|\nu, R>$
        vR.calculate_grad_term();   // $<\mu, 0|\nabla|\nu, R>$
        return vR;
    }

    inline double lorentz_delta(const double freq_diff, const double eta)
    {
        return eta / (freq_diff * freq_diff + eta * eta) / M_PI;
    }

    template<typename T>
    double LR::LR_Spectrum<T>::mean_square_transition_velocity(const int istate)
    {
        const TD_current& vR = get_velocity_matrix_R(ucell, gd_, pmat, two_center_bundle_);
        elecstate::DensityMatrix<T, T> DM_trans(&this->pmat, 1, this->kv.kvec_d, this->nk);
        LR_Util::initialize_DMR(DM_trans, this->pmat, this->ucell, gd_, orb_cutoff_);

        const int offset_b = istate * ldim;    //start index of band istate
        std::vector<std::complex<double>> transition_velocity(3, 0.0);    // $=\sum_{uvR} v(R) D(R) = \sum_{iak}X_{iak}<ck|v|vk>$
        double mean_transition_velocity_norm2 = 0.0;    // $= |v|^2/3$
        for (int i = 0; i < 3; i++)
        {
            for (int is = 0;is < this->nspin_x; ++is)
            {
                const int offset_x = offset_b + is * nk * pX[0].get_local_size();
                //1. transition density 
#ifdef __MPI
                std::vector<container::Tensor>  dm_trans_2d = cal_dm_trans_pblas(X + offset_x, this->pX[is], psi_ks[is], this->pc, this->naos, this->nocc[is], this->nvirt[is], this->pmat);
                // if (this->tdm_sym) for (auto& t : dm_trans_2d) LR_Util::matsym(t.data<T>(), naos, pmat);
#else
                std::vector<container::Tensor>  dm_trans_2d = cal_dm_trans_blas(X + offset_x, this->psi_ks[is], this->nocc[is], this->nvirt[is]);
                // if (this->tdm_sym) for (auto& t : dm_trans_2d) LR_Util::matsym(t.data<T>(), naos);
#endif
                for (int ik = 0;ik < this->nk;++ik) { DM_trans.set_DMK_pointer(ik, dm_trans_2d[ik].data<T>()); }
                DM_trans.cal_DMR();
                transition_velocity[i] += LR_Util::dot_R_matrix(*vR.get_current_term_pointer(i), *DM_trans.get_DMR_pointer(is + 1), ucell.nat);
            }   // end for spin_x, only matter in open-shell system
            mean_transition_velocity_norm2 += std::norm(transition_velocity[i]);
        }   // end for direction
        return mean_transition_velocity_norm2 / 3.; // for non-polarized light)
    }

    template<typename T>
    void LR::LR_Spectrum<T>::optical_absorption_velocity(const std::vector<double>& freq, const double eta, const std::string& spintype)
    {
        ModuleBase::TITLE("LR::LR_Spectrum", "optical_absorption_velocity");
        std::ofstream ofs(PARAM.globalv.global_out_dir + "absorption_" + spintype + ".dat");
        if (GlobalV::MY_RANK == 0) { ofs << "Frequency (eV) | wave length(nm) | Absorption (a.u.)" << std::endl; }
        const double fac = 4 * M_PI * M_PI * ModuleBase::e2 / ucell.omega;
        for (int f = 0;f < freq.size();++f)
        {
            double abs_value = 0.0;
            for (int i = 0;i < nstate;++i)
            {
                abs_value += mean_square_transition_velocity(i) * lorentz_delta((freq[f] - eig[i]), eta) / (eig[i] * eig[i]);
                if (GlobalV::MY_RANK == 0) { ofs << freq[f] * ModuleBase::Ry_to_eV << "\t" << 91.126664 / freq[f] << "\t" << abs_value << std::endl; }
            }
            abs_value *= fac;
        }
    }
}
template class LR::LR_Spectrum<double>;
template class LR::LR_Spectrum<std::complex<double>>;