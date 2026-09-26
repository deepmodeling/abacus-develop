#ifndef HAMILTLIP_H
#define HAMILTLIP_H

#include "source_hamilt/hamilt_hs_adapter.h"
#include "source_pw/module_pwdft/hamilt_pw.h"
#ifdef __EXX
#include "source_lcao/module_ri/exx_lip.h"
#endif

namespace hamilt
{

    template <typename T>
    class HamiltLIP : public HamiltPW<T, base_device::DEVICE_CPU>
    {
    public:
      HamiltLIP(elecstate::Potential* pot_in,
                ModulePW::PW_Basis_K* wfc_basis,
                K_Vectors* p_kv,
                pseudopot_cell_vnl* nlpp,
                const UnitCell* ucell)
          : HamiltPW<T, base_device::DEVICE_CPU>(pot_in, wfc_basis, p_kv, nlpp, nullptr, ucell, nullptr, false){};
#ifdef __EXX
      HamiltLIP(elecstate::Potential* pot_in,
                ModulePW::PW_Basis_K* wfc_basis,
                K_Vectors* p_kv,
                pseudopot_cell_vnl* nlpp,
                const UnitCell* ucell,
                Exx_Lip<T>& exx_lip_in)
          : HamiltPW<T, base_device::DEVICE_CPU>(pot_in, wfc_basis, p_kv, nlpp, nullptr, ucell, nullptr, false),
            exx_lip(exx_lip_in){};
      Exx_Lip<T>& exx_lip;
#endif
    };

    /// HamiltLIP seen through hsolver::HSOperator. Besides H and S it feeds the
    /// EXX term into the subspace Hamiltonian and hands the subspace
    /// eigenvectors back to Exx_Lip, which is what HSolverLIP needs.
    template <typename T>
    class HamiltLIPHSOperator : public HamiltHSOperator<T, base_device::DEVICE_CPU>
    {
      public:
        HamiltLIPHSOperator(HamiltLIP<T>* hm,
                            const ModulePW::PW_Basis_K* wfc_basis,
                            const bool cal_exx,
                            const double hybrid_alpha)
            : HamiltHSOperator<T, base_device::DEVICE_CPU>(hm, wfc_basis), hm_lip_(hm), cal_exx_(cal_exx),
              hybrid_alpha_(hybrid_alpha){};

#ifdef __EXX
        void add_to_subspace_h(T* hcc, const int naos) const override
        {
            if (!cal_exx_)
            {
                return;
            }
            const int ik = this->ik_;
            for (int n = 0; n < naos; ++n)
            {
                for (int m = 0; m < naos; ++m)
                {
                    hcc[n * naos + m] += (T)hybrid_alpha_ * hm_lip_->exx_lip.get_exx_matrix()[ik][m][n];
                }
            }
        }

        void export_subspace_vec(const T* vcc, const int naos, const int nbands) const override
        {
            if (cal_exx_)
            {
                hm_lip_->exx_lip.set_hvec(this->ik_, vcc, naos, nbands);
            }
        }
#endif

      private:
        HamiltLIP<T>* hm_lip_ = nullptr;
        const bool cal_exx_;
        const double hybrid_alpha_;
    };

} // namespace hamilt

#endif
