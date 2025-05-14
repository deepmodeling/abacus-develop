#ifndef HAMILTLIP_H
#define HAMILTLIP_H

#include "module_hamilt_pw/hamilt_pwdft/hamilt_pw.h"
#ifdef __EXX
#include "module_ri/exx_lip.h"
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
          : HamiltPW<T, base_device::DEVICE_CPU>(pot_in, wfc_basis, p_kv, nlpp,ucell){};
#ifdef __EXX
      HamiltLIP(elecstate::Potential* pot_in,
                ModulePW::PW_Basis_K* wfc_basis,
                K_Vectors* p_kv,
                pseudopot_cell_vnl* nlpp,
                const UnitCell* ucell,
                Exx_Lip<T>& exx_lip_in)
          : HamiltPW<T, base_device::DEVICE_CPU>(pot_in, wfc_basis, p_kv, nlpp,ucell), exx_lip(exx_lip_in){};
      Exx_Lip<T>& exx_lip;
#endif
    };

} // namespace hamilt

#endif