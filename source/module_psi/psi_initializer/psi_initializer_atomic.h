#include "module_psi/psi_initializer/psi_initializer.h"

template<typename T, typename Device>
class psi_initializer_atomic : public psi_initializer<T, Device>
{
    public:
        psi_initializer_atomic(UnitCell* p_ucell_in, 
                               ModulePW::PW_Basis_K* pw_wfc_in,
                               Parallel_Kpoints* p_parakpts_in,
                               Representation<T, Device>* p_rep_in);
        ~psi_initializer_atomic();
        void initialize(const psi::Psi<T, Device>& psi_out);
};