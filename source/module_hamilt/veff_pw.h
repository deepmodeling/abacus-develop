#ifndef __VEFFPW
#define __VEFFPW

#include "operator.h"
#include "module_base/matrix.h"
#include "module_pw/pw_basis_k.h"

namespace hamilt
{

class VeffPW : public Operator
{
    public:
    VeffPW(
        int max_npw_in,
        int npol_in,
        const int* ngk_in,
        const int* isk_in,
        const ModuleBase::matrix* veff_in,
        ModulePW::PW_Basis_K* wfcpw_in
    );

    void act(const std::complex<double> *psi_in, std::complex<double> *hpsi, const size_t size) const override;

    private:

    int max_npw = 0;

    int npol = 0;

    const int* ngk = nullptr;

    const int* isk = nullptr;

    const ModuleBase::matrix* veff = nullptr;

    ModulePW::PW_Basis_K* wfcpw = nullptr;
};

} // namespace hamilt

#endif