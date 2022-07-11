#ifndef __METAPW
#define __METAPW

#include "operator.h"
#include "module_base/matrix.h"
#include "module_pw/pw_basis_k.h"

namespace hamilt
{

class MetaPW : public Operator
{
    public:
    MetaPW(
        int max_npw_in,
        int npol_in,
        double tpiba2_in,
        const int* ngk_in, 
        const int* isk_in,
        const ModuleBase::matrix* vk,
        ModulePW::PW_Basis_K* wfcpw
    );

    void act(const std::complex<double> *psi_in, std::complex<double> *hpsi, const size_t size) const override;

    private:

    int max_npw = 0;

    int npol = 0;

    double tpiba2 = 0.0;

    const int* ngk = nullptr;

    const int* isk = nullptr;

    const ModuleBase::matrix* vk = nullptr;

    ModulePW::PW_Basis_K* wfcpw = nullptr;
};

} // namespace hamilt

#endif