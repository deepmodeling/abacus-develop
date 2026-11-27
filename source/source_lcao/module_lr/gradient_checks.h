#ifndef ABACUS_LR_GRADIENT_CHECKS_H
#define ABACUS_LR_GRADIENT_CHECKS_H
#include "utils/lr_util_print.h"
namespace LR
{
// check C_uaC_va-C_uiC_vi of lumo-homo, nocc=1,  nk=1
template<typename T>
inline void test_dm_diff_H2(const T* dm, const psi::Psi<T>& c, const int nbasis)
{
    std::cout << "difference dm cal: " << std::endl;
    LR_Util::print_value(dm, nbasis, nbasis);
    std::cout << "difference dm ref: " << std::endl;
    for (int i = 0;i < nbasis;++i)
    {
        for (int j = 0;j < nbasis;++j)
        {
            std::cout << c(0, 1, i) * c(0, 1, j) - c(0, 0, i) * c(0, 0, j) << " ";
        }
        std::cout << std::endl;
    }
}

// check e_aC_uaC_va-e_iC_uiC_vi of lumo-homo, nocc=1,  nk=1
template<typename T>
inline void test_edm_H2(const T* const edm, const double* const eig_ks, const psi::Psi<T>& c, const int nbasis)
{
    std::cout << "edm cal: " << std::endl;
    LR_Util::print_value(edm, nbasis, nbasis);
    std::cout << "edm ref: " << std::endl;
    for (int i = 0;i < nbasis;++i)
    {
        for (int j = 0;j < nbasis;++j)
        {
            std::cout << eig_ks[1] * c(0, 1, i) * c(0, 1, j) - eig_ks[0] * c(0, 0, i) * c(0, 0, j) << " ";
        }
        std::cout << std::endl;
    }
}

}
#endif
