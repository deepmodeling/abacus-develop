#pragma once
#include "lr_util.h"
namespace LR_Util
{
    template<typename T>
    void grad(const T* rhor,
        ModuleBase::Vector3<T>* gradrho,
        const ModulePW::PW_Basis& rho_basis,
        const double& tpiba)
    {
        std::vector<typename ToComplex<T>::type> rhog(rho_basis.npw);
        rho_basis.real2recip(rhor, rhog.data());
        XC_Functional::grad_rho(rhog.data(), gradrho, &rho_basis, tpiba);
    }
    template<typename T>
    void grad(const std::vector<T>& rhor,
        std::vector<ModuleBase::Vector3<T>>& gradrho,
        const ModulePW::PW_Basis& rho_basis,
        const double& tpiba)
    {
        grad(rhor.data(), gradrho.data(), rho_basis, tpiba);
    }

    template<typename T>
    void lapl(const T* rhor, T* lapn,
        const ModulePW::PW_Basis& rho_basis,
        const double& tpiba2)
    {
        ModuleBase::GlobalFunc::ZEROS(lapn, rho_basis.nrxx);
        std::vector<typename ToComplex<T>::type> rhog(rho_basis.npw);
        std::vector<T> tmp_rhor(rho_basis.nrxx);
        rho_basis.real2recip(rhor, rhog.data());
        for (int i = 0;i < 3;++i)
        {
            for (int ig = 0; ig < rho_basis.npw; ig++) { rhog[ig] *= pow(rho_basis.gcar[ig][i], 2); }
            rho_basis.recip2real(rhog.data(), tmp_rhor.data());
            for (int ir = 0; ir < rho_basis.nrxx; ir++) { lapn[ir] -= tmp_rhor[ir] * tpiba2; }
        }
    }
    template<typename T>
    void lapl(const std::vector<T>& rhor,
        std::vector<T>& lapn,
        const ModulePW::PW_Basis& rho_basis,
        const double& tpiba2)
    {
        lapl(rhor.data(), lapn.data(), rho_basis, tpiba2);
    }
}
