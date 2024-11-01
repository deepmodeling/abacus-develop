#include "lr_util.h"
namespace LR_Util
{
    void grad(const double* rhor,
        ModuleBase::Vector3<double>* gdr,
        const ModulePW::PW_Basis& rho_basis,
        const double& tpiba)
    {
        std::vector<std::complex<double>> rhog(rho_basis.npw);
        rho_basis.real2recip(rhor, rhog.data());
        XC_Functional::grad_rho(rhog.data(), gdr, &rho_basis, tpiba);
    }
    void laplace(const double* rhor, double* lapn,
        const ModulePW::PW_Basis& rho_basis,
        const double& tpiba2)
    {
        ModuleBase::GlobalFunc::ZEROS(lapn, rho_basis.nrxx);
        std::vector<std::complex<double>> rhog(rho_basis.npw);
        std::vector<double> tmp_rhor(rho_basis.nrxx);
        rho_basis.real2recip(rhor, rhog.data());
        for (int i = 0;i < 3;++i)
        {
            for (int ig = 0; ig < rho_basis.npw; ig++) { rhog[ig] *= pow(rho_basis.gcar[ig][i], 2); }
            rho_basis.recip2real(rhog.data(), tmp_rhor.data());
            for (int ir = 0; ir < rho_basis.nrxx; ir++) { lapn[ir] -= tmp_rhor[ir] * tpiba2; }
        }
    }
}
